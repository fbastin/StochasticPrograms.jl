# MIT License
#
# Copyright (c) 2018 Martin Biel
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

abstract type RandomVariable end
abstract type RandomVector end

struct RawStoch{N}
    name::String
    random_variables::NTuple{N, Vector{RandomVariable}}
    random_vectors::NTuple{N, Vector{RandomVector}}

    function RawStoch(name::String, ran_vars::NTuple{N, Vector{RandomVariable}}, ran_vectors::NTuple{N, Vector{RandomVector}}) where N
        return new{N}(name, ran_vars, ran_vectors)
    end
end

struct IndepDiscrete{T <: AbstractFloat} <: RandomVariable
    rowcol::RowCol
    inclusion::InclusionType
    support::Vector{T}
    probabilities::Vector{T}

    function IndepDiscrete(::Type{T}, rowcol::RowCol, inclusion::InclusionType) where T <: AbstractFloat
        inclusion in INCLUSIONS || error("Unknown inclusion $inclusion")
        return new{T}(rowcol, inclusion, Vector{T}(), Vector{T}())
    end
end

struct IndepDistribution{T <: AbstractFloat} <: RandomVariable
    rowcol::RowCol
    inclusion::InclusionType
    distribution::DistributionType
    parameters::Pair{T}

    function IndepDistribution(rowcol::RowCol, inclusion::InclusionType, distribution::DistributionType, parameters::Pair{T}) where T <: AbstractFloat
        inclusion in INCLUSIONS || error("Unknown inclusion $inclusion")
        distribution in DISTRIBUTIONS || error("Unknown distribution $distribution")
        return new{T}(rowcol, inclusion, distribution, parameters)
    end
end

struct BlockDiscrete{T <: AbstractFloat} <: RandomVector
    inclusion::InclusionType
    support::Vector{Dict{RowCol,T}}
    probabilities::Vector{T}

    function BlockDiscrete(::Type{T}, inclusion::InclusionType) where T
        inclusion in INCLUSIONS || error("Unknown inclusion $inclusion")
        return new{T}(inclusion, Vector{Dict{RowCol,T}}(), Vector{T}())
    end
end

function parse_sto(::Type{T}, tim::RawTime, cor::RawCor, filename::AbstractString; N::Integer = 2) where T <: AbstractFloat
    # Initialize auxiliary variables
    name      = "SLP"
    mode      = INDEP
    dist      = DISCRETE
    inclusion = REPLACE
    N         = num_stages(tim)
    stage     = 2
    v₁        = zero(T)
    v₂        = zero(T)
    block     = :none
    # Prepare sto data
    random_vars = ntuple(Val(N-1)) do i
        Vector{RandomVariable}()
    end
    random_vecs = ntuple(Val(N-1)) do i
        Vector{RandomVector}()
    end
    # Bookkeeping
    ranvar_indices = ntuple(Val(N-1)) do i
        Dict{RowCol, Int}()
    end
    ranvec_indices = ntuple(Val(N-1)) do i
        Dict{Block, Int}()
    end
    root_scenarios = Dict{Symbol, Int}()   # the scenarios branching from the root => their outcome
    # Parse the file
    data = open(filename) do io
        firstline = split(first_data_line(io))
        if !isempty(firstline) && Symbol(firstline[1]) == :STOCH
            name = join(firstline[2:end], " ")
        else
            throw(ArgumentError("`STOCH` field is expected on the first line."))
        end
        for line in eachline(io)
            if mode == END
                # Parse finished
                break
            end
            is_data(line) || continue
            words = split(line)
            first_word = Symbol(words[1])
            if first_word in STO_MODES
                mode = first_word
                if length(words) == 2
                    dist = DistributionType(words[2])
                    dist in DISTRIBUTIONS || error("Unknown distribution $dist.")
                elseif length(words) == 3
                    dist = DistributionType(words[2])
                    dist in DISTRIBUTIONS || error("Unknown distribution $dist.")
                    inclusion = InclusionType(words[3])
                    inclusion in INCLUSIONS || error("Unknown inclusion $inclusion.")
                end
                 continue
            end
            if mode == INDEP
                length(words) == 4 || length(words) == 5 || error("Malformed sto file at line $words.")
                col = Col(words[1])
                col = col == cor.rhsname ? RHS : col
                row = Row(words[2])
                if length(words) == 4
                    v₁ = parse(Float64, words[3])
                    v₂ = parse(Float64, words[4])
                elseif length(words) == 5
                    v₁ = parse(Float64, words[3])
                    v₂ = parse(Float64, words[5])
                    stage_name = Symbol(words[4])
                    haskey(tim.stages, stage_name) || error("Stage name $stage_name not specified in .tim file.")
                    stage = tim.stages[stage_name]
                else
                    error("Malformed sto file at line $words")
                end
                if dist == DISCRETE
                    # Check if random variable exists already
                    idx = get(ranvar_indices[stage - 1], (row, col), 0)
                    ran_var = if idx > 0
                        ran_var = random_vars[stage - 1][idx]
                        # Sanity check
                        ran_var isa IndepDiscrete || error("Random variable at $((row, col)) has more than one specified distribution.")
                        ran_var
                    else
                        # Create new INDEP random variable
                        ran_var = IndepDiscrete(T, (row,col), inclusion)
                        push!(random_vars[stage - 1], ran_var)
                        # Bookkeep
                        ranvar_indices[stage - 1][(row, col)] = length(random_vars[stage - 1])
                        ran_var
                    end
                    # Update values
                    push!(ran_var.support, v₁)
                    push!(ran_var.probabilities, v₂)
                else
                    haskey(ranvar_indices[stage - 1], (row, col)) && error("Random variable at $((row, col)) has more than one specified distribution.")
                    push!(random_vars[stage - 1], IndepDistribution((row, col), inclusion, dist, Pair(v₁, v₂)))
                    # Bookkeep
                    ranvar_indices[(row, col)] = length(random_vars[stage - 1])
                end
            elseif mode == BLOCKS
                if dist == DISCRETE
                    if words[1] == "BL"
                        length(words) == 4 || error("Malformed sto file at line $words.")
                        stage_name = Symbol(words[3])
                        haskey(tim.stages, stage_name) || error("Stage name $stage_name not specified in .tim file.")
                        stage = tim.stages[stage_name]
                        prob = parse(Float64, words[4])
                        # The block of this outcome, known by its name
                        block = Block(words[2])
                        idx = get(ranvec_indices[stage - 1], block, 0)
                        ran_vec = if idx > 0
                            ran_vec = random_vecs[stage - 1][idx]
                            # Sanity check
                            ran_vec isa BlockDiscrete || error("Random vector at $block has more than one specified distribution.")
                            ran_vec
                        else
                            # Create new DISCRETE random vector
                            ran_vec = BlockDiscrete(T, inclusion)
                            push!(random_vecs[stage - 1], ran_vec)
                            # Bookkeep
                            ranvec_indices[stage - 1][block] = length(random_vecs[stage - 1])
                            ran_vec
                        end
                        # Update probability
                        push!(ran_vec.support, Dict{RowCol,T}())
                        push!(ran_vec.probabilities, prob)
                        continue
                    else
                        read_entries!(random_vecs[stage - 1][ranvec_indices[stage - 1][block]], words, cor)
                    end
                else
                    error("Distribution $dist is not supported.")
                end
            elseif mode == SCENARIOS
                if words[1] == "SC"
                    length(words) == 5 || error("Malformed+ sto file at line $words.")
                    stage_name = Symbol(words[5])
                    haskey(tim.stages, stage_name) || error("Stage name $stage_name not specified in .tim file.")
                    # A scenario branching from the root at the first period (SAPHIR's files)
                    # sets the data of the second: the first stage is not random
                    stage = max(tim.stages[stage_name], 2)
                    # In a two-stage model, a scenario may also branch from one that branches
                    # from the root (at the second period): it takes its parent's values, then
                    # sets its own
                    parent = words[3] in ("ROOT", "'ROOT'") ? nothing : Symbol(words[3])
                    parent === nothing || (N == 2 && haskey(root_scenarios, parent)) ||
                        error("SCEN format currently only supported for two-stage models.")
                    prob = parse(Float64, words[4])
                    # The scenarios are the outcomes of a single block
                    block = SCENARIOS
                    idx = get(ranvec_indices[stage - 1], block, 0)
                    ran_vec = if idx > 0
                        ran_vec = random_vecs[stage - 1][idx]
                        # Sanity check
                        ran_vec isa BlockDiscrete || error("Random vector at $block has more than one specified distribution.")
                        ran_vec
                    else
                        # Create new DISCRETE random vector
                        ran_vec = BlockDiscrete(T, inclusion)
                        push!(random_vecs[stage - 1], ran_vec)
                        # Bookkeep
                        ranvec_indices[stage - 1][block] = length(random_vecs[stage - 1])
                        ran_vec
                    end
                    # Update probability
                    push!(ran_vec.support, parent === nothing ? Dict{RowCol,T}() : copy(ran_vec.support[root_scenarios[parent]]))
                    push!(ran_vec.probabilities, prob)
                    parent === nothing && (root_scenarios[Symbol(words[2])] = length(ran_vec.support))
                    continue
                else
                    read_entries!(random_vecs[stage - 1][ranvec_indices[stage - 1][block]], words, cor)
                end
            else
                throw(ArgumentError("$(mode) is not a valid sto file mode."))
            end
        end
    end
    # Return raw data
    return RawStoch(name, random_vars, random_vecs)
end
parse_sto(tim::RawTime, cor::RawCor, filename::AbstractString) = parse_sto(Float64, tim, cor, filename)

# The entries of a data line of a BLOCKS or SCENARIOS section, one or two (row, value) pairs of a
# column, in the latest outcome of `ran_vec`
function read_entries!(ran_vec::BlockDiscrete{T}, words, cor::RawCor) where T <: AbstractFloat
    length(words) in (3, 5) || error("Malformed sto file at line $words.")
    col = Col(words[1])
    col = col == cor.rhsname ? RHS : col
    for k in 2:2:length(words)
        ran_vec.support[end][(Row(words[k]), col)] = convert(T, parse(Float64, words[k + 1]))
    end
    return ran_vec
end

function uncertainty_template(sto::RawStoch{N}, cor::RawCor, model::LPData{T, Matrix{T}}, stage::Integer) where {N, T}
    1 <= stage <= N || error("$(stage.id + 1) not in range 2 to $(N + 1).")
    map = model.indexmap
    # Prepare data
    Δc  = zero(model.c₁)
    ΔA  = zero(model.A)
    Δb  = zero(model.b)
    Δd₁ = zero(model.d₁)
    ΔC  = zero(model.C)
    Δd₂ = zero(model.d₂)
    # Collect rows, columns, and inclusions
    uncertainty_structure = Dict{RowCol, InclusionType}()
    for ran_var in sto.random_variables[stage]
        uncertainty_structure[ran_var.rowcol] = ran_var.inclusion
    end
    for ran_vec in sto.random_vectors[stage]
        isempty(ran_vec.support[1]) && error("Block $block has empty support.")
        for (rowcol, val) in ran_vec.support[1]
            rowcol = (row, col)
            uncertainty_structure[rowcol] = ran_vec.inclusion
        end
    end
    # Loop over all scenarios
    for (rowcol, inclusion) in uncertainty_structure
        (row, col) = rowcol
        rowsense = cor.rows[row][3]
        if row == cor.objname
            # Objective
            j = map[rowcol][2]
            if inclusion == REPLACE
                Δc[j] = -model.c[j]
            else
                Δc[j] = model.c[j]
            end
        elseif col == RHS
            # RHS
            i = map[rowcol][1]
            if rowsense == EQ
                if inclusion == REPLACE
                    Δb[i] = -model.b[i]
                else
                    Δb[i] = model.b[i]
                end
            elseif rowsense == LEQ
                if inclusion == REPLACE
                    Δd₂[i] = -model.d₂[i]
                    haskey(cor.ranges, row) && (Δd₁[i] = -model.d₁[i])
                else
                    Δd₂[i] = model.d₂[i]
                    haskey(cor.ranges, row) && (Δd₁[i] = model.d₁[i])
                end
            elseif rowsense == GEQ
                if inclusion == REPLACE
                    Δd₁[i] = -model.d₁[i]
                    haskey(cor.ranges, row) && (Δd₂[i] = -model.d₂[i])
                else
                    Δd₁[i] = model.d₁[i]
                    haskey(cor.ranges, row) && (Δd₂[i] = model.d₂[i])
                end
            end
        else
            # Matrix coeffs
            (i, j, _) = map[rowcol]
            if rowsense == EQ
                if inclusion == REPLACE
                    ΔA[i,j] = -model.A[i,j]
                else
                    ΔA[i] = model.A[i,j]
                end
            else
                if inclusion == REPLACE
                    ΔC[i,j] = -model.C[i,j]
                else
                    ΔC[i,j] = model.C[i,j]
                end
            end
        end
    end
    # Return template
    return LPData(T, Δc, zero(T), ΔA, Δb, Δd₁, ΔC, Δd₂, model.lb, model.ub, model.is_binary, model.is_integer, model.indexmap)
end

function uncertainty_template(sto::RawStoch{N}, cor::RawCor, model::LPData{T, SparseMatrixCSC{T,Int}}, stage::Integer) where {N, T}
    1 <= stage <= N || error("$(stage.id + 1) not in range 2 to $(N + 1).")
    map = model.indexmap
    # Get sizes
    n   = size(model.A, 2)
    m₁  = size(model.A, 1)
    m₂  = size(model.C, 1)
    # Prepare vectors
    Δc  = zero(model.c₁)
    Δb  = zero(model.b)
    Δd₁ = zero(model.d₁)
    Δd₂ = zero(model.d₂)
    # Prepare matrices
    ΔAᵢ = Vector{Int}()
    ΔAⱼ = Vector{Int}()
    ΔAᵥ = Vector{T}()
    ΔCᵢ = Vector{Int}()
    ΔCⱼ = Vector{Int}()
    ΔCᵥ = Vector{T}()
    # Collect rows, columns, and inclusions
    uncertainty_structure = Dict{RowCol, InclusionType}()
    for ran_var in sto.random_variables[stage]
        uncertainty_structure[ran_var.rowcol] = ran_var.inclusion
    end
    for ran_vec in sto.random_vectors[stage]
        isempty(ran_vec.support[1]) && error("Block $block has empty support.")
        for (rowcol, val) in ran_vec.support[1]
            uncertainty_structure[rowcol] = ran_vec.inclusion
        end
    end
    # Loop over all scenarios
    for (rowcol, inclusion) in uncertainty_structure
        (row, col) = rowcol
        rowsense = cor.rows[row][3]
        if row == cor.objname
            # Objective
            j = map[rowcol][2]
            if inclusion == REPLACE
                Δc[j] = -model.c₁[j]
            else
                Δc[j] = model.c₁[j]
            end
        elseif col == RHS
            # RHS
            i = map[rowcol][1]
            if rowsense == EQ
                if inclusion == REPLACE
                    Δb[i] = -model.b[i]
                else
                    Δb[i] = model.b[i]
                end
            elseif rowsense == LEQ
                if inclusion == REPLACE
                    Δd₂[i] = -model.d₂[i]
                    haskey(cor.ranges, row) && (Δd₁[i] = -model.d₂[i])
                else
                    Δd₂[i] = model.d₂[i]
                    haskey(cor.ranges, row) && (Δd₁[i] = model.d₂[i])
                end
            elseif rowsense == GEQ
                if inclusion == REPLACE
                    Δd₁[i] = -model.d₁[i]
                    haskey(cor.ranges, row) && (Δd₂[i] = -model.d₁[i])
                else
                    Δd₁[i] = model.d₁[i]
                    haskey(cor.ranges, row) && (Δd₂[i] = model.d₁[i])
                end
            end
        else
            # Matrix coeffs
            (i, j, _) = map[rowcol]
            if rowsense == EQ
                push!(ΔAᵢ, i)
                push!(ΔAⱼ, j)
                if inclusion == REPLACE
                    push!(ΔAᵥ, -model.A[i,j])
                else
                    push!(ΔAᵥ, model.A[i,j])
                end
            else
                push!(ΔCᵢ, i)
                push!(ΔCⱼ, j)
                if inclusion == REPLACE
                    push!(ΔCᵥ, -model.C[i,j])
                else
                    push!(ΔCᵥ, model.C[i,j])
                end
            end
        end
    end
    # Create sparse structures
    ΔA  = sparse(ΔAᵢ, ΔAⱼ, ΔAᵥ, m₁, n)
    ΔC  = sparse(ΔCᵢ, ΔCⱼ, ΔCᵥ, m₂, n)
    # Return template
    return LPData(T, Δc, zero(T), ΔA, Δb, Δd₁, ΔC, Δd₂, model.lb, model.ub, model.is_binary, model.is_integer, model.indexmap)
end

# The SMPS files of SPLIB and of Felt's test set, as they are written: several blocks, comment
# lines, two entries per line, scenarios branching at the first period, a single scenario. And
# what the reader refuses, as it always has: scenarios branching from a scenario, problem names
# that differ between files, probabilities that do not sum up to one. Each instance is `simple`
# with a file of its own.

"""The path of the instance made of `simple`'s files, its stoch file replaced by `sto` and its
core file by `cor`."""
function simple_with(sto::AbstractString; cor::AbstractString = read("io/smps/simple.cor", String))
    dir = mktempdir()
    write(joinpath(dir, "instance.cor"), cor)
    cp("io/smps/simple.tim", joinpath(dir, "instance.tim"))
    write(joinpath(dir, "instance.sto"), sto)
    touch(joinpath(dir, "instance.smps"))
    return joinpath(dir, "instance.smps")
end

support_of(path) = StochasticPrograms.SMPS.full_support(read(path, SMPSSampler))

@info "Running SMPS tests..."
@testset "SMPS files" begin
    simple = support_of("io/smps/simple.smps")
    @testset "Comment lines" begin
        cor = "* a comment line before the name\n" *
              replace(read("io/smps/simple.cor", String), "COLUMNS\n" => "COLUMNS\n* a comment line\n\n")
        sto = replace(read("io/smps/simple.sto", String), "BLOCKS" => "* a comment line\nBLOCKS")
        commented = support_of(simple_with(sto; cor = cor))
        @test probability.(commented) ≈ probability.(simple)
        @test [(ξ.Δq, ξ.Δd₂) for ξ in commented] == [(ξ.Δq, ξ.Δd₂) for ξ in simple]
    end
    @testset "Independent blocks" begin
        # simple's costs and bounds, drawn independently: their four combinations
        blocks = support_of(simple_with("""
            STOCH         SIMPLE
            BLOCKS        DISCRETE
             BL COSTS     STAGE2    0.4
                Y1        OBJ       -24.0
                Y2        OBJ       -28.0
             BL COSTS     STAGE2    0.6
                Y1        OBJ       -28.0
                Y2        OBJ       -32.0
             BL BOUNDS    STAGE2    0.5
                RHS       Y1UP      500.0
                RHS       Y2UP      100.0
             BL BOUNDS    STAGE2    0.5
                RHS       Y1UP      300.0
                RHS       Y2UP      300.0
            ENDATA
            """))
        @test length(blocks) == 4
        @test sort(probability.(blocks)) ≈ [0.2, 0.2, 0.3, 0.3]
        @test Set((ξ.Δq, ξ.Δd₂) for ξ in blocks) == Set((a.Δq, b.Δd₂) for a in simple, b in simple)
    end
    @testset "Two entries on a line" begin
        paired = support_of(simple_with("""
            STOCH         SIMPLE
            BLOCKS        DISCRETE
             BL BLOCK1    STAGE2    0.4
                Y1        OBJ       -24.0
                Y2        OBJ       -28.0
                RHS       Y1UP      500.0          Y2UP      100.0
             BL BLOCK1    STAGE2    0.6
                Y1        OBJ       -28.0
                Y2        OBJ       -32.0
                RHS       Y1UP      300.0          Y2UP      300.0
            ENDATA
            """))
        @test probability.(paired) ≈ probability.(simple)
        @test [(ξ.Δq, ξ.Δd₂) for ξ in paired] == [(ξ.Δq, ξ.Δd₂) for ξ in simple]
    end
    @testset "Scenarios branching from 'ROOT'" begin
        scenarios = support_of(simple_with(read("io/smps/simple_scen.sto", String)))
        @test probability.(scenarios) ≈ probability.(simple)
        @test [(ξ.Δq, ξ.Δd₂) for ξ in scenarios] == [(ξ.Δq, ξ.Δd₂) for ξ in simple]
    end
    @testset "Scenarios branching at the first period" begin
        sto = replace(read("io/smps/simple_scen.sto", String), "STAGE2" => "STAGE1")
        scenarios = support_of(simple_with(sto))
        @test probability.(scenarios) ≈ probability.(simple)
        @test [(ξ.Δq, ξ.Δd₂) for ξ in scenarios] == [(ξ.Δq, ξ.Δd₂) for ξ in simple]
    end
    @testset "Scenarios branching from a scenario" begin
        # as SAPHIR's: refused, the reader takes scenarios branching from the root only
        path = simple_with("""
            STOCH         SIMPLE
            SCENARIOS     DISCRETE
             SC SCEN1     ROOT      0.4            STAGE1
                Y1        OBJ       -24.0
                Y2        OBJ       -28.0
                RHS       Y1UP      500.0
                RHS       Y2UP      100.0
             SC SCEN2     SCEN1     0.6            STAGE2
                Y1        OBJ       -28.0
                Y2        OBJ       -32.0
            ENDATA
            """)
        @test_throws ErrorException support_of(path)
    end
    @testset "Names that differ" begin
        path = simple_with(replace(read("io/smps/simple.sto", String), "STOCH         SIMPLE" => "STOCH         OTHER"))
        @test_throws ErrorException read(path, SMPSSampler)
    end
    @testset "Probabilities that do not sum up to one" begin
        path = simple_with("""
            STOCH         SIMPLE
            INDEP         DISCRETE
                RHS       Y1UP      300.0          STAGE2    0.33333
                RHS       Y1UP      400.0          STAGE2    0.33333
                RHS       Y1UP      500.0          STAGE2    0.33333
            ENDATA
            """)
        @test_throws DomainError read(path, SMPSSampler)
    end
    @testset "A single scenario" begin
        one = support_of(simple_with("""
            STOCH         SIMPLE
            BLOCKS        DISCRETE
             BL BLOCK1    STAGE2    1.0
                Y1        OBJ       -24.0
            ENDATA
            """))
        @test length(one) == 1
        @test probability(only(one)) ≈ 1.0
    end
end

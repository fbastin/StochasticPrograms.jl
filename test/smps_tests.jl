# The SMPS files of SPLIB and of Felt's test set, as they are written: several blocks, comment
# lines, two entries per line, problem names that differ between files, probabilities rounded
# to a few digits, a single scenario. Each instance is `simple` with a file of its own.

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
    @testset "Names that differ" begin
        path = simple_with(replace(read("io/smps/simple.sto", String), "STOCH         SIMPLE" => "STOCH         OTHER"))
        sampler = @test_logs (:warn, r"name different problems") read(path, SMPSSampler)
        @test length(StochasticPrograms.SMPS.full_support(sampler)) == 2
    end
    @testset "Rounded probabilities" begin
        path = simple_with("""
            STOCH         SIMPLE
            INDEP         DISCRETE
                RHS       Y1UP      300.0          STAGE2    0.33333
                RHS       Y1UP      400.0          STAGE2    0.33333
                RHS       Y1UP      500.0          STAGE2    0.33333
            ENDATA
            """)
        sampler = @test_logs (:warn, r"sum up to 0.9999") read(path, SMPSSampler)
        scenarios = StochasticPrograms.SMPS.full_support(sampler)
        @test length(scenarios) == 3
        @test sum(probability, scenarios) ≈ 0.99999
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

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

@info "Running functionality tests..."
@testset "Stochastic Programs: Functionality" begin
    tol = 1e-2
    for (model, _scenarios, res, name) in problems
        sp = instantiate(model, _scenarios, optimizer = GLPK.Optimizer)
        @testset "SP Constructs: $name" begin
            optimize!(sp, cache = true)
            @test termination_status(sp) == MOI.OPTIMAL
            @test isapprox(optimal_decision(sp), res.x̄, rtol = tol)
            for i in 1:num_scenarios(sp)
                @test isapprox(optimal_recourse_decision(sp, i), res.ȳ[i], rtol = tol)
            end
            @test isapprox(objective_value(sp), res.VRP, rtol = tol)
            @test isapprox(EWS(sp), res.EWS, rtol = tol)
            @test isapprox(EVPI(sp), res.EVPI, rtol = tol)
            @test isapprox(VSS(sp), res.VSS, rtol = tol)
            @test isapprox(EV(sp), res.EV, rtol = tol)
            @test isapprox(EEV(sp), res.EEV, rtol = tol)
        end
        @testset "Inequalities: $name" begin
            @test EWS(sp) <= VRP(sp)
            @test VRP(sp) <= EEV(sp)
            @test VSS(sp) >= 0
            @test EVPI(sp) >= 0
            @test VSS(sp) <= EEV(sp) - EV(sp)
            @test EVPI(sp) <= EEV(sp) - EV(sp)
        end
        @testset "Copying: $name" begin
            sp_copy = copy(sp, optimizer = GLPK.Optimizer)
            add_scenarios!(sp_copy, scenarios(sp))
            @test num_scenarios(sp_copy) == num_scenarios(sp)
            generate!(sp_copy)
            @test num_subproblems(sp_copy) == num_subproblems(sp)
            optimize!(sp)
            optimize!(sp_copy)
            @test termination_status(sp_copy) == MOI.OPTIMAL
            @test isapprox(optimal_decision(sp_copy), optimal_decision(sp), rtol = tol)
            for i in 1:num_scenarios(sp)
                @test isapprox(optimal_recourse_decision(sp_copy, i), optimal_recourse_decision(sp, i), rtol = tol)
            end
            @test isapprox(objective_value(sp_copy), objective_value(sp), rtol = tol)
            @test isapprox(EWS(sp_copy), EWS(sp), rtol = tol)
            @test isapprox(EVPI(sp_copy), EVPI(sp), rtol = tol)
            @test isapprox(VSS(sp_copy), VSS(sp), rtol = tol)
            @test isapprox(EV(sp_copy), EV(sp), rtol = tol)
            @test isapprox(EEV(sp_copy), EEV(sp), rtol = tol)
        end
    end
    @testset "Sampling" begin
        sampled_sp = instantiate(simple, sampler, 100, optimizer = GLPK.Optimizer)
        generate!(sampled_sp)
        @test num_scenarios(sampled_sp) == 100
        @test isapprox(stage_probability(sampled_sp), 1.0)
        StochasticPrograms.sample!(sampled_sp, sampler, 100)
        generate!(sampled_sp)
        @test num_scenarios(sampled_sp) == 200
        @test isapprox(stage_probability(sampled_sp), 1.0)
    end
    @testset "Instant" begin
        optimize!(simple_sp)
        @test termination_status(simple_sp) == MOI.OPTIMAL
        @test isapprox(optimal_decision(simple_sp), simple_res.x̄, rtol = tol)
        for i in 1:num_scenarios(simple_sp)
            @test isapprox(optimal_recourse_decision(simple_sp, i), simple_res.ȳ[i], rtol = tol)
        end
        @test isapprox(objective_value(simple_sp), simple_res.VRP, rtol = tol)
        @test isapprox(EWS(simple_sp), simple_res.EWS, rtol = tol)
        @test isapprox(EVPI(simple_sp), simple_res.EVPI, rtol = tol)
        @test isapprox(VSS(simple_sp), simple_res.VSS, rtol = tol)
        @test isapprox(EV(simple_sp), simple_res.EV, rtol = tol)
        @test isapprox(EEV(simple_sp), simple_res.EEV, rtol = tol)
    end
    @testset "SMPS" begin
        simple_smps = read("io/smps/simple.smps", StochasticProgram, optimizer = GLPK.Optimizer)
        optimize!(simple_smps)
        @test termination_status(simple_smps) == MOI.OPTIMAL
        @test isapprox(optimal_decision(simple_smps), simple_res.x̄, rtol = tol)
        for i in 1:num_scenarios(simple_smps)
            @test isapprox(optimal_recourse_decision(simple_smps, i), simple_res.ȳ[i], rtol = tol)
        end
        @test isapprox(objective_value(simple_smps), simple_res.VRP, rtol = tol)
        @test isapprox(EWS(simple_smps), simple_res.EWS, rtol = tol)
        @test isapprox(EVPI(simple_smps), simple_res.EVPI, rtol = tol)
        @test isapprox(VSS(simple_smps), simple_res.VSS, rtol = tol)
        @test isapprox(EV(simple_smps), simple_res.EV, rtol = tol)
        @test isapprox(EEV(simple_smps), simple_res.EEV, rtol = tol)
    end
    @testset "Shadow prices agree with JuMP" begin
        # Capacity x <= 4 bought at 1 per unit; the demand ξ is covered by
        # y[1] <= x, free, and y[2] at 3 per unit, at least one unit of which
        # is compulsory through the balance y[2] - s == 1 (s costs 1). Every
        # row has a nonzero multiplier in some scenario, and the balance row
        # takes both signs. `sign = -1` writes the same problem as a
        # maximization.
        for sign in (1.0, -1.0)
            sense = sign > 0 ? MOI.MIN_SENSE : MOI.MAX_SENSE
            shadow = @stochastic_model begin
                @stage 1 begin
                    @decision(model, x >= 0)
                    @constraint(model, capacity, x <= 4)
                    @objective(model, sense, sign * x)
                end
                @stage 2 begin
                    @known(model, x)
                    @uncertain ξ
                    @recourse(model, y[1:2] >= 0)
                    @recourse(model, s >= 0)
                    @constraint(model, link, y[1] <= x)
                    @constraint(model, demand, y[1] + y[2] >= ξ)
                    @constraint(model, balance, y[2] - s == 1)
                    @objective(model, sense, sign * (3y[2] + s))
                end
            end
            demands = [4.0, 6.0]
            sp = instantiate(shadow, [@scenario(ξ = d, probability = 0.5) for d in demands],
                             optimizer = GLPK.Optimizer)
            optimize!(sp)
            @test termination_status(sp) == MOI.OPTIMAL
            # the extensive form, written out as a plain JuMP model
            m = Model(GLPK.Optimizer)
            @variable(m, x >= 0)
            capacity = @constraint(m, x <= 4)
            @variable(m, y[1:2, 1:2] >= 0)
            @variable(m, s[1:2] >= 0)
            link = [@constraint(m, y[1, k] <= x) for k in 1:2]
            demand = [@constraint(m, y[1, k] + y[2, k] >= demands[k]) for k in 1:2]
            balance = [@constraint(m, y[2, k] - s[k] == 1) for k in 1:2]
            @objective(m, sense, sign * (x + sum(0.5 * (3y[2, k] + s[k]) for k in 1:2)))
            optimize!(m)
            @test isapprox(shadow_price(sp[1, :capacity]), shadow_price(capacity), atol = 1e-8)
            for (name, refs) in ((:link, link), (:demand, demand), (:balance, balance)), k in 1:2
                @test isapprox(shadow_price(sp[2, name], k), shadow_price(refs[k]), atol = 1e-8)
                # relaxing a row never worsens the objective
                @test sign * shadow_price(sp[2, name], k) <= 1e-8
            end
            @test shadow_price(sp[1, :capacity]) ≈ -sign       # binding
            @test shadow_price(sp[2, :balance], 1) ≈ -1.5sign  # dual of either sign
            @test shadow_price(sp[2, :balance], 2) ≈ -0.5sign
        end
    end
end

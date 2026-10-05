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

function MOIU.map_indices(index_map::Function, change::DecisionCoefficientChange)
    return DecisionCoefficientChange(index_map(change.decision), change.new_coefficient)
end
MOIU.map_indices(index_map::Function, change::DecisionStateChange) = change
MOIU.map_indices(index_map::Function, change::KnownValuesChange) = change

# `modify_function!` must modify `f` in place: MOI's generic `modify_function`,
# which MOIU.UniversalFallback applies to the objective of a cached model, calls
# it on a copy of `f` and returns that copy, ignoring what it returns. Building
# and returning a new function instead silently drops the change; on the
# L-shaped master, it dropped θ from the objective whenever the first cuts came
# before the master had been attached to its optimizer.
function MOIU.modify_function!(f::AffineDecisionFunction{T}, change::MOI.ScalarConstantChange) where T
    MOIU.modify_function!(f.variable_part, MOI.ScalarConstantChange(T(change.new_constant)))
    return f
end
function MOIU.modify_function!(f::AffineDecisionFunction{T}, change::MOI.ScalarCoefficientChange) where T
    MOIU.modify_function!(f.variable_part, MOI.ScalarCoefficientChange(change.variable, T(change.new_coefficient)))
    return f
end

function MOIU.modify_function!(f::QuadraticDecisionFunction{T}, change::MOI.ScalarConstantChange) where T
    MOIU.modify_function!(f.variable_part, MOI.ScalarConstantChange(T(change.new_constant)))
    return f
end
function MOIU.modify_function!(f::QuadraticDecisionFunction{T}, change::MOI.ScalarCoefficientChange) where T
    MOIU.modify_function!(f.variable_part, MOI.ScalarCoefficientChange(change.variable, T(change.new_coefficient)))
    return f
end

function MOIU.modify_function!(f::VectorAffineDecisionFunction{T}, change::MOI.VectorConstantChange) where T
    MOIU.modify_function!(f.variable_part, MOI.VectorConstantChange(T.(change.new_constant)))
    return f
end
function MOIU.modify_function!(f::VectorAffineDecisionFunction{T}, change::MOI.MultirowChange) where T
    MOIU.modify_function!(f.variable_part, MOI.MultirowChange(change.variable, Tuple{Int64, T}[(idx, T(val)) for (idx, val) in change.new_coefficients]))
    return f
end

function MOIU.modify_function!(f::AffineDecisionFunction{T}, change::DecisionCoefficientChange) where T
    MOIU.modify_function!(f.decision_part, MOI.ScalarCoefficientChange(change.decision, T(change.new_coefficient)))
    return f
end

function MOIU.modify_function!(f::QuadraticDecisionFunction{T}, change::DecisionCoefficientChange) where T
    MOIU.modify_function!(f.decision_part, MOI.ScalarCoefficientChange(change.decision, T(change.new_coefficient)))
    return f
end

function MOIU.modify_function!(f::VectorAffineDecisionFunction{T}, change::DecisionMultirowChange) where T
    MOIU.modify_function!(f.decision_part, MOI.MultirowChange(change.decision, Tuple{Int64, T}[(idx, T(val)) for (idx, val) in change.new_coefficients]))
    return f
end

function MOIU.modify_function!(f::Union{SingleDecision, AffineDecisionFunction, QuadraticDecisionFunction, VectorAffineDecisionFunction}, change::Union{DecisionStateChange, KnownValuesChange})
    # Nothing to do here, handled in bridges
    return f
end

# Can rely on scalarize bridge if modifications are passed along
function MOI.modify(model::MOI.ModelLike,
                    bridge::MOIB.Constraint.ScalarizeBridge{T, AffineDecisionFunction{T}},
                    change::DecisionModification) where T
    for constraint in bridge.scalar_constraints
        MOI.modify(model, constraint, change)
    end
    return nothing
end

# Can rely on vectorize bridge if modifications are passed along
function MOI.modify(model::MOI.ModelLike,
                    bridge::MOIB.Constraint.VectorizeBridge{T, VectorAffineDecisionFunction{T}},
                    change::DecisionModification) where T
    MOI.modify(model, bridge.vector_constraint, change)
    return nothing
end

# Can rely on flipsign bridges if modifications are passed along
function MOI.modify(model::MOI.ModelLike,
                    bridge::MOIB.Constraint.FlipSignBridge,
                    change::DecisionModification)
    MOI.modify(model, bridge.constraint, change)
    return nothing
end

# Can rely on norm-bridges if modifications are passed along
function MOI.modify(model::MOI.ModelLike,
                    bridge::MOIB.Constraint.NormInfinityBridge{T, VectorAffineDecisionFunction{T}},
                    change::DecisionModification) where T
    MOI.modify(model, bridge.constraint, change)
    return nothing
end
function MOI.modify(model::MOI.ModelLike,
                    bridge::MOIB.Constraint.NormOneBridge{T, VectorAffineDecisionFunction{T}},
                    change::DecisionModification) where T
    MOI.modify(model, bridge.nn_index, change)
    return nothing
end

function MOI.modify(uf::MOIU.UniversalFallback, ::CI{F,S}, ::Union{DecisionStateChange,KnownValuesChange}) where {F, S}
    return nothing
end

function MOI.get(model::MOIU.CachingOptimizer,
                 attr::DecisionIndex,
                 ci::CI)
    return MOI.get(model.optimizer, attr, model.model_to_optimizer_map[ci])
end

function MOI.get(b::MOIB.AbstractBridgeOptimizer,
                 attr::DecisionIndex, ci::MOI.ConstraintIndex)
    return MOIB.call_in_context(b, ci, bridge -> MOI.get(b, attr, bridge))
end

function MOIU.shift_constant(f::AffineDecisionFunction{T}, offset) where T
    return typeof(f)(MOIU.shift_constant(f.variable_part, offset), copy(f.decision_part))
end

function MOIU.shift_constant(f::QuadraticDecisionFunction{T}, offset) where T
    return typeof(f)(MOIU.shift_constant(f.variable_part, offset), copy(f.decision_part), copy(f.cross_terms))
end

function MOIU.shift_constant(f::VectorAffineDecisionFunction{T}, offset) where T
    return typeof(f)(MOIU.shift_constant(f.variable_part, offset), copy(f.decision_part))
end

# Move the constant of a scalar decision function to the set, as MOI does for
# affine and quadratic functions. MOI's fallback leaves the function alone, and
# scalar constraints with a nonzero constant are refused by `add_constraint`:
# the ScalarizeBridge relies on this to split a vector constraint into rows.
# The constant of a decision function lives in its variable part.
function MOIU.normalize_constant(func::Union{AffineDecisionFunction{T}, QuadraticDecisionFunction{T}},
                                 set::MOI.AbstractScalarSet;
                                 allow_modify_function::Bool = false) where T
    if MOIU.supports_shift_constant(typeof(set))
        set = MOIU.shift_constant(set, -MOI.constant(func))
        if !allow_modify_function
            func = copy(func)
        end
        func.variable_part.constant = zero(T)
    end
    return func, set
end

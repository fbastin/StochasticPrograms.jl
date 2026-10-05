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

# Report a syntax error in a macro call, as `JuMP._macro_error` did before it
# was removed from JuMP: where the call is, the call itself, and what is wrong.
function _macro_error(macroname, args, source, str...)
    location = source isa LineNumberNode ? "At $(source.file):$(source.line): " : ""
    error(location, "`@$macroname($(join(args, ", ")))`: ", str...)
end

# The same for an expression `x` met inside the body of `@macroname`: a macro
# call is reported as itself, anything else as part of the enclosing call.
function _macro_error(x, macroname::Symbol, source, str...)
    if Meta.isexpr(x, :macrocall)
        name = Symbol(String(x.args[1])[2:end])
        return _macro_error(name, prettify.(x.args[3:end]), x.args[2], str...)
    end
    return _macro_error(macroname, [prettify(x)], source, str...)
end

include("scenario.jl")
include("decisions.jl")
include("define_scenario.jl")
include("sampler.jl")
include("stage.jl")
include("model.jl")

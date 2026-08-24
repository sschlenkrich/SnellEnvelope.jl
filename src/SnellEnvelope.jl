module SnellEnvelope

using Distributions
using Interpolations
using QuadGK
using SpecialFunctions


@enum ModelType begin
    normal
    lognormal
end

include("ImpliedVolatility.jl")
include("RollBack.jl")
include("BackwardInduction.jl")



end

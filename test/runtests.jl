using SnellEnvelope
using Test

@testset "SnellEnvelope.jl" begin

    @info "Running tests for SnellEnvelope.jl"

    include("implied_volatility.jl")
    include("roll_back.jl")
    include("backward_induction.jl")
    include("snell_envelope.jl")

end

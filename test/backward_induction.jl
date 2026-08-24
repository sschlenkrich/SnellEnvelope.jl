using SnellEnvelope
using Test

@testset "Backward induction tests" begin

    @testset "Normal model single step" begin
        F = 1.0
        for K in range(0.5, 1.5, length=11)
            U(x) = max((F + x) - K, 0.0)  # payoff function
            σ = 0.2 * sqrt(10.0)  # standard deviation
            λ = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.backward_induction([U], [σ], SnellEnvelope.normal, λ, n_points)
            V = SnellEnvelope.bachelier_price(K, F, σ, +1)
            @test isapprox(H0, V; atol=1e-5)
        end
    end
    
    @testset "Normal model multi-step" begin
        F = 1.0
        K = 1.0
        U(x) = max((F + x) - K, 0.0)  # payoff function
        Zero(x) = 0.0
        λ = 5.0
        n_points = 101
        #
        σ1 = 0.2 * sqrt(2.0)  # standard deviation to exercise
        σ2 = 0.2 * sqrt(5.0)
        σ3 = 0.2 * sqrt(10.0)
        #
        H0 = SnellEnvelope.backward_induction([Zero, Zero, U], [σ1, σ2, σ3], SnellEnvelope.normal, λ, n_points)
        V = SnellEnvelope.bachelier_price(K, F, σ3, +1)
        @test isapprox(H0, V; atol=1e-5)
    end

    @testset "Log-normal model single step" begin
        F = 1.0
        for K in range(0.5, 1.5, length=11)
            U(x) = max((F * exp(x)) - K, 0.0)  # payoff function
            σ = 0.2 * sqrt(10.0)  # standard deviation
            λ = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.backward_induction([U], [σ], SnellEnvelope.lognormal, λ, n_points)
            V = SnellEnvelope.black_price(K, F, σ, +1)
            @test isapprox(H0, V; atol=1e-5)
        end
    end

    @testset "Log-normal model multi-step" begin
        F = 1.0
        K = 1.0
        U(x) = max((F * exp(x)) - K, 0.0)  # payoff function
        Zero(x) = 0.0
        λ = 5.0
        n_points = 101
        #
        σ1 = 0.2 * sqrt(2.0)  # standard deviation to exercise
        σ2 = 0.2 * sqrt(5.0)
        σ3 = 0.2 * sqrt(10.0)
        #
        H0 = SnellEnvelope.backward_induction([Zero, Zero, U], [σ1, σ2, σ3], SnellEnvelope.lognormal, λ, n_points)
        V = SnellEnvelope.black_price(K, F, σ3, +1)
        @test isapprox(H0, V; atol=1e-5)
    end
end
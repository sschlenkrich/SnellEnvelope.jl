using SnellEnvelope
using Test

@testset "Conditional expectation via roll-back" begin

    @testset "Normal model single step" begin
        F = 1.0
        for K in range(0.5, 1.5, length=11)
            H(x) = 0.0
            U(x) = max((F + x) - K, 0.0)  # payoff function
            μ = 0.0
            σ = 0.2 * sqrt(10.0)  # standard deviation
            λ = 5.0
            x0 = 0.0
            #
            H0 = SnellEnvelope.roll_back(x0, H, U, μ, σ, λ)
            V = SnellEnvelope.bachelier_price(K, F, σ, +1)
            @test isapprox(H0, V; atol=1e-5)
        end
    end

    @testset "Normal model multi-step" begin
        F = 1.0
        K = 1.0
        H(x) = 0.0
        U(x) = max((F + x) - K, 0.0)  # payoff function
        μ = 0.0
        λ = 5.0
        x0 = 0.0
        #
        σ1 = 0.2 * sqrt(2.0)  # standard deviation
        σ2 = 0.2 * sqrt(3.0)
        σ3 = 0.2 * sqrt(5.0)
        #
        n_points = 101
        H2 = SnellEnvelope.roll_back(x0, sqrt(σ1^2 + σ2^2), H, U, μ, σ3, λ, n_points)
        H1 = SnellEnvelope.roll_back(x0, σ1, H, H2, μ, σ2, λ, n_points)
        H0 = SnellEnvelope.roll_back(x0, H, H1, μ, σ1, λ)
        V = SnellEnvelope.bachelier_price(K, F, sqrt(σ1^2 + σ2^2 + σ3^2), +1)
        @test isapprox(H0, V; atol=1e-5)
    end

    @testset "Log-normal model single step" begin
        F = 1.0
        for K in range(0.5, 1.5, length=11)
            H(x) = 0.0
            U(x) = max((F * exp(x)) - K, 0.0)  # payoff function
            σ = 0.2 * sqrt(10.0)  # standard deviation
            μ = -0.5 * σ^2  # drift
            λ = 5.0
            x0 = 0.0
            #
            H0 = SnellEnvelope.roll_back(x0, H, U, μ, σ, λ)
            V = SnellEnvelope.black_price(K, F, σ, +1)
            @test isapprox(H0, V; atol=1e-5)
        end
    end

    @testset "Log-normal model multi-step" begin
        F = 1.0
        K = 1.0
        H(x) = 0.0
        U(x) = max((F * exp(x)) - K, 0.0)  # payoff function
        λ = 5.0
        x0 = 0.0
        #
        σ1 = 0.2 * sqrt(2.0)  # standard deviation
        σ2 = 0.2 * sqrt(3.0)
        σ3 = 0.2 * sqrt(5.0)
        μ1 = -0.5 * σ1^2
        μ2 = -0.5 * σ2^2
        μ3 = -0.5 * σ3^2
        #
        n_points = 101
        H3 = SnellEnvelope.roll_back(x0 + μ1 + μ2, sqrt(σ1^2 + σ2^2), H, U, μ3, σ3, λ, n_points)
        H2 = SnellEnvelope.roll_back(x0 + μ1, σ1, H, H3, μ2, σ2, λ, n_points)
        H1 = SnellEnvelope.roll_back(x0, H, H2, μ1, σ1, λ)
        V = SnellEnvelope.black_price(K, F, sqrt(σ1^2 + σ2^2 + σ3^2), +1)
        @test isapprox(H1, V; atol=1e-5)
    end

end
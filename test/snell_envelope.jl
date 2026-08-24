using SnellEnvelope
using Test

@testset "Snell envelope tests" begin

    @testset "Normal model single step" begin
        A = 5.0
        F = 0.05
        σ = 0.01 * sqrt(10.0)
        cp = 1.0
        for K in range(0.00, 0.10, length=11)
            V0 = A * SnellEnvelope.bachelier_price(K, F, σ, cp)
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope([F], [K], [A], [V0], cp, SnellEnvelope.normal, lambda, n_points)
            @test isapprox(H0, V0; atol=1e-5)
        end
    end

    @testset "Normal model multi-step" begin
        F = 0.05
        A = [0.01, 0.01, 5.0]
        σ = 0.01 .* sqrt.([2.0, 3.0, 10.0])
        cp = 1.0
        for K in range(0.00, 0.10, length=11)
            V0 = [
                A_ * SnellEnvelope.bachelier_price(K, F, σ_, cp)
                for (A_, σ_) in zip(A, σ)
            ]
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope([F, F, F], [K, K, K], A, V0, cp, SnellEnvelope.normal, lambda, n_points)
            @test isapprox(H0, V0[end]; atol=1e-5)
        end
    end

    @testset "Log-normal model single step" begin
        A = 5.0
        F = 0.05
        σ = 0.25 * sqrt(10.0)
        cp = 1.0
        for K in range(0.01, 0.10, length=10)
            V0 = A * SnellEnvelope.black_price(K, F, σ, cp)
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope([F], [K], [A], [V0], cp, SnellEnvelope.lognormal, lambda, n_points)
            @test isapprox(H0, V0; atol=1e-5)
        end
    end

    @testset "Log-normal model multi-step" begin
        F = 0.05
        A = [0.01, 0.01, 5.0]
        σ = 0.25 .* sqrt.([2.0, 3.0, 10.0])
        cp = 1.0
        for K in range(0.01, 0.10, length=10)
            V0 = [
                A_ * SnellEnvelope.black_price(K, F, σ_, cp)
                for (A_, σ_) in zip(A, σ)
            ]
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope([F, F, F], [K, K, K], A, V0, cp, SnellEnvelope.lognormal, lambda, n_points)
            @test isapprox(H0, V0[end]; atol=1e-5)
        end
    end

    @testset "Positive switch value normal model" begin
        F = 0.05
        A = [5.0, 4.0, 3.0, 2.0, 1.0]
        σ = 0.01 .* sqrt.([1.0, 2.0, 3.0, 4.0, 5.0])
        cp = 1.0
        for K in range(0.01, 0.10, length=10)
            V0 = [
                A_ * SnellEnvelope.bachelier_price(K, F, σ_, cp)
                for (A_, σ_) in zip(A, σ)
            ]
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope(F.*ones(5), K.*ones(5), A, V0, cp, SnellEnvelope.normal, lambda, n_points)
            SW = SnellEnvelope.swaption_switch_value(F.*ones(5), K.*ones(5), A, V0, cp, SnellEnvelope.normal, lambda, n_points)
            @test H0 > maximum(V0)
            @test SW == H0 - maximum(V0)
        end
    end

    @testset "Positive switch value log-normal model" begin
        F = 0.05
        A = [5.0, 4.0, 3.0, 2.0, 1.0]
        σ = 0.25 .* sqrt.([1.0, 2.0, 3.0, 4.0, 5.0])
        cp = -1.0
        for K in range(0.01, 0.10, length=10)
            V0 = [
                A_ * SnellEnvelope.black_price(K, F, σ_, cp)
                for (A_, σ_) in zip(A, σ)
            ]
            #
            lambda = 5.0
            n_points = 101
            #
            H0 = SnellEnvelope.swaption_snell_envelope(F.*ones(5), K.*ones(5), A, V0, cp, SnellEnvelope.lognormal, lambda, n_points)
            SW = SnellEnvelope.swaption_switch_value(F.*ones(5), K.*ones(5), A, V0, cp, SnellEnvelope.lognormal, lambda, n_points)
            @test H0 > maximum(V0)
            @test SW == H0 - maximum(V0)
        end
    end
end
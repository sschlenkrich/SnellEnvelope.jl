using SnellEnvelope
using Test

@testset "Implied Volatility" begin

    @testset "Normal volatility" begin
        K = 1.0
        F = range(0.5, 1.5, length=11)
        σ = 0.2
        T = 10.0
        call_put = +1
        bachelier_price = [
            SnellEnvelope.bachelier_price(K, F_, σ, T, call_put)
            for F_ in F
        ]
        σ_iv = [
            SnellEnvelope.bachelier_implied_volatility(price, K, F_, T, call_put)
            for (price, F_) in zip(bachelier_price, F)
        ]
        @test all(abs.(σ_iv .- σ) .< 1e-16)
    end

    @testset "Log-normal volatility" begin
        K = 1.0
        F = range(0.5, 1.5, length=11)
        σ = 0.2
        T = 10.0
        call_put = +1
        black_price = [
            SnellEnvelope.black_price(K, F_, σ, T, call_put)
            for F_ in F
        ]
        σ_iv = [
            SnellEnvelope.black_implied_volatility(price, K, F_, T, call_put)
            for (price, F_) in zip(black_price, F)
        ]
        @test all(abs.(σ_iv .- σ) .< 1e-15)
    end

end

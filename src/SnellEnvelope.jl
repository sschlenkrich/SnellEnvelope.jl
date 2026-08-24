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


"""
    swaption_snell_envelope(
        swap_rates::Vector{Float64},
        strike_rates::Vector{Float64},
        annuities::Vector{Float64},
        prices::Vector{Float64},
        call_put::Float64,
        model_type::ModelType,
        lambda::Float64,
        n_points::Int64,
        )

Calculate the Snell envelope for a sequence of European swaptions with
given swap rates, strike rates, annuities, and prices.

The `call_put` parameter indicates whether the options are calls (+1)
or puts (-1).

The `model_type` specifies whether to use the normal (Bachelier) or
lognormal (Black) model.

The `lambda` parameter specifies the number of standard deviations for
integration, and `n_points` specifies the number of points for interpolation.
"""
function swaption_snell_envelope(
    swap_rates::Vector{Float64},
    strike_rates::Vector{Float64},
    annuities::Vector{Float64},
    prices::Vector{Float64},
    call_put::Float64,
    model_type::ModelType,
    lambda::Float64,
    n_points::Int64,
    )
    #
    @assert length(swap_rates) == length(strike_rates) == length(annuities) == length(prices)
    if model_type == normal
        implied_volatility = bachelier_implied_volatility
    elseif model_type == lognormal
        implied_volatility = black_implied_volatility
    else
        error("Unsupported model type: $model_type")
    end
    std_devs = [
        implied_volatility(price / annuity, strike, forward, 1.0, call_put)
        for (price, strike, forward, annuity) in zip(prices, strike_rates, swap_rates, annuities)
    ]
    #
    if model_type == normal
        payoffs = [
            x -> A * max(call_put * ((F + x) - K), 0.0)
            for (A, F, K) in zip(annuities, swap_rates, strike_rates)
        ]
    elseif model_type == lognormal
        payoffs = [
            x -> A * max(call_put * ((F * exp(x)) - K), 0.0)
            for (A, F, K) in zip(annuities, swap_rates, strike_rates)
        ]
    else
        error("Unsupported model type: $model_type")
    end
    #
    H0 = backward_induction(payoffs, std_devs, model_type, lambda, n_points)
    return H0
end


"""
    swaption_switch_value(
        swap_rates::Vector{Float64},
        strike_rates::Vector{Float64},
        annuities::Vector{Float64},
        prices::Vector{Float64},
        call_put::Float64,
        model_type::ModelType,
        lambda::Float64,
        n_points::Int64,
        )

Calculate the switch value for a sequence of European swaptions.

For details on the inputs, see `swaption_snell_envelope`.

The switch value is defined as the difference between the Snell envelope and
the maximum of the prices of the European swaptions.

Switch value is floored at zero to capture small deviations due to numerical
errors for out-of-the-money swaptions.
"""
function swaption_switch_value(
    swap_rates::Vector{Float64},
    strike_rates::Vector{Float64},
    annuities::Vector{Float64},
    prices::Vector{Float64},
    call_put::Float64,
    model_type::ModelType,
    lambda::Float64,
    n_points::Int64,
    )
    H0 = swaption_snell_envelope(swap_rates, strike_rates, annuities, prices, call_put, model_type, lambda, n_points)
    return max(H0 - maximum(prices), 0.0)
end


end # module

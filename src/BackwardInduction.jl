

"""
    backward_induction(
    payoffs::Vector,
    std_devs::Vector{Float64},
    model_type::ModelType,
    lambda::Float64,
    n_points::Int,
    )

Calcuate the Snell envelope using backward induction.
# Arguments
- `payoffs`: A vector of payoff functions, one for each exercise date.
- `std_devs`: A vector of standard deviations, one for each exercise date.
- `model_type`: The model type, either `ModelType.normal` or `ModelType.lognormal`.
- `lambda`: The number of standard deviations to use for the integration range.
- `n_points`: The number of points to use for the integration.
"""
function backward_induction(
    payoffs::Vector,
    std_devs::Vector{Float64},
    model_type::ModelType,
    lambda::Float64,
    n_points::Int,
    )
    #
    @assert length(payoffs) > 0 "At least one payoff function must be provided."
    @assert length(payoffs) == length(std_devs) "Number of payoffs and standard deviations must match."
    @assert all(std_devs[2:end] - std_devs[1:end-1] .> 0.0) "Standard deviations must be strictly increasing."
    #
    σ = vcat(std_devs[begin], sqrt.(std_devs[2:end].^2 .- std_devs[1:end-1].^2))
    if model_type == normal
        μ = zeros(length(σ))
    elseif model_type == lognormal
        μ = -0.5 .* σ.^2
    else
        error("Unsupported model type: $model_type")
    end
    #
    H(x) = 0.0  # continuation value after last exercise
    for i in length(payoffs):-1:2
        if model_type == normal
            x0 = 0.0
        elseif model_type == lognormal
            x0 = -0.5 * std_devs[i-1]^2
        else
            error("Unsupported model type: $model_type")
        end
        H = roll_back(x0, std_devs[i-1], H, payoffs[i], μ[i], σ[i], lambda, n_points)
    end
    H0 = roll_back(0.0, H, payoffs[1], μ[1], σ[1], lambda)
    return H0
end

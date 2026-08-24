

"""

Calculate the conditional expectation.
"""
function roll_back(x, H, U, μ, σ, λ)
    p(u) = pdf(Normal(x + μ, σ), u)
    f(u) = max(H(u), U(u)) * p(u)
    x_min = x + μ - λ * σ
    x_max = x + μ + λ * σ
    integral, _ = quadgk(f, x_min, x_max, rtol=1e-8, atol=1e-8)
    return integral
end


"""
Calculate conditional expectation as a function of state variable.

We use polynomial interpolation in state space.
"""
function roll_back(x0, σ0, H, U, μ, σ, λ, n_points)
    x_range = range(-λ*σ0, +λ*σ0, n_points) .+ x0
    v0 = [
        roll_back(x, H, U, μ, σ, λ)
        for x in x_range
    ]
    H0 = extrapolate(
        interpolate(x_range, v0, AkimaMonotonicInterpolation()),
        Line(),
    )
    return H0
end


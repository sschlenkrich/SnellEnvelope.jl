# For reference, see https://github.com/frame-consulting/DiffFusion.jl/tree/main/src/utils

"""
    _d_1(strike, forward, nu)

Calculate the Black formula d_1 term.

We implement this once in order to fine-tune corner cases later.
"""
_d_1(strike, forward, nu) = log.(forward ./ strike) ./ nu .+ nu / 2.0


"""
    black_price(strike, forward, nu, call_put)

Calculate Vanilla option price V in Black model.

Argument `strike` (``K``) represents the option strike, `forward` (``F``) is the
forward price of the underlying ``S``, i.e. ``F = E[S(T)]`` in the respective
pricing measure.

Argument nu represents the standard deviation of the forward price. For
annualised lognormal volatility σ, we have ``ν = σ √T``. Finally, `call_put`
is the call (+1) or put (-1) option flag.

Note that we allow negative standard deviation ν. This allows modelling negative
skew for shifted-lognormal models.

Note that we do not test for vanishing volatility. Zero volatility may cause
floating point errors.

We allow broadcasting for arguments.
"""
function black_price(strike, forward, nu, call_put)
    d1 = _d_1(strike, forward, nu)
    d2 = d1 .- nu
    option_value = call_put * (
        forward .* cdf.(Normal(), call_put .* d1) .-
        strike  .* cdf.(Normal(), call_put .* d2)
        )
    intrinsic_value = max.(call_put * (forward .- strike), 0.0)
    return max.(option_value, intrinsic_value)
end

"""
    black_price(strike, forward, σ, T, call_put)

Calculate Vanilla option price ``V`` in Black model with volatility parameter.
"""
black_price(strike, forward, σ, T, call_put) = black_price(strike, forward, σ * sqrt(T), call_put)

"""
    bachelier_price(strike, forward, nu, call_put)

Calculate Vanilla option price ``V`` in Bachelier model.

Argument `strike` (``K``) represents the option strike, `forward` (``F``) is the
forward price of the underlying ``S``, i.e. ``F = E[S(T)]`` in the respective
pricing measure.

Argument `nu` represents the standard deviation of the forward price. For
annualised normal volatility σ, we have `ν = σ √T`. Finally, `call_put`
is the call (+1) or put (-1) option flag.

Note that we do not test for vanishing volatility. Zero volatility may cause
floating point errors.
"""
function bachelier_price(strike, forward, nu, call_put)
    u = call_put .* (forward .- strike)
    h = u ./ nu
    option_value = u .* cdf.(Normal(), h) .+ nu .* pdf.(Normal(), h)
    intrinsic_value = max.(u, 0.0)
    return max.(option_value, intrinsic_value)
end

"""
    bachelier_price(strike, forward, σ, T, call_put)

Calculate Vanilla option price ``V`` in Bachelier model with volatility parameter.
"""
bachelier_price(strike, forward, σ, T, call_put) = bachelier_price(strike, forward, σ * sqrt(T), call_put)


# Implied volatility calculations



_Φ_tilde(x) = cdf(Normal(), x) + pdf(Normal(), x) / x

"""
    bachelier_implied_volatility(price, strike, forward, T, call_put)

Calculate the Bachelier or normal implied volatility.

`price` is the forward or undiscounted option price, `strike` is the option strike,
`forward` is the forward asset price or expectation of factor, `T` is time to option
expiry, and `call_put` encodes call (`+1`) or put (`-1`) options.
"""
function bachelier_implied_volatility(price, strike, forward, T, call_put)
    # P. Jäckel, Implied Normal Volatility, 2017
    if isapprox(strike, forward)  # ATM option
        atm_vega = sqrt(T) / sqrt(2*π)
        return price / atm_vega
    end
    time_value = price - max(call_put*(forward-strike), 0.0)
    if isapprox(time_value, 0.0)  # option = intrinsic value
        return 0.0
    end
    @assert	time_value > 0.0  "time_value: " * string(time_value)
    Φ_tilde_star = -abs(time_value) / abs(forward - strike)
    # solve Φ_tilde(x_star) = Φ_tilde_star
    if Φ_tilde_star < -0.001882039271
        g = 1.0 / (Φ_tilde_star - 0.5)
        ξ_bar = 
            (0.032114372355 -
             g * g *
                 (0.016969777977 - g * g * (2.6207332461E-3 - 9.6066952861E-5 * g * g))) /
            (1.0 -
             g * g * (0.6635646938 - g * g * (0.14528712196 - 0.010472855461 * g * g)))
        x_bar = g * (0.3989422804014326 + ξ_bar * g * g)
    else
        h = sqrt(-log(-Φ_tilde_star))
        x_bar =
            (9.4883409779 - h * (9.6320903635 - h * (0.58556997323 + 2.1464093351 * h))) /
            (1.0 - h * (0.65174820867 + h * (1.5120247828 + 6.6437847132E-5 * h)))
    end
    q = (_Φ_tilde(x_bar) - Φ_tilde_star) / pdf(Normal(), x_bar)
    x_star = x_bar +
        3.0 * q * x_bar * x_bar * (2.0 - q * x_bar * (2.0 + x_bar * x_bar)) /
        (6.0 + q * x_bar * (-12.0 + x_bar * (6.0 * q + x_bar *
        (-6.0 + q * x_bar * (3.0 + x_bar * x_bar)))))
    σ = abs((forward - strike)/x_star) / sqrt(T)
    return σ
end



# This file contains methods to calculate implied log-normal volatilities from
# https://github.com/jherekhealy/AQFED.jl.
#
# We include the code here directly to avoid the dependency on the unregistered
# AQFED.jl package.
#
# For code reference, see
#
# https://github.com/jherekhealy/AQFED.jl/blob/master/src/black/iv_solver_common.jl
# https://github.com/jherekhealy/AQFED.jl/blob/master/src/black/iv_solver_householder.jl
#


"""
    black_implied_volatility(price, strike, forward, T, call_put)

Calculate the Black or log-normal implied volatility.

`price` is the forward or undiscounted option price, `strike` is the option strike,
`forward` is the forward asset price or expectation of factor, `T` is time to option
expiry, and `call_put` encodes call (`+1`) or put (`-1`) options.
"""
function black_implied_volatility(price, strike, forward, T, call_put)
    isCall = (call_put==1.0)
    df = 1.0
    ftolrel = 0.0  # see AQFED.jl
    maxEval = 64   # see AQFED.jl
    return impliedVolatilitySRHalley(isCall, price, forward, strike, T, df, ftolrel, maxEval)
end


function normalizePrice(isCall::Bool, price::TP, f::T, strike::T, df::T) where {T,TP}
    c = price / f / df
    ex = f / strike

    if !isCall
        if ex <= 1
            c = c + 1 - 1 / ex # put call parity
        else
            #else duality + put call parity
            c = ex * c
            ex = 1 / ex
        end
    else
        if ex > 1
            # use c(-x0, v)
            c = (f * (c - 1) + strike) / strike # in out duality, c = ex*c + 1 - ex  //not as good numericall
            ex = 1 / ex
        end
    end
    return c, ex
end

const SqrtEpsilon = sqrt(eps())

# f = 1.0, strike = 1/ex, ex=ey = 1.0/strike=f/strike, y = log(f/strike) = log(ey), alpha = price*ex
#polya is 2 / pi, polya_factor = sqrt(pi/8) is aludaat,
function impliedVolatilitySRGuessUndiscountedCall(price::Real, ey::Real, y::Real, polya_factor::Real = 2/pi)::Real
    alpha = price *ey
    r = 2 * alpha - ey + 1
    em2piy = exp(-polya_factor * y)
    A = (ey * em2piy - 1.0 / (ey * em2piy))^2
    r2 = r * r
    B = 4 * (em2piy + 1.0 / em2piy) - 2 / ey * (ey * em2piy + 1.0 / (ey * em2piy)) * ((ey * ey) + 1 - r2)
    beta = 0.0
    if abs(alpha) < SqrtEpsilon
        C = -16 * (1 - 1.0 / ey) * alpha - 16 * (1 - 3 / ey + 1 / (ey * ey)) * alpha * alpha
        beta = C / B
        #    println(" alpha ",alpha, " C ",C," beta ",beta)
    elseif abs(ey - alpha) < SqrtEpsilon
        C = -16 * (1 + 1.0 / ey) * (alpha - ey) - 16 * (1 + 3 / ey + 1 / (ey * ey)) * (alpha - ey) * (alpha - ey)
        if C == 0
            beta = B / A #this is wrong sign but gives order of magnitude
        else
            beta = C / B
        end
        #    println(" alpha ",alpha, " C ",C," beta ",beta)
    elseif abs(y) < SqrtEpsilon
        a2 = alpha * alpha
        a4 = a2 * a2
        C = 16 * ((a2 - a4) + (2 * a4 + 2 * alpha * a2 - a2 - alpha) * y)
        B = 16 * (a2 - y * (alpha + a2))
        beta = C / B
    else
        eym1 = ey - 1
        eyp1 = ey + 1
        C = 1.0 / (ey * ey) * (r2 - eym1 * eym1) * (eyp1 * eyp1 - r2)
        beta = 2 * C / (B + sqrt(B * B + 4 * A * C))
        # println("beta ",beta, " ",C, " ",B, " ",A)
    end
    if beta < 0
        beta = SqrtEpsilon
    end
    gamma = -log(beta) / polya_factor
    if y >= 0
        Asqrty = 0.5 * (1 + sqrt(1 - em2piy * em2piy))
        c0 = (Asqrty - 0.5/ey)
        gmy = max(gamma - y,0)      #machine epsilon issues
        if price <= c0
            return sqrt(gamma + y) - sqrt(gmy)
        end
        return (sqrt(gamma + y) + sqrt(gmy))
    end
    Asqrty = 0.5 * (1 - sqrt(1 - 1.0 / (em2piy * em2piy)))
    c0 = (0.5 - Asqrty/ey)
    gpy = gamma + y
    if gpy < 0
        #machine epsilon issues
        gpy = 0.0
    end
    if price <= c0
        return (-sqrt(gpy) + sqrt(gamma - y))
    end
    return (sqrt(gpy) + sqrt(gamma - y))
end

function impliedVolatilitySRGuess(isCall::Bool, price::Real, f::Real, strike::Real, tte::Real, df::Real)::Real
    c,ex = normalizePrice(isCall, price, f, strike, df)
    scaledVol = impliedVolatilitySRGuessUndiscountedCall(c, ex, log(ex))
    return scaledVol / sqrt(tte)
end

function objectiveHouseholderLog(x::T, ex::T, v::TP, logc::TP) where {T,TP}
    v = abs(v)
    h = x / v
    t = v / 2
    h2 = h^2
    t2 = t^2
    sqrt2 = sqrt(T(2))
    Np = erfcx(-(h + t) / sqrt2)
    Nm = erfcx(-(h - t) / sqrt2)
    eh2t2 = exp(-(h2 + t2) / 2)
    norm = 1 / (2 * sqrt(ex)) * eh2t2
    cEstimate = norm * (Np - Nm)
    logcEstimate = log(cEstimate)
    twopi = 2 * T(pi)
    logvega = (2 /sqrt(twopi)) / (Np - Nm)  #u'/u
    volgaOverVega = (h + t) * (h - t) / v
    logvolgaOverVega = volgaOverVega - logvega  #u'' u - u'2 /u^2
    h2mt2 = h2 - t2
    c3OverVega = (-3 * h2 - t2 + h2mt2 * h2mt2) / (v^2)   #
    logc3overVega = c3OverVega - 3 * logvega * volgaOverVega + 2 * logvega * logvega
    return logcEstimate - logc, (logcEstimate - logc) / logvega, logvolgaOverVega, logc3overVega
end

function impliedVolatilitySRHalley(
    isCall::Bool,
    price::TP,
    f::T,
    strike::T,
    tte::T,
    df::T,
    ftolrel::T,
    maxEval::Int,
    # solver::Householder,  # no need to specialise here
)::TP where {T,TP}
    c, ex = normalizePrice(isCall, price, f, strike, df)
    if c >= 1 / ex || c <= 0
        throw(DomainError(c, string("Price out of range, must be < ", 1 / ex, " and > 0. denormalized price=",price)))
    end
    x = log(ex)
    guess = TP(impliedVolatilitySRGuessUndiscountedCall(c, ex,x))
    b = guess
    xtolrel = 32 * eps(T)
    xtol = T(0)
    ftolrel = max(ftolrel, eps(T))
    ftol = ftolrel * max(1, c)
    logc = log(c)
    fb, fbOverfpb, fp2bOverfpb, fp3bOverfpb = objectiveHouseholderLog(x, ex, b, logc)
    if abs(fb) < ftol
        return b / sqrt(tte)
    end
    for iteration = 1:maxEval
        x0 = b
        lf = (fbOverfpb) * (fp2bOverfpb)
        hn = -fbOverfpb
        num = (1 + fp2bOverfpb * hn / 2)
        denom = (1 + fp2bOverfpb * hn + fp3bOverfpb / 6 * hn * hn)
        x1 = x0 + hn * num / denom
        a = x0
        fb, fbOverfpb, fp2bOverfpb, fp3bOverfpb = objectiveHouseholderLog(x, ex, x1, logc)
        b = x1
        xtol_ = xtol + max(1, abs(b)) * xtolrel
        if abs(b - a) <= xtol_ || abs(fb) < ftol
            break
        end
    end
    return b / sqrt(tte)
end
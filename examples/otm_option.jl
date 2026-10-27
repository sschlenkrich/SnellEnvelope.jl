using SnellEnvelope

using DataFrames
using PlotlyBase
using PlotlyJS


# snell envelope integration parameters
# we use a higher number of points and a larger lambda to get
# get reasonable accuracy for far out-of-the-money options
lambda = 10.0
n_points = 1_001

function annuity(rate, years)
    if rate == 0
        return years
    end
    q = 1 / (1.0 + rate)
    annuity_ = q * (1 - q^years) / (1 - q)
    return annuity_
end



function model_inputs(swap_rate, strike_rate, maturity_years, no_call_years, call_put, volatility, vanilla_model_function)
    expiry_years = no_call_years .+ collect(0:(maturity_years-no_call_years-1))
    swap_years = maturity_years .- expiry_years
    #
    swap_rates     = [swap_rate   for _ in expiry_years]
    strike_rates   = [strike_rate for _ in expiry_years]
    #
    nu = volatility .* sqrt.(expiry_years)
    fwd_annuities = [annuity(swap_rate, y) for y in swap_years]
    df_expiry = (1 + swap_rate).^(-1.0.*expiry_years)
    annuities = df_expiry .* fwd_annuities
    #
    european_npvs = annuities .* vanilla_model_function(
        strike_rates,
        swap_rates,
        nu,
        call_put,
    )
    return swap_rates, strike_rates, annuities, european_npvs, call_put
end


rates = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]  # in bp
maturity_years = 10
no_call_years = 1

# collect all results
results = []

# constant lognormal volatility
volatility = 0.35
vanilla_model_function = SnellEnvelope.black_price
for rate in rates
    swap_rate = 0.10
    strike_rate = rate * 1e-4
    call_put = -1.0
    p = model_inputs(swap_rate, strike_rate, maturity_years, no_call_years, call_put, volatility, vanilla_model_function)
    se_normal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.normal, lambda, n_points)
    se_lognormal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.lognormal, lambda, n_points)
    #
    o = SnellEnvelope.black_price(strike_rate, swap_rate, volatility * sqrt(maturity_years), call_put)
    v_b76 = SnellEnvelope.black_implied_volatility(o, strike_rate, swap_rate, 1.0*maturity_years, call_put)
    moneyness = log(strike_rate/swap_rate)
    #
    test_case = "P_$(rate)_L"
    #
    push!(results, (
        test_case=test_case,
        v_b76=v_b76,
        moneyness=moneyness,
        max_european=maximum(p[4]),
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Test case: $test_case, d: $moneyness, E: $(maximum(p[4])), SN: $se_normal, SL: $se_lognormal")
end
println("")
for rate in rates
    strike_rate = 0.10
    swap_rate = rate * 1e-4
    call_put = 1.0
    p = model_inputs(swap_rate, strike_rate, maturity_years, no_call_years, call_put, volatility, vanilla_model_function)
    se_normal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.normal, lambda, n_points)
    se_lognormal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.lognormal, lambda, n_points)
    #
    o = SnellEnvelope.black_price(strike_rate, swap_rate, volatility * sqrt(maturity_years), call_put)
    v_b76 = SnellEnvelope.black_implied_volatility(o, strike_rate, swap_rate, 1.0*maturity_years, call_put)
    moneyness = log(strike_rate/swap_rate)
    #
    test_case = "C_$(rate)_L"
    #
    push!(results, (
        test_case=test_case,
        v_b76=v_b76,
        moneyness=moneyness,
        max_european=maximum(p[4]),
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Test case: $test_case, d: $moneyness, E: $(maximum(p[4])), SN: $se_normal, SL: $se_lognormal")
end
println("")

# constant normal volatility
volatility = 0.0050  # bp
vanilla_model_function = SnellEnvelope.bachelier_price
for rate in rates
    swap_rate = 0.10
    strike_rate = rate * 1e-4
    call_put = -1.0
    p = model_inputs(swap_rate, strike_rate, maturity_years, no_call_years, call_put, volatility, vanilla_model_function)
    se_normal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.normal, lambda, n_points)
    se_lognormal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.lognormal, lambda, n_points)
    #
    o = SnellEnvelope.bachelier_price(strike_rate, swap_rate, volatility * sqrt(maturity_years), call_put)
    v_b76 = SnellEnvelope.black_implied_volatility(o, strike_rate, swap_rate, 1.0*maturity_years, call_put)
    moneyness = log(strike_rate/swap_rate)
    #
    test_case = "P_$(rate)_N"
    #
    push!(results, (
        test_case=test_case,
        v_b76=v_b76,
        moneyness=moneyness,
        max_european=maximum(p[4]),
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Test case: $test_case, d: $moneyness, E: $(maximum(p[4])), SN: $se_normal, SL: $se_lognormal")
end
println("")
for rate in rates
    strike_rate = 0.10
    swap_rate = rate * 1e-4
    call_put = 1.0
    p = model_inputs(swap_rate, strike_rate, maturity_years, no_call_years, call_put, volatility, vanilla_model_function)
    se_normal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.normal, lambda, n_points)
    se_lognormal = SnellEnvelope.swaption_snell_envelope(p..., SnellEnvelope.lognormal, lambda, n_points)
    #
    o = SnellEnvelope.bachelier_price(strike_rate, swap_rate, volatility * sqrt(maturity_years), call_put)
    v_b76 = SnellEnvelope.black_implied_volatility(o, strike_rate, swap_rate, 1.0*maturity_years, call_put)
    moneyness = log(strike_rate/swap_rate)
    #
    test_case = "C_$(rate)_N"
    #
    push!(results, (
        test_case=test_case,
        v_b76=v_b76,
        moneyness=moneyness,
        max_european=maximum(p[4]),
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Test case: $test_case, d: $moneyness, E: $(maximum(p[4])), SN: $se_normal, SL: $se_lognormal")
end
println("")



result_table = DataFrame(results)
# We calculate relative switch values to compare methods
result_table[!, "rel_normal"] = result_table.se_normal ./ result_table.max_european .- 1.0
result_table[!, "rel_lognormal"] = result_table.se_lognormal ./ result_table.max_european .- 1.0

println(result_table)


function plot_relative_switch(result_table::DataFrame, vol_input_type)
    table = filter(row -> row.test_case[end] == vol_input_type, result_table)
    #
    traces = PlotlyBase.AbstractTrace[]
    #
    table_c = filter(row -> row.test_case[begin] == 'C', table)
    table_p = filter(row -> row.test_case[begin] == 'P', table)
    #
    trace = scatter(
        name = "P - Normal Approx.",
        x = table_p.moneyness,
        y = table_p.rel_normal * 100,
        mode = "markers",
        marker_color = "green",
        marker_line_color = "green",
        marker_line_width = 1.0,
        marker_symbol = "x-thin",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "P - Lognormal Approx.",
        x = table_p.moneyness,
        y = table_p.rel_lognormal * 100,
        mode = "markers",
        marker_color = "green",
        marker_line_color = "green",
        marker_line_width = 1.0,
        marker_symbol = "cross-thin",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "C - Normal Approx.",
        x = table_c.moneyness,
        y = table_c.rel_normal * 100,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "x-thin",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "C - Lognormal Approx.",
        x = table_c.moneyness,
        y = table_c.rel_lognormal * 100,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "cross-thin",
    )
    push!(traces, trace)
    #
    if vol_input_type == 'L'
        title = "Relative switch value for constant lognormal volatility input"
    elseif vol_input_type == 'N'
        title = "Relative switch value for constant normal volatility input"
    else
        error("Invalid vol_input_type: $vol_input_type. Must be 'L' or 'N'.")
    end
    layout = Layout(
        title = title,
        xaxis_title="Moneyness log(K/F)",
        yaxis_title="Relative switch value (%)",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end


p = plot_relative_switch(result_table, 'L')
savefig(p, "./data/otm_options_lognormal.svg")

p = plot_relative_switch(result_table, 'N')
savefig(p, "./data/otm_options_normal.svg")

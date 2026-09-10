using SnellEnvelope

using CSV
using DataFrames
using PlotlyBase
using PlotlyJS

file_name = "./data/callable_floater_analysis.csv"
table  = CSV.read(file_name, DataFrame)
table = filter(row -> row.european_npv > 0.0, table)  # 0.0 causes error with implied vol

# snell envelope integration parameters
lambda = 5.0
n_points = 101

println(table)

results = []
for strike_rate in unique(table.strike_rate)
    table_tc = filter(row -> row.strike_rate == strike_rate, table)
    #
    se_normal = SnellEnvelope.swaption_snell_envelope(
        table_tc.fair_rate,
        table_tc.strike_rate,
        table_tc.annuity,
        table_tc.european_npv,
        float(table_tc.payer_or_receiver[begin]),
        #
        SnellEnvelope.normal,
        lambda,
        n_points,
    )
    se_lognormal = SnellEnvelope.swaption_snell_envelope(
        table_tc.fair_rate,
        table_tc.strike_rate,
        table_tc.annuity,
        table_tc.european_npv,
        float(table_tc.payer_or_receiver[begin]),
        #
        SnellEnvelope.lognormal,
        lambda,
        n_points,
    )
    #
    push!(results, (
        strike_rate=strike_rate,
        max_european=maximum(table_tc.european_npv),
        bermudan=table_tc.bermudan_npv[begin],
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Strike rate: $strike_rate, E: $(maximum(table_tc.european_npv)), B: $(table_tc.bermudan_npv[begin]), SN: $se_normal, SL: $se_lognormal")
end

result_table = DataFrame(results)

# We calculate relative switch values to compare methods
result_table[!, "rel_switch"] = result_table.bermudan ./ result_table.max_european .- 1.0
result_table[!, "rel_normal"] = result_table.se_normal ./ result_table.max_european .- 1.0
result_table[!, "rel_lognormal"] = result_table.se_lognormal ./ result_table.max_european .- 1.0

println(result_table)


# volatility smile

function safe_black_implied_volatility(price, strike, forward, T, call_put)
    try
        return SnellEnvelope.black_implied_volatility(price, strike, forward, T, call_put)
    catch e
        @warn "Failed to compute implied volatility: $e"
        return NaN
    end
end


function plot_smiles(expiries, table)
    traces = PlotlyBase.AbstractTrace[]
    for expiry in expiries
        table_expiry = filter(row -> row.expiry == expiry, table)
        implied_vols = [
            safe_black_implied_volatility(price, strike, forward, T, call_put)
            for (price, strike, forward, T, call_put) in zip(
                table_expiry.european_npv ./ table_expiry.annuity,
                table_expiry.strike_rate,
                table_expiry.fair_rate,
                table_expiry.exercise_time,
                table_expiry.payer_or_receiver,
            )
        ]
        trace = scatter(
            name = "T = " * expiry,
            x = table_expiry.strike_rate * 100,
            y = implied_vols * 100,
            mode = "lines",
            line_shape = "spline",
            line_smoothing = 1.3,
        )
        push!(traces, trace)
    end
    title = "Implied credit spread volatility"
    layout = Layout(
        title = title,
        xaxis_title="Funding spread rate (%)",
        yaxis_title="Black volatility (%)",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end


# volatility term structure

function plot_term_structures(strikes, table)
    traces = PlotlyBase.AbstractTrace[]
    for strike in strikes
        table_strike = filter(row -> row.strike_rate == strike, table)
        implied_vols = [
            safe_black_implied_volatility(price, strike, forward, T, call_put)
            for (price, strike, forward, T, call_put) in zip(
                table_strike.european_npv ./ table_strike.annuity,
                table_strike.strike_rate,
                table_strike.fair_rate,
                table_strike.exercise_time,
                table_strike.payer_or_receiver,
            )
        ]
        trace = scatter(
            name = "K = " * string(strike * 100) * "%",
            x = table_strike.exercise_time,
            y = implied_vols * 100,
            mode = "lines",
            line_shape = "spline",
            line_smoothing = 1.3,
        )
        push!(traces, trace)
    end
    title = "Implied credit spread volatility"
    layout = Layout(
        title = title,
        xaxis_title="Expiry",
        yaxis_title="Black volatility (%)",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end


# Option prices

function plot_prices(result_table::DataFrame)
    traces = PlotlyBase.AbstractTrace[]
    #
    trace = scatter(
        name = "AMC Bermudan",
        x = result_table.strike_rate * 100,
        y = result_table.bermudan,
        mode = "markers",
        marker_color = "red",
        marker_symbol = "diamond",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "MC Max-European",
        x = result_table.strike_rate * 100,
        y = result_table.max_european,
        mode = "markers",
        marker_color = "red",
        marker_symbol = "arrow-up",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "Normal Approx.",
        x = result_table.strike_rate * 100,
        y = result_table.se_normal,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "x-thin",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "Lognormal Approx.",
        x = result_table.strike_rate * 100,
        y = result_table.se_lognormal,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "cross-thin",
    )
    push!(traces, trace)
    #
    title = "Option prices (AMC vs. Normal and Lognormal Approximations)"
    layout = Layout(
        title = title,
        xaxis_title="Funding spread rate (%)",
        yaxis_title="Option price",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end

# Switch values

function plot_relative_switch(result_table::DataFrame)
    traces = PlotlyBase.AbstractTrace[]
    #
    trace = scatter(
        name = "AMC",
        x = result_table.strike_rate * 100,
        y = result_table.rel_switch * 100,
        mode = "markers",
        marker_color = "red",
        marker_symbol = "diamond",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "Normal Approx.",
        x = result_table.strike_rate * 100,
        y = result_table.rel_normal * 100,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "x-thin",
    )
    push!(traces, trace)
    #
    trace = scatter(
        name = "Lognormal Approx.",
        x = result_table.strike_rate * 100,
        y = result_table.rel_lognormal * 100,
        mode = "markers",
        marker_color = "blue",
        marker_line_color = "blue",
        marker_line_width = 1.0,
        marker_symbol = "cross-thin",
    )
    push!(traces, trace)
    #
    title = "Relative switch value (AMC vs. Normal and Lognormal Approximations)"
    layout = Layout(
        title = title,
        xaxis_title="Funding spread rate (%)",
        yaxis_title="Relative switch value (%)",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end


p = plot_smiles(["2y", "5y", "9y"], table)
savefig(p, "./data/callable_floater_smiles.svg")

p = plot_prices(result_table)
savefig(p, "./data/callable_floater_option_prices.svg")

p = plot_relative_switch(result_table)
savefig(p, "./data/callable_floater_relative_switch_values.svg")

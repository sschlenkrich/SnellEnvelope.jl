using SnellEnvelope

using CSV
using DataFrames
using PlotlyBase
using PlotlyJS

file_name = "./data/bermudan_versus_european_analysis.csv"
table  = CSV.read(file_name, DataFrame)

# snell envelope integration parameters
lambda = 5.0
n_points = 101

results = []
for tc in unique(table.test_case)
    table_tc = filter(row -> row.test_case == tc, table)
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
        test_case=tc,
        max_european=maximum(table_tc.european_npv),
        bermudan=table_tc.bermudan_npv[begin],
        se_normal=se_normal,
        se_lognormal=se_lognormal,
        )
    )
    println("Test case: $tc, E: $(maximum(table_tc.european_npv)), B: $(table_tc.bermudan_npv[begin]), SN: $se_normal, SL: $se_lognormal")
end

result_table = DataFrame(results)

# We calculate relative switch values to compare methods
result_table[!, "rel_switch"] = result_table.bermudan ./ result_table.max_european .- 1.0
result_table[!, "rel_normal"] = result_table.se_normal ./ result_table.max_european .- 1.0
result_table[!, "rel_lognormal"] = result_table.se_lognormal ./ result_table.max_european .- 1.0


function plot_test_case(tc::String, input_table::DataFrame, result_table::DataFrame)
    t = input_table[input_table.test_case .== tc, :]
    #
    option_type = (tc[1] == 'P') ? "Pay" : "Rec"
    strike = string(t.strike_rate[1]*100) * "%"
    interest_rate = string(parse(Float64, tc[3:4])) * "%"
    mean_reversion = string(parse(Float64, tc[7:8])) * "%"
    title = "$option_type $strike, at interest rate $interest_rate, mean reversion $mean_reversion"
    #
    labels = t[!, "expiry"] .* "-" .* t[!, "swap_term"]
    trace1 = bar(x=labels, y=t.european_npv)
    #
    t = result_table[result_table.test_case .== tc, :]
    labels = ["Bermudan", "Normal Approx.", "Lognormal Approx."]
    values = [t.bermudan[1], t.se_normal[1], t.se_lognormal[1]]
    trace2 = bar(x=labels, y=values)
    layout = Layout(
        title = title,
        xaxis_title="Option",
        yaxis_title="Value (per 10k notional)",
        width = 800,
        height = 500,
        showlegend = false,
    )
    return plot([trace1, trace2], layout)
end


function plot_relative_switch(payer_or_receiver::Char, result_table::DataFrame)
    idx_list = [ s[1] == payer_or_receiver for s in result_table.test_case ]
    t = result_table[idx_list, :]
    option_type = (payer_or_receiver == 'P') ? "Pay" : "Rec"
    strike = "5%"
    mean_reversion_labels = ["10", "05", "01"]
    mean_reversion_symbols = ["arrow-up", "diamond", "arrow-down"]
    traces = PlotlyBase.AbstractTrace[]
    for (mr_label, mr_symbol) in zip(mean_reversion_labels, mean_reversion_symbols)
        mean_reversion = string(parse(Float64, mr_label)) * "%"
        idx_list = [ s[7:8] == mr_label for s in t.test_case ]
        t_mr = t[idx_list, :]
        x_values = [ parse(Float64, s[3:4]) for s in t_mr.test_case ]  # flat rates in %
        y_values = t_mr.rel_switch
        trace = scatter(
            name = "a =" * mean_reversion,
            x = x_values,
            y = t_mr.rel_switch * 100,
            mode = "markers",
            marker_color = "red",
            marker_symbol = mr_symbol,
        )
        push!(traces, trace)
        if mr_label == mean_reversion_labels[end]
            trace = scatter(
                name = "Normal Approx.",
                x = x_values,
                y = t_mr.rel_normal * 100,
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
                x = x_values,
                y = t_mr.rel_lognormal * 100,
                mode = "markers",
                marker_color = "blue",
                marker_line_color = "blue",
                marker_line_width = 1.0,
                marker_symbol = "cross-thin",
            )
            push!(traces, trace)
        end
    end
    title = "Relative switch value, " * option_type * " " * strike
    layout = Layout(
        title = title,
        xaxis_title="Flat interest rate (%)",
        yaxis_title="Relative switch value (%)",
        width = 800,
        height = 500,
        showlegend = true,
    )
    return plot(traces, layout)
end

p = plot_test_case("P_03_a05", table, result_table)
savefig(p, "./data/bermudan_versus_european_analysis_P_03_a05.svg")

p = plot_test_case("P_07_a05", table, result_table)
savefig(p, "./data/bermudan_versus_european_analysis_P_07_a05.svg")

p = plot_relative_switch('P', result_table)
savefig(p, "./data/relative_switch_value_Payer.svg")

p = plot_relative_switch('R', result_table)
savefig(p, "./data/relative_switch_value_Receiver.svg")
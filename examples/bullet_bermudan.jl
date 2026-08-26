using SnellEnvelope

using CSV
using DataFrames

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


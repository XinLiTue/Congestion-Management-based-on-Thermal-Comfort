

# compute ray f(x) through two adjacent median PPD points
f(x, (x1, y1), (x2, y2)) = y1 + (y2 - y1) / (x2 - x1) * (x - x1)

function add_ppd(model::Model, grid::Grid, df::DataFrame)
    # @objective(model, Min, J_gen)
    Ta = df.Ta
    PPD = df.Q2
    # medians = [(Ta[i], PPD[i]) for i in eachindex(Ta)]
    medians = [(Ta[i], PPD[i]) for i in [1, 5, 8, 11, 15]]

    # sets
    sets = get_sets(grid)
    H_HP, T = sets.H_HP, sets.T

    # we need an auxilary variable to represent the PPD
    @variable(model, PPD[i in H_HP, t in T])

    # T_i can be indexed by [H_HP, T]
    T_i = model[:Te][:, :, :i]

    # prepare linear functions for outer approximation, i.e. f(x) = a + b * x
    # PPD ≥ a_i * ​x + b_i    ∀i=1,…,n, where n is the number of medians
    @constraint(model,
        PPDOuterApprox[t in T, i in H_HP, m in 1:length(medians)-1],
        PPD[i, t] ≥ f(T_i[i, t], medians[m], medians[m+1])
    )
end
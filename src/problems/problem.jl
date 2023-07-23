"""
Problem type and associated methods
"""


abstract type AbstractCollocationProblem end

struct CollocationProblem <: AbstractCollocationProblem
    model::Model
    hs_defect::Function
    get_target_state::Union{Function,Nothing}

    function CollocationProblem(model, hs_defect, get_target_state=nothing)
        new(model, hs_defect, get_target_state)
    end
end



"""
    solve!(problem::CollocationProblem)

Solve a collocation problem
"""
function solve!(problem::CollocationProblem)
    # solve the problem
    optimize!(problem.model)
    # flush the memory
    empty!(memoize_cache(problem.hs_defect))
    if isnothing(problem.get_target_state) == false
        empty!(memoize_cache(problem.get_target_state))
    end
    return
end

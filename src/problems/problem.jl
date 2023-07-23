"""
Problem type and associated methods
"""


abstract type AbstractCollocationProblem end

struct CollocationProblem <: AbstractCollocationProblem
    model::Model
    hs_defect::Function

    function CollocationProblem(model, hs_defect)
        new(model, hs_defect)
    end
end



"""
Solve a collocation problem
"""
function solve!(problem::CollocationProblem)
    # solve the problem
    optimize!(problem.model)
    # flush the memory
    memoize_cache(problem.hs_defect)
    return
end
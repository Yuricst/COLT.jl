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
Construct Hermite-Simpson collocation JuMP problem from partially initialized model.
"""
struct HermiteSimpsonProblem <: AbstractCollocationProblem
    model::Model
    hs_defect::Function
    get_target_state::Union{Function,Nothing}

    function HermiteSimpsonProblem(
        model,
        dynamics_orbitraise::Function,
        state_variables::Vector,
        control_variables::Vector,
        ode_parameters::Union{Vector,Nothing},
        get_target_state::Union{Function,Nothing} = nothing,
    )
        @macroexpand @assert length(state_variables[1]) == length(control_variables[1]) "state and control must have same length!"
        
        # get dimensions of state and control variables
        nx, nu = length(state_variables), length(control_variables)

        # get dimension of nodes
        N = length(state_variables[1]) - 1

        # safety assertions
        for i in range(1,nx)
            @macroexpand @assert length(state_variables[i]) - 1 == N "state variables must have same length!"
        end
        for i in range(1,nu)
            @macroexpand @assert length(control_variables[i]) - 1 == N "control variables must have same length!"
        end

        # get function for defects
        hs_defect = get_hs_defect_function(
            dynamics_orbitraise,
            nx,
            nu,
            ode_parameters,
        )

        # NEW API FOR NONLINEAR PROGRAM (JuMP v1.x)
        txu0_txu1_uc_tf = []
        for i = 1:N
            ix1, ix2 = i, i+1
            iu1, iu2 = i, i+1
            push!(
                txu0_txu1_uc_tf,
                [
                    model[:t][ix1], [_x[ix1] for _x in state_variables]..., [_u[iu1] for _u in control_variables]...,
                    model[:t][ix2], [_x[ix2] for _x in state_variables]..., [_u[iu2] for _u in control_variables]...,
                    model[:tf]
                ]
            )
        end

        defects = @expression(model, [i = 1:N], hs_defect(txu0_txu1_uc_tf[i]))
        @constraint(model, [i = 1:N], defects[i] .== 0.0)

        new(model, hs_defect, get_target_state)
    end
end


"""
    solve!(problem::CollocationProblem)

Solve a collocation problem
"""
function solve!(problem::AbstractCollocationProblem)
    # solve the problem
    optimize!(problem.model)
    # flush the memory
    try
        empty!(memoize_cache(problem.hs_defect))
    catch
        # do nothing
    end
    if isnothing(problem.get_target_state) == false
        empty!(memoize_cache(problem.get_target_state))
    end
    return
end

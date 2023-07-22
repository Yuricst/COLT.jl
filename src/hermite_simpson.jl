"""
Functions for Hermite-Simpson collocation methods
"""

function get_hs_defect_function(
    state_variables,
    control_variables,
    dynamics::Function,
    params,
)
    # number of states and controls 
    nx, nu = length(state_variables), length(control_variables)

    # prepare function that creates Hermite-Simpson defect constraints
    @memoize function hs_defect(txu0_txu1_uc...)
        # unsplat the input
        txu0 = txu0_txu1_uc[1:1+nx+nu]
        txu1 = txu0_txu1_uc[2+nx+nu:2*(1+nx+nu)]
        tc = (txu0[1] + txu1[1])/2
        
        T = txu1[1] - txu0[1]
        
        # compute state & state derivative at center
        x0 = collect(txu0[2:2+nx-1])
        x1 = collect(txu1[2:2+nx-1])
        u0 = collect(txu0[2+nx:end])
        u1 = collect(txu1[2+nx:end])
        uc = (u0 + u1)/2
        f0 = dynamics(txu0[1], x0, u0, params)
        f1 = dynamics(txu1[1], x1, u1, params)
        
        xc = (x0 + x1)/2 + T*(f0 + f1)/8
        
        # compute dynamics at center
        fc = dynamics(tc, xc, uc, params)
        return (x0 - x1) + T/6*(f0 + 4fc + f1)
    end
    return hs_defect
end
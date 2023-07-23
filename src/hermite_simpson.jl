"""
Functions for Hermite-Simpson collocation methods
"""

"""
Get the Hermite-Simpson defect function for a given dynamics function
"""
function get_hs_defect_function(
    dynamics::Function,
    nx::Int,
    nu::Int,
    params,
)
    # prepare function that creates Hermite-Simpson defect constraints
    @memoize function hs_defect(txu0_txu1_uc_tf...)
        # unsplat the input
        txu0 = txu0_txu1_uc_tf[1:1+nx+nu]
        txu1 = txu0_txu1_uc_tf[2+nx+nu:2*(1+nx+nu)]
        tf = txu0_txu1_uc_tf[2*(1+nx+nu)+1]            # final time 
        
        # recreate dimensional time
        t0 = txu0[1]*tf
        t1 = txu1[1]*tf
        tc = (t0 + t1)/2
        T = t1 - t0
        
        # compute state & state derivative at center
        x0 = collect(txu0[2:2+nx-1])
        x1 = collect(txu1[2:2+nx-1])
        u0 = collect(txu0[2+nx:end])
        u1 = collect(txu1[2+nx:end])
        uc = (u0 + u1)/2
        f0 = dynamics(t0, x0, u0, params)
        f1 = dynamics(t1, x1, u1, params)
        
        xc = (x0 + x1)/2 + T*(f0 + f1)/8
        
        # compute dynamics at center
        fc = dynamics(tc, xc, uc, params)
        return (x0 - x1) + T/6*(f0 + 4fc + f1)
    end
    return hs_defect
end
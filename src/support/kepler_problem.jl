"""Functions for solving Kepler's problem

Taken from AstrodynamicsBase.jl (https://github.com/Yuricst/AstrodynamicsBase.jl)
"""

"""
Compute semi-major axis from Cartesian state
"""
function get_semiMajorAxis(state::Array{<:Real,1}, mu::Real)
	# Initialize Cartesian Polistion and Velocity
    r = state[1:3]
    v = state[4:6]
    a     = 1.0/(2.0/norm(r) - norm(v)^2/mu)     # Semi-major axis
	return a
end


"""
Compute eccentricity from Cartesian state
"""
function get_eccentricity(state::Array{<:Real,1}, mu::Real)
	# Initialize Cartesian Polistion and Velocity
    r = state[1:3]
    v = state[4:6]
    h = cross(r, v)  # Angular momentum vector
    p = norm(h)*norm(h)/mu                         # Semi-latus rectum
    a = 1.0/(2.0/norm(r) - norm(v)^2/mu)     # Semi-major axis
    # Numerical stability hack for circular and near-circular orbits
    # Ensures that (1-p/a) is always positive
    if isapprox(a, p, atol=1e-9, rtol=1e-8)
        p = a
    end
    e = sqrt(1 - p/a)    # Eccentricity
	return e
end


"""
    hypertrig_s(z::Union{Real,ForwardDiff.Dual,QuadExpr}, epsilon::Real=1e-6)

Hyperbolic sine function
"""
function hypertrig_s(z::Union{Real,ForwardDiff.Dual,QuadExpr}, epsilon::Real=1e-6)
    # if z > epsilon
    #     s = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3
    # elseif z < -epsilon
    #     s = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3
    # else
    #     s = 1/6 - z/120 + z^2/5040 - z^3/362880
    # end
    s = 1/6 - z/120 + z^2/5040 - z^3/362880
    return s
end



"""
    hypertrig_c(z::Union{Real,ForwardDiff.Dual,QuadExpr}, epsilon::Real=1e-6)

Hyperbolic cosine function
"""
function hypertrig_c(z::Union{Real,ForwardDiff.Dual,QuadExpr}, epsilon::Real=1e-6)
    # if z > epsilon
    #     c = (1.0 - cos(sqrt(z))) / z
    # elseif z < -epsilon
    #     c = (cosh(-z) - 1) / (-z)
    # else
    #     c = 1/2 - z/24 + z^2/720 - z^3/40320
    # end
    c = 1/2 - z/24 + z^2/720 - z^3/40320
    return c
end


"""Lagrange parameter functions"""
function universal_functions(x::Union{Real,ForwardDiff.Dual,AffExpr}, alpha::Float64)
    """Function computes U0 ~ U3 from G. Der 1996"""
    # evaluate hypertrig function
    S = hypertrig_s(alpha * x^2)
    C = hypertrig_c(alpha * x^2)
    # parameters
    u0 = 1 - alpha * x^2 * C
    u1 = x * (1 - alpha * x^2 * S)
    u2 = x^2 * C
    u3 = x^3 * S
    return u0, u1, u2, u3
end


"""Function computes Kepler's time equation and its derivatives in G. Der 1996 form"""
function kepler_der_residual(
    x, #::Union{Real,ForwardDiff.Dual},
    alpha::Float64,
    t, #::Union{Real,ForwardDiff.Dual},
    t0::Union{Real,ForwardDiff.Dual},
    mu::Float64,
    r0::Float64,
    sigma0::Float64,
)
    u0, u1, u2, u3 = universal_functions(x, alpha)
    Fun = r0 * u1 + sigma0 * u2 + u3 - sqrt(mu) * (t - t0)
    dF = r0 * u0 + sigma0 * u1 + u2
    d2F = sigma0 * u0 + (1 - alpha * r0) * u1
    return Fun, dF, d2F
end


"""Function computes Lagrange coefficients (as functionined in G. Der 1996)"""
function lagrange_coefficients_der(
    mu::Real,
    alpha::Real,
    r0::Real,
    v0::Real,
    sigma0::Union{Real,ForwardDiff.Dual},
    u0::Union{Real,ForwardDiff.Dual},
    u1::Union{Real,ForwardDiff.Dual},
    u2::Union{Real,ForwardDiff.Dual},
    u3::Union{Real,ForwardDiff.Dual},
    r::Union{Real,ForwardDiff.Dual},
    sigma::Union{Real,ForwardDiff.Dual},
)
    # scalar function
    f = 1.0 - u2 / r0
    g = r0 * u1 / sqrt(mu) + sigma0 * u2 / sqrt(mu)
    fdot = -sqrt(mu) / (r * r0) * u1
    gdot = 1.0 - u2 / r
    return f, g, fdot, gdot
end



"""
    keplerder_nostm_colt(mu::Float64, state0::Vector{Float64}, t0::Float64, t::Float64, tol::Float64, maxiter::Int)

Function computes position at some future time t, without computing STM
Formulation from G. Der formulation.

# Arguments
    `mu::Float64`: gravitational parameter
    `state0::Vector{Float64}`: initial state [x,y,z,vx,vy,vz]
    `t0::Float64`: initial time at state0
    `t::Float64`: final time
    `tol::Float64`: tolerance on Laguerre-correction, suggest 1.e-14
    `maxiter::Int`: max allowed iteration allowed for Laguerre-correction, suggest 10

# Returns
    `(array)`: final state
"""
function keplerder_nostm_colt(
    mu::Real,
    state0::Array{<:Real,1},
    t0::Real,
    t,
    tol::Float64,
    maxiter::Int,
)
    # ------------------------------------------ #
    # SET-UP PROBLEM
    r0, v0 = state0[1:3], state0[4:6]
    sma = get_semiMajorAxis(state0, mu)
    alpha = 1.0 / sma
    sigma0 = dot(r0, v0) / sqrt(mu)

    # ------------------------------------------ #
    # ITERATION WITH LAGUERRE-CONWAY METHOD
    # initial guess based on eccentricity
    ecc = norm(get_eccentricity(state0, mu))
    # initial guess for circular or elliptical case
    x0 = alpha * sqrt(mu) * (t - t0)
    # if ecc < 1
    #     # initial guess for circular or elliptical case
    #     x0 = alpha * sqrt(mu) * (t - t0)
    # else
    #     # initial guess for parabola or hyperbola
    #     x0 = sqrt(mu) * (t - t0) / (10 * norm(r0))
    # end
    
    # initialize final guess
    x1 = x0
    # initialize iteratation
    count = 0
    Fun = 1.0  # FIXME - need storage
    while count < maxiter
        # evaluate function
        Fun, dF, d2F = kepler_der_residual(x0, alpha, t, t0, mu, norm(r0), sigma0)
        if abs(Fun) < tol
            break
        end
        # Laguerre-Conway iteration
        x1 = x0 - 5 * Fun / (dF + dF / abs(dF) * sqrt(abs(16 * dF^2 - 20 * Fun * d2F)))
        # else
        #     # Newton iteration
        #     x1 = x0 - Fun / dF
        # end
        count += 1
        x0 = x1
    end

    # ------------------------------------------ #
    # COMPUTE FINAL POSITION
    # convert back to position
    u0, u1, u2, u3 = universal_functions(x1, alpha)
    r_scal = norm(r0) * u0 + sigma0 * u1 + u2
    sigma = sigma0 * u0 + (1 - alpha * norm(r0)) * u1

    # get lagrange coefficients
    f, g, fdot, gdot = lagrange_coefficients_der(
        mu,
        alpha,
        norm(r0),
        norm(v0),
        sigma0,
        u0,
        u1,
        u2,
        u3,
        r_scal,
        sigma,
    )
    # create map for state   # --- FIXME! can probably speed this up!
    rmap = hcat(f * Matrix(I, 3, 3), g * Matrix(I, 3, 3))
    vmap = hcat(fdot * Matrix(I, 3, 3), gdot * Matrix(I, 3, 3))
    fullmap = vcat(rmap, vmap)
    state1 = fullmap * state0
    return state1
end
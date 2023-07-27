module COLT

    using JuMP
    using Ipopt 
    using Memoize
    using LRUCache
    using LinearAlgebra
    import ForwardDiff

    include("hermite_simpson.jl")
    include("support/kepler_problem.jl")
    include("problems/problem.jl")
    include("problems/orbit_raising.jl")
    include("problems/twobody_rendez_vous.jl")

end # module COLT

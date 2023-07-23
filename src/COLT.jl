module COLT

    using JuMP
    using Ipopt 
    using Memoize

    include("hermite_simpson.jl")
    include("problems/problem.jl")
    include("problems/orbit_raising.jl")

end # module COLT

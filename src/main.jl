module SlaterKosterTB
include("structure.jl")
include("klib.jl")
include("hamilton.jl")
include("param.jl")
include("draw.jl")
using .Structure, .Klib, .Params, .Hamiltonian, .Draw
end 
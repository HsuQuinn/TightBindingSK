module test4
# https://zhuanlan.zhihu.com/p/565526694
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw
using Plots

# structure args
lat = [[1.0, 0.0, 0.0], [-0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 10.0]]
a = Structure.CreateAtom("a", [2/3,1/3,0], ["dyz","dxz"])
b = Structure.CreateAtom("b", [1/3,2/3,0], ["dyz","dxz"])
atoms = [a,b]
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0,"e_p" => 0),
    "ii" => Dict("V_dds"=>  0.0,
                "V_ddp"=>   0.4096,
                "V_ddd"=>   0.0259)
)
cutoff = 0.8
# kpath args
Γ = [0.0, 0.0, 0.0]
M = [0.5, 0.0, 0.0]
K = [1/3, 1/3, 0.0]
kptnum_perline = 100
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Γ, M, K, Γ],
    "labels"    =>  ["Γ", "M", "K", "Γ"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)
# diag and plot
evals = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict,evals,"honey")
end
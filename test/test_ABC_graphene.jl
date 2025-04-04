module test6
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
lat = [[1.0, 0.0, 0.0], [0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 50.0]]
atoms = []
for i in 1:50
    a = Structure.CreateAtom("C", [i/3.0, i/3.0, 0.01 * i], ["pz"])
    b = Structure.CreateAtom("C", [(i+1)/3.0, (i+1)/3.0, 0.01 * i], ["pz"])
    push!(atoms, a)
    push!(atoms, b)
end
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_p" => 0),
    "ii" => Dict("V_ppp" => 0.2,"V_pps" => 0.18)
)
cutoff = 0.7
# kpath args
Γ = [0.0, 0.0, 0.0]
M = [0.5, 0.5, 0.0]
K = [2/3, 1/3, 0.0]
kptnum_perline = 15
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Γ, K, M, Γ],
    "labels"    =>  ["Γ", "K", "M", "Γ"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)
# diag and plot
H = Hamiltonian.gen_ham_k(Γ, lattice, atoms, interactions, cutoff)
@show H[2,3]
@show H[1,3]
evals,evecs = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict, evals,"ABC_graphene")
# Draw.plot_contour(kdict, lattice,atoms,interactions,cutoff, [-1,1], [-1,1],"graphene", index=2)
end
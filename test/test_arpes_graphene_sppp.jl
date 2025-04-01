module test0
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
# Note: CreateAtom.pos should be the relative coords!!! different from chinook!

alatt = 1.0
lat = [[-0.5*alatt, sqrt(3/4.)*alatt, 0.0], [0.5*alatt, sqrt(3/4.)*alatt, 0.0], [0.0, 0.0, 100.0]]
a = Structure.CreateAtom("C", [0,0,0], ["s","py","px","pz"])
b = Structure.CreateAtom("C", [1/3,1/3,0], ["s","py","px","pz"])
atoms = [a,b]
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => -8.81,"e_p" => -0.44),
    "ii" => Dict("V_sss"=>  -5.279,
                "V_sps"=>   5.618,
                "V_pps"=>   6.05,
                "V_ppp"=>   -3.07)
)
cutoff = alatt * 0.7
# kpath args
G = [0.0, 0.0, 0.0]
K = [1.0/3,2.0/3,0]
M = [0,0.5,0.0]
kptnum_perline = 100
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [G,K,M,G],
    "labels"    =>  ["G","K","M","G"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)
# diag and plot
@show H = Hamiltonian.gen_ham_k(G, lattice, atoms, interactions, cutoff)
evals,evecs = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict,evals,"arpes_graphene")
Draw.plot_contour(kdict, lattice,atoms,interactions,cutoff, [-1,1], [-1,1],"arpes_graphene", index=1)
opdict = Dict(
    "orbital"   =>  "pz" 
)
Draw.plot_orbital_projection(kdict,atoms,evals,evecs,opdict)
end
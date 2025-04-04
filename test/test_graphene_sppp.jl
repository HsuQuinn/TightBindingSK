module test2
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
lat = [[1.0, 0.0, 0.0], [0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 100.0]] 
a = Structure.CreateAtom("C", [1.0/3.0,1.0/3.0,0], ["s","px","py","pz"])
b = Structure.CreateAtom("C", [2.0/3.0,2.0/3.0,0], ["s","px","py","pz"])
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
cutoff = 0.8 * lattice_constant 
# kpath args
Γ = [0.0, 0.0, 0.0]
M = [0.5, 0.5, 0.0]
K = [2/3, 1/3, 0.0]
kptnum_perline = 80
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Γ, K, M, Γ],
    "labels"    =>  ["Γ", "K", "M", "Γ"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)
# diag and plot
evals, evecs = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict,evals,"graphene_4orb")
# orbital projection
opdict = Dict(
    "orbital"   =>  "pz" 
)
Draw.plot_orbital_projection(kdict,atoms,evals,evecs,opdict)
end
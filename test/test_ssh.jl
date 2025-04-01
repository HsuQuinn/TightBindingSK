module test0
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw
using Plots

# structure args
lat = [[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
a = Structure.CreateAtom("a", [0,0,0], ["s","px"])
atoms = [a]
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0,"e_p" => 0),
    "ii" => Dict("V_sss"=>  -0.2,
                "V_sps"=>   -0.05,
                "V_pps"=>   0.2,
                "V_ppp"=>   0.0)
)
cutoff = 1.2
# kpath args
Γ = [0.0, 0.0, 0.0]
Mx = [0.5, 0.0, 0.0]
My = [-0.5, 0.0, 0.0]
kptnum_perline = 100
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Mx, Γ, My],
    "labels"    =>  ["Mx", "Γ", "My"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)

# diag and plot
evals = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict,evals,"ssh")
end
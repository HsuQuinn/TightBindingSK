module test3
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
# Note: CreateAtom.pos should be the relative coords!!! different from chinook!

alatt = 5.0
lat = [[sqrt(0.5)*alatt, sqrt(0.5)*alatt, 0.0], [sqrt(0.5)*alatt, -sqrt(0.5)*alatt, 0.0], [0.0, 0.0, 5.0]]
a = Structure.CreateAtom("Sb", [0,0,0], ["py","px","pz"])
b = Structure.CreateAtom("Sb", [1/2,1/2,0], ["py","px","pz"])
atoms = [a,b]
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0.0,"e_p" => 0.7),
    "ii" => Dict("V_sss"=>  0.0,
                "V_sps"=>   0.0,
                "V_pps"=>   0.25,
                "V_ppp"=>   -1.0)
)
cutoff = alatt * 0.72
# Ep = 0.7
# Vpps = 0.25
# Vppp = -1.0
# VSK = {'051':Ep,'005511S':Vpps,'005511P':Vppp}
# cutoff = 0.72*a
# kpath args
G = [0.0, 0.0, 0.0]
Mx = [0.5, 0.5, 0.0]
My = [0.5, -0.5, 0.0]
kptnum_perline = 100
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Mx,G,My],
    "labels"    =>  ["Mx","G","My"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)
# diag and plot
@show H = Hamiltonian.gen_ham_k(G, lattice, atoms, interactions, cutoff)
evals,evecs = Hamiltonian.SolveHk(kdict, lattice, atoms, interactions, cutoff)
Draw.plot_bandstructre(kdict,evals,"chinook")
Draw.plot_contour(kdict, lattice,atoms,interactions,cutoff, [-1,1], [-1,1],"chinook", index=1)
opdict = Dict(
    "orbital"   =>  "pz" 
)
Draw.plot_orbital_projection(kdict,atoms,evals,evecs,opdict)
end
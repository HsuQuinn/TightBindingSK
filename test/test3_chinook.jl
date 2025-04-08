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
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
a = Structure.CreateAtom("Sb", [0,0,0], ["pz","px","py"])
b = Structure.CreateAtom("Sb", [1/2,1/2,0], ["py","px","pz"])
atoms = [a,b]
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0.0,"e_p" => 0.7),
    "ii" => Dict("V_pps"=>   0.25,
                "V_ppp"=>   -1.0),
    "cutoff" => alatt * 0.72
)
spin = Dict(
    "soc"  =>  true,
    "λ"    =>  Dict("Sb"  =>  0.5),
)
rdict = Dict(
    "lattice"       =>  lattice,
    "atoms"         =>  atoms,
    "interactions"  =>  interactions,
    "spin"          =>  spin
)
# kpath args
G = [0.0, 0.0, 0.0]
Mx = [0.5, 0.5, 0.0]
My = [0.5, -0.5, 0.0]
kptnum_perline = 100
kdict = Dict(
    "pts"       =>  [Mx,G,My],
    "labels"    =>  ["Mx","G","My"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  Structure.ReciprocalLattice(lattice)
)

# diag and plot
basis = Hamiltonian.gen_basis(rdict) 
# @show H = Hamiltonian.gen_hamk(G, rdict, basis)
evals,evecs = Hamiltonian.SolveHk(kdict, rdict, basis)
Draw.plot_bandstructre(kdict,evals,"chinook")
Draw.plot_contour(kdict,rdict, [-1,1], [-1,1],"chinook", index=1)
opdict = Dict(
    "orbital"   =>  "px" 
)
Draw.plot_orbital_projection(kdict,basis,evals,evecs,opdict,"chinook")
end
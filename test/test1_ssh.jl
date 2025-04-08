module test0
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
lat = [[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
a = Structure.CreateAtom("a", [0,0,0], ["s","px"])
atoms = [a]
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0.1,"e_p" => 0.1),
    "ii" => Dict("V_sss"=>  -0.2,
                "V_sps"=>   -0.05,
                "V_pps"=>   0.2,
                "V_ppp"=>   0.0),
    "cutoff" => 1.2
)
spin = Dict(
    "soc"  =>  true,
    "λ"    => Dict("a"  =>  0.5)
)
rdict = Dict(
    "lattice"       =>  lattice,
    "atoms"         =>  atoms,
    "interactions"  =>  interactions,
    "spin"          =>  spin
)
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
basis = Hamiltonian.gen_basis(rdict) 
@show H = Hamiltonian.gen_hamk(Γ, rdict, basis)
evals,evecs = Hamiltonian.SolveHk(kdict,rdict, basis)
Draw.plot_bandstructre(kdict,evals,"ssh")
Draw.plot_contour(kdict,rdict, [-1,1], [-1,1],"ssh", index=1)
# orbital projection
opdict = Dict(
    "orbital"   =>  "s" 
)
Draw.plot_orbital_projection(kdict,basis,evals,evecs,opdict,"ssh")
end
module test0
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

# structure args
lat = [[1.0, 0.0, 0.0], [0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 50.0]]
atoms = []
for i in 1:20
    a = Structure.CreateAtom("C", [i/3.0, i/3.0, 0.01 * i], ["pz"])
    b = Structure.CreateAtom("C", [(i+1)/3.0, (i+1)/3.0, 0.01 * i], ["pz"])
    push!(atoms, a)
    push!(atoms, b)
end
lattice_constant = 1.0
lattice = Structure.CreateCell(lat, lattice_constant)
# tight-binding args
interactions = Dict(
    "i" => Dict("e_s" => 0.0,"e_p" => 0.0),
    "ii" => Dict("V_pps"=>   0.18,
                "V_ppp"=>   0.2),
    "cutoff" => 0.7
)
spin = Dict(
    "soc"  =>  false,
    "λ"    => Dict("C"  =>  0.0)
)
rdict = Dict(
    "lattice"       =>  lattice,
    "atoms"         =>  atoms,
    "interactions"  =>  interactions,
    "spin"          =>  spin
)
# kpath args
Γ = [0.0, 0.0, 0.0]
M = [0.5, 0.5, 0.0]
K = [2/3, 1/3, 0.0]
kptnum_perline = 60
b_matrix = Structure.ReciprocalLattice(lattice)
kdict = Dict(
    "pts"       =>  [Γ, K, M, Γ],
    "labels"    =>  ["Γ", "K", "M", "Γ"],
    "grid"      =>  kptnum_perline,
    "rvec"      =>  b_matrix
)

# diag and plot
basis = Hamiltonian.gen_basis(rdict) 
evals,evecs = Hamiltonian.SolveHk(kdict,rdict, basis)
Draw.plot_bandstructre(kdict,evals,"ABC")
# opdict = Dict(
#     "orbital"   =>  "s" 
# )
# Draw.plot_orbital_projection(kdict,basis,evals,evecs,opdict,"ssh")
end


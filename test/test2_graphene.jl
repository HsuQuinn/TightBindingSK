module test1
include("../src/klib.jl")
include("../src/param.jl")
include("../src/structure.jl")
include("../src/hamilton.jl")
include("../src/draw.jl")
using .Structure, .Hamiltonian, .Params, .Klib, .Draw

let
    # Step1: structure args
    lat = [[1.0, 0.0, 0.0], [0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 100.0]] 
    a = Structure.CreateAtom("C", [1.0/3.0,1.0/3.0,0], ["s","px","py","pz"])
    b = Structure.CreateAtom("C", [2.0/3.0,2.0/3.0,0], ["s","px","py","pz"])
    atoms = [a,b]
    lattice_constant = 1.0
    lattice = Structure.CreateCell(lat, lattice_constant)
    # Step2: tight-binding args
    interactions = Dict(
        "i" => Dict("e_s" => -8.81,"e_p" => -0.44),
        "ii" => Dict("V_sss"=>  -5.279,
                    "V_sps"=>   5.618,
                    "V_pps"=>   6.05,
                    "V_ppp"=>   -3.07),
        "cutoff" => 0.8 * lattice_constant 
    )
    spin = Dict(
        "soc"  =>  true,
        "λ"    => Dict("C"  =>  1.5)
    )
    rdict = Dict(
        "lattice"       =>  lattice,
        "atoms"         =>  atoms,
        "interactions"  =>  interactions,
        "spin"          =>  spin
    )
    # Step3: kpath args
    Γ = [0.0, 0.0, 0.0]
    M = [0.5, 0.5, 0.0]
    K = [2/3, 1/3, 0.0]
    kptnum_perline = 100
    b_matrix = Structure.ReciprocalLattice(lattice)
    kdict = Dict(
        "pts"       =>  [Γ, K, M, Γ],
        "labels"    =>  ["Γ", "K", "M", "Γ"],
        "grid"      =>  kptnum_perline,
        "rvec"      =>  b_matrix
    )
    # Step4: diag and plot
    basis = Hamiltonian.gen_basis(rdict) 
    evals,evecs = Hamiltonian.SolveHk(kdict,rdict, basis)
    Draw.plot_bandstructre(kdict,evals,"graphene_sppp")
    # Draw.plot_contour(kdict,rdict, [-1,1], [-1,1],"graphene_sppp", index=1)
    opdict = Dict(
        "orbital"   =>  "pz" 
    )
    Draw.plot_orbital_projection(kdict,basis,evals,evecs,opdict,"graphene_sppp")
end


# Only pz orbital
let
    # Step1: structure args
    lat = [[1.0, 0.0, 0.0], [0.5, sqrt(3)/2, 0.0], [0.0, 0.0, 100.0]] 
    a = Structure.CreateAtom("C", [1.0/3.0,1.0/3.0,0], ["pz"])
    b = Structure.CreateAtom("C", [2.0/3.0,2.0/3.0,0], ["pz"])
    atoms = [a,b]
    lattice_constant = 1.0
    lattice = Structure.CreateCell(lat, lattice_constant)
    # Step2: tight-binding args
    interactions = Dict(
        "i" => Dict("e_p" => 0),
        "ii" => Dict("V_ppp" => -0.2),
        "cutoff" => 0.8
    )
    spin = Dict(
        "soc"  =>  true,
        "λ"    => Dict("C"  =>  1.5)
    )
    rdict = Dict(
        "lattice"       =>  lattice,
        "atoms"         =>  atoms,
        "interactions"  =>  interactions,
        "spin"          =>  spin
    )
    # Step3: kpath args
    Γ = [0.0, 0.0, 0.0]
    M = [0.5, 0.5, 0.0]
    K = [2/3, 1/3, 0.0]
    kptnum_perline = 100
    b_matrix = Structure.ReciprocalLattice(lattice)
    kdict = Dict(
        "pts"       =>  [Γ, K, M, Γ],
        "labels"    =>  ["Γ", "K", "M", "Γ"],
        "grid"      =>  kptnum_perline,
        "rvec"      =>  b_matrix
    )
    # Step4: diag and plot
    basis = Hamiltonian.gen_basis(rdict) 
    evals,evecs = Hamiltonian.SolveHk(kdict,rdict, basis)
    Draw.plot_bandstructre(kdict,evals,"graphene_p")
    # Draw.plot_contour(kdict,rdict, [-1,1], [-1,1],"graphene_sppp", index=1)
    opdict = Dict(
        "orbital"   =>  "pz" 
    )
    Draw.plot_orbital_projection(kdict,basis,evals,evecs,opdict,"graphene_p")
end

end
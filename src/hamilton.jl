module Hamiltonian
include("param.jl")
include("klib.jl")
using LinearAlgebra
using .Params,.Klib

const ORBITALS_ALL = ["s", "px", "py", "pz", "dxy", "dyz", "dxz", "dx2-y2", "dz2", "S"]

struct Basis
    atom_index::Int64
    element::String
    orbital::String 
    spin::String
    pos::Vector{Float64}
end

function gen_basis(rdict)
    BasisArray = []
    lattice = rdict["lattice"]
    atoms = rdict["atoms"]
    spin = rdict["spin"]
    if spin["soc"] == true
        for (ia, atom) in enumerate(atoms)
            for orb in atom.orbitals
                push!(BasisArray,Basis(ia, atom.element, orb, "up", atom.pos))
                push!(BasisArray,Basis(ia, atom.element, orb, "down", atom.pos))
            end
        end
    else
        for (ia, atom) in enumerate(atoms)
            for orb in atom.orbitals
                push!(BasisArray,Basis(ia, atom.element, orb, "none", atom.pos))
            end
        end
    end
    return BasisArray
end

function expand_to_supercell_333(lattice, pos)
    expanded_positions = []
    for i in -1:1
        for j in -1:1
            for k in -1:1
                new_frac = pos .+ [i, j, k]
                new_cart = lattice.matrix * new_frac
                push!(expanded_positions, new_cart)
            end
        end
    end
    return expanded_positions
end

function map_orbital_to_index(orbital::String)::Int
    index = findfirst(x -> x == orbital, ORBITALS_ALL)
    if index === nothing
        error("Orbital '$orbital' is not in the list of valid orbitals: $ORBITALS_ALL")
    end
    return index
end

function gen_hamk(k, rdict, basis)
    lattice = rdict["lattice"]
    atoms = rdict["atoms"]
    interactions = rdict["interactions"]
    num_orbitals = sum(length(basis))
    H = zeros(ComplexF64, num_orbitals, num_orbitals)
    for (i1,o1) in enumerate(basis) 
        for (i2,o2) in enumerate(basis) 
            if o1.spin != o2.spin
                continue
            end
            pos1 = lattice.matrix * o1.pos
            expand_pos = expand_to_supercell_333(lattice, o2.pos)
            for pos2 in expand_pos
                distance = norm(pos1 - pos2)
                if distance > interactions["cutoff"]
                    continue
                end 
                if distance < 1E-6
                    if o1.orbital == o2.orbital 
                        H[i1, i2] += Params.onsite(o1.orbital,interactions)
                    end
                    continue
                end
                r = pos1 - pos2
                l = r[1]/norm(r)
                m = r[2]/norm(r)
                n = r[3]/norm(r)
                t_matrix = Params.get_hop_int(interactions,l,m,n)
                H[i1, i2] += t_matrix[map_orbital_to_index(o1.orbital), map_orbital_to_index(o2.orbital)] * exp(1.0im * dot(k, r))
            end
        end
    end
    if rdict["spin"]["soc"] == true 
        for (i1,o1) in enumerate(basis) 
            for (i2,o2) in enumerate(basis) 
                if o1.atom_index == o2.atom_index
                    if haskey(rdict["spin"]["λ"], o1.element) && first(o1.orbital) == 'p' && first(o2.orbital) == 'p'
                        H[i1, i2] += rdict["spin"]["λ"][o1.element] / 2 * Params.hsoc(o1.orbital, o2.orbital, o1.spin, o2.spin)
                    end
                end
            end
        end
    end
    return H
end

function SolveHk(kdict, rdict, basis)
    kpath, _ = Klib.gen_K(kdict)
    eigenvalues = []
    eigenstates = []
    for k in kpath
        Hk = gen_hamk(k, rdict, basis)
        evals, evecs = eigen(Hk)
        push!(eigenvalues, evals)
        push!(eigenstates, evecs)
    end
    return eigenvalues, eigenstates
end
end # module Hamiltonian
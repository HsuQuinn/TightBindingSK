module Hamiltonian
include("param.jl")
include("klib.jl")
using LinearAlgebra
using .Params,.Klib

const ORBITALS_ALL = ["s", "px", "py", "pz", "dxy", "dyz", "dxz", "dx2-y2", "dz2", "S"]

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


function gen_ham_k(k, lattice, atoms, interactions, cutoff::Float64)
    num_orbitals = sum(length(atom.orbitals) for atom in atoms)
    H = zeros(ComplexF64, num_orbitals, num_orbitals)
    orbital_start_indices = cumsum(vcat(0, [length(atom.orbitals) for atom in atoms]))
    for α in 1:length(atoms)
        for β in 1:length(atoms)
            orbital_α = atoms[α].orbitals
            orbital_β = atoms[β].orbitals
            subH = zeros(ComplexF64, length(orbital_α), length(orbital_β))
            pos_α = lattice.matrix * atoms[α].pos
            expand_pos_β = expand_to_supercell_333(lattice, atoms[β].pos)
            for pos_β in expand_pos_β
                distance = norm(pos_α - pos_β)
                if distance > cutoff
                    continue
                end
                if distance < 1E-6
                    for o in 1:length(orbital_α)
                        subH[o, o] += Params.onsite(orbital_α[o],interactions)
                    end
                    continue
                end
                r = pos_α - pos_β
                l = r[1]/norm(r)
                m = r[2]/norm(r)
                n = r[3]/norm(r)
                t_matrix = Params.get_hop_int(interactions,l,m,n)
                for oi in 1:length(orbital_α)
                    for oj in 1:length(orbital_β)
                        ii = map_orbital_to_index(orbital_α[oi])
                        jj = map_orbital_to_index(orbital_β[oj])
                        subH[oi, oj] += t_matrix[ii, jj] * exp(1.0im * dot(k, r))               
                    end
                end
            end
            start_row = orbital_start_indices[α] + 1
            end_row = orbital_start_indices[α + 1]
            start_col = orbital_start_indices[β] + 1
            end_col = orbital_start_indices[β + 1]
            H[start_row:end_row, start_col:end_col] += subH
        end
    end
    return H
end


function SolveHk(kdict, lattice, atoms, interactions, cutoff::Float64)
    kpath, _ = Klib.gen_K(kdict)
    eigenvalues = []
    for k in kpath
        Hk = gen_ham_k(k, lattice, atoms, interactions, cutoff::Float64)
        evals = eigen(Hk).values
        push!(eigenvalues, evals)
    end
    return eigenvalues
end
end # module Hamiltonian
module Structure
using LinearAlgebra

struct Lattice
    matrix::Matrix{Float64}   
    a::Float64               
    b::Float64               
    c::Float64               
    alpha::Float64            # b, c 
    beta::Float64             # c, a 
    gamma::Float64            # a, b 
end

const ORBITALS_ALL = ["s", "px", "py", "pz", "dxy", "dyz", "dxz", "dx2-y2", "dz2", "S"]

mutable struct  Atom
    element::String
    pos::Vector{Float64}
    orbitals::Union{Nothing, Vector{String}}
end

function CreateCell(matrix::Union{AbstractMatrix{T}, Vector{Vector{T}}}, lat_const::T) where T<:Real
    if isa(matrix, Vector{Vector{T}})
        matrix = hcat(matrix...)
    end
    scaled_matrix = matrix * lat_const
    a_vec = matrix[1, :]
    b_vec = matrix[2, :]
    c_vec = matrix[3, :]
    a = norm(a_vec)
    b = norm(b_vec)
    c = norm(c_vec)
    alpha = acos(dot(b_vec, c_vec) / (norm(b_vec) * norm(c_vec)))
    beta = acos(dot(c_vec, a_vec) / (norm(c_vec) * norm(a_vec)))
    gamma = acos(dot(a_vec, b_vec) / (norm(a_vec) * norm(b_vec)))
    return Lattice(Matrix{Float64}(scaled_matrix),
                   float(a), float(b), float(c),
                   float(alpha), float(beta), float(gamma))
end

function ReciprocalLattice(lattice::Lattice)::Matrix{Float64}
    a1 = lattice.matrix[:, 1]
    a2 = lattice.matrix[:, 2]
    a3 = lattice.matrix[:, 3]

    volume = dot(a1, cross(a2, a3))
    if iszero(volume)
        error("The lattice vectors are linearly dependent, cannot compute reciprocal lattice.")
    end

    b1 = 2π * cross(a2, a3) / volume
    b2 = 2π * cross(a3, a1) / volume
    b3 = 2π * cross(a1, a2) / volume

    return hcat(b1, b2, b3)
end

function CreateAtom(element::String, pos, orbitals::Union{Nothing, Vector{String}}=nothing)
    pos_array = Float64.(pos)
    if orbitals !== nothing
        if !(Set(orbitals) ⊆ Set(ORBITALS_ALL))
            error("wrong orbitals")
        end
    end
    return Atom(element, pos_array, orbitals)
end

end
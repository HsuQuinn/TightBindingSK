module Klib
include("structure.jl")
using LinearAlgebra
using .Structure

function gen_K(kdict::Dict)
    pts = kdict["pts"] 
    labels = kdict["labels"] 
    grid = kdict["grid"]
    rvec = kdict["rvec"]

    b1, b2, b3 = rvec[:,1], rvec[:,2], rvec[:,3]
    kpath = []
    kpath_rlt = []
    kdist = []
    kcoords = 0.0
    for i in 1:(length(pts)-1)
        k_start = pts[i]
        k_end = pts[i + 1]
        segment_relative = [k_start + t * (k_end - k_start) for t in range(0, 1, length=grid)]
        segment_absolute = [k[1] * b1 + k[2] * b2 + k[3] * b3 for k in segment_relative]
        append!(kpath, segment_absolute)
        append!(kpath_rlt, segment_relative)
    end
    for i in 1:(length(kpath_rlt)-1)
        deltal = norm(kpath_rlt[i+1]-kpath_rlt[i])
        append!(kdist, kcoords)
        kcoords += deltal
    end
    append!(kdist,kcoords)
    return kpath,kdist
end


end
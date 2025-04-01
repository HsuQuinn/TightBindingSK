module Draw
include("hamilton.jl")
include("klib.jl")
using LinearAlgebra, Plots
using .Klib, .Hamiltonian

function plot_bandstructre(kdict,evals,filename)
    kpath, kdist = Klib.gen_K(kdict)
    num_bands = length(evals[1])  
    num_kpoints = length(kpath)  
    label_string = join(kdict["labels"], "->")

    bands = [Float64[] for _ in 1:num_bands]
    for eval in evals
        for i in 1:num_bands
            push!(bands[i], real(eval[i]))
        end
    end

    plot(kdist, bands[1], label="Band 1", xlabel="k-path:$label_string", ylabel="Energy (eV)", lw=2)
    for i in 2:num_bands
        plot!(kdist, bands[i], label="Band $i", lw=2)
    end

    xticks = [0.0]
    for i in 2:length(kdict["pts"])
        push!(xticks, kdist[(i - 1) * kdict["grid"]])
    end
    vline!(xticks, label="", lw=1, lc=:black, ls=:dash)
    savefig("band_$filename.png")
end

function plot_contour(kdict::Dict, lattice,atoms,interactions,cutoff,xrange,yrange,filename; index=1)
    e = []
    rvec = kdict["rvec"]
    b1, b2, b3 = rvec[:,1], rvec[:,2], rvec[:,3]
    x_vals = range(xrange[1], xrange[2], length=40)
    y_vals = range(yrange[1], yrange[2], length=40)
    for i in x_vals
        for j in y_vals
            k_frac = [i, j, 0]
            k = k_frac[1] * b1 + k_frac[2] * b2 + k_frac[3] * b3 
            Hk = Hamiltonian.gen_ham_k(k, lattice, atoms, interactions, cutoff)
            evals = real(eigen(Hk).values)
            push!(e, evals[index])
        end
    end
    e_matrix = reshape(e, (40,40))
    fig2 = contourf(x_vals, y_vals, e_matrix, levels=10, color=:viridis)
    title!("Energy contour plot")
    savefig("contour_$filename.png")
end
end
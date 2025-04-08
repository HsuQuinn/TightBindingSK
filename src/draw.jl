module Draw
include("hamilton.jl")
include("klib.jl")
using LinearAlgebra, Plots, Colors, ColorSchemes, PyPlot
using .Klib, .Hamiltonian

cdict = Dict(
    "red" => [
        (0.0,   1.0, 1.0),   # 0~0.5：白 → 深蓝（红通道从1降到0.0）
        (0.5,   0.0, 0.0),   
        (1.0,   0.0, 0.0)    # 0.5~1.0：深蓝（红通道保持0）
    ],
    "green" => [
        (0.0,   1.0, 1.0),   # 0~0.5：白 → 深蓝（绿通道从1降到0.0）
        (0.5,   0.0, 0.0),   
        (1.0,   0.0, 0.0)    # 0.5~1.0：深蓝（绿通道保持0）
    ],
    "blue" => [
        (0.0,   1.0, 1.0),   # 0~0.5：白 → 深蓝（蓝通道从1保持1）
        (0.5,   1.0, 1.0),   
        (1.0,   0.5, 0.5)    # 0.5~1.0：深蓝（蓝通道从1降到0.5）
    ],
    "alpha" => [
        (0.0,   0.0, 0.0),   # 透明度从0（完全透明）
        (0.5,   0.5, 0.5),   # 渐变到0.5（半透明）
        (1.0,   1.0, 1.0)    # 到1（完全不透明）
    ]
)
cmap_bwr = matplotlib.colors.LinearSegmentedColormap("TransparentWhiteToDeepBlue", cdict)

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
    for i in 1:num_bands
        PyPlot.plot(kdist, bands[i], label="Band $i", linewidth=1.5)
    end
    xticks = [0.0]
    for i in 2:length(kdict["pts"])
        push!(xticks, kdist[(i - 1) * kdict["grid"]])
    end
    for x in xticks
        PyPlot.axvline(x=x, color="black", linestyle="--", linewidth=0.8)
    end
    PyPlot.xlabel("k-path: $label_string")
    PyPlot.ylabel("Energy (eV)")
    PyPlot.title("Band Structure")
    PyPlot.legend()
    PyPlot.tight_layout()
    PyPlot.savefig("band_$filename.png")
    PyPlot.close()
end

function plot_orbital_projection(kdict, basis,evals,evecs,opdict,filename)
    orbs = []
    op = opdict["orbital"]
    for (i,b) in enumerate(basis)
        if b.orbital == op 
            push!(orbs,i)
        end
    end
    label_string = join(kdict["labels"], "->")
    kpath, kdist = Klib.gen_K(kdict)
    num_bands = length(evals[1])  
    num_kpoints = length(kpath)  
    bands = [Float64[] for _ in 1:num_bands]
    weights = [Float64[] for _ in 1:num_bands]
    for k_index in 1:num_kpoints
        eval = evals[k_index]
        evec = evecs[k_index]
        for i in 1:num_bands
            push!(bands[i], real(eval[i]))
            push!(weights[i], sum(abs2.(vec(evec[orbs, i]))))          
        end
    end
    for i in 1:num_bands
        PyPlot.scatter(kdist, bands[i], c=weights[i], cmap=cmap_bwr, s=10, edgecolors="none", norm=matplotlib.colors.Normalize(vmin=0, vmax=1))
    end
    for i in 1:num_bands
        PyPlot.plot(kdist, bands[i], label="Band $i", color="grey", linewidth=1.0, alpha=0.8)  
    end
    colorbar(label="Orbital Weight")
    xlabel("k-path: $label_string")
    ylabel("Energy (eV)")
    title("Orbital Projection: $op")
    PyPlot.savefig("$op orbital projection_$filename.png")
    close()
end


function plot_contour(kdict::Dict, rdict , xrange,yrange,filename; index=1)
    e = []
    basis = Hamiltonian.gen_basis(rdict) 
    rvec = kdict["rvec"]
    b1, b2, b3 = rvec[:,1], rvec[:,2], rvec[:,3]
    x_vals = range(xrange[1], xrange[2], length=40)
    y_vals = range(yrange[1], yrange[2], length=40)
    for i in x_vals
        for j in y_vals
            k_frac = [i, j, 0]
            k = k_frac[1] * b1 + k_frac[2] * b2 + k_frac[3] * b3 
            Hk = Hamiltonian.gen_hamk(k, rdict, basis)
            evals = real(eigen(Hk).values)
            push!(e, evals[index])
        end
    end
    e_matrix = reshape(e, (40,40))
    fig2 = Plots.contourf(x_vals, y_vals, e_matrix, levels=10, color=:viridis)
    title!("Energy contour plot")
    Plots.savefig("contour_$filename.png")
end
end
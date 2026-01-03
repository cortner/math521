
# using QuadGK, LinearAlgebra, FFTW 

# -------------------------------------------- 

"""
    function eval_p1(x, U)

evaluates a piecewise affine function on the interval (-pi, pi) 
with uniform grid spacing and nodal values U at the point x.
"""
function eval_p1(x, U)
    N = div( (length(U) - 1), 2 )
    h = pi / N 
    i1 = floor(Int, (x+pi) * 2*N / (2*pi) ) + 1
    if i1 < 1; i1 = 1; end
    i2 = i1+1 
    if i2 > length(U)
        return U[end] 
    end 
    x1 = -pi + h * i1
    return U[i1] + (x - x1) * (U[i2] - U[i1]) / h
end

# -------------------------------------------- 
#
#   Makie & Ferrite Plotting Tools 
#

using GeometryBasics: TriangleFace

# this is a little wrapper around Makie.mesh to plot 
# P1 functions on unstructured meshes
function trisurf!(ax, grid, U; kwargs...)
    x = [n.x[1] for n in grid.nodes]
    y = [n.x[2] for n in grid.nodes]
    z = U
    faces = TriangleFace.([ c.nodes for c in grid.cells ])
    mesh!(ax, x, y, z, faces, color=z, kwargs...)
end

function trisurf(grid, U; size = (300, 300), kwargs...)
    fig = Figure(size = size)
    ax = Axis3(fig[1,1]; )
    trisurf!(ax, grid, U; kwargs...) 
    return fig
end

function trimesh!(ax, grid; kwargs...)
    x = Float64[]
    y = Float64[] 
    for c in grid.cells 
        for (i, j) in zip(c.nodes, circshift(c.nodes, 1))
            push!(x, grid.nodes[i].x[1])
            push!(x, grid.nodes[j].x[1])
            push!(y, grid.nodes[i].x[2])
            push!(y, grid.nodes[j].x[2])
        end
    end 
    return linesegments!(ax, x, y; kwargs...)
end

function trimesh(grid; size = (300, 300), kwargs...)
    fig = Figure(size = size)
    ax = Axis(fig[1,1]; )
    hidedecorations!(ax); hidespines!(ax)
    trimesh!(ax, grid; kwargs...)
    return fig 
end



"""
simple utility function to plot a Ferrite.jl P1 FEM solution
"""
function trisurf_ferrite(u, dh, pgrid = dh.grid; 
                    sym = :u, size = (400, 400), kwargs...) 
    fig = Figure(size = size)
    ax = Axis3(fig[1,1]; )
    
    # evaluate the function at the nodes of the plotting grid 
    # sounds stupid when it's the same grid, but it is just 
    # a simple trick to get the right ordering which ferrite 
    # messes up internally ... 
    pnodes = [ Tensors.Vec(n.x...) for n in pgrid.nodes ]
    ph = PointEvalHandler(dh.grid, pnodes)
    pu = evaluate_at_points(ph, dh, u, :u)
    
    trisurf!(ax, pgrid, pu; kwargs...) 
    return fig 
end


# -------------------------------------------- 
#
#   Trigonometric interpolation tools 
#


#=

"""
This is a very naive method to compute the 
Fourier coefficients. We will soon learn about 
much more efficient schemes based on the FFT.
"""
function compute_fcoeffs_quad(f_fun, N)
    fhat = zeros(ComplexF64, 2*N+1)
    kgrid = -N:N
    for (i, k) in enumerate(kgrid)
        g = x -> f_fun(x) * exp(-im * k * x)
        fhat[i] = quadgk(g, -pi, pi; rtol = 1e-6, atol=1e-6)[1] / (2*pi)
    end
    return fhat, kgrid
end


"interpolation nodes"
xgrid(N) = [ j * π / N  for j = 0:2N-1 ]

"fourier coefficient indices"
kgrid(N) = [ 0:N; -N+1:-1 ]


"""
construct the coefficients of the trigonometric interpolant
"""
function triginterp_fft(f::Function, N)
    X = xgrid(N)
    # nodal values at interpolation nodes
    F = f.(X) 
    return fft(F) / (2*N)    
end 


"""
to evaluate a trigonometric polynomial just sum coefficients * basis
we the take the real part because we assume the function we are 
approximating is real.
"""
evaltrig(x, F̂) = sum( real(F̂k * exp(im * x * k))
                      for (F̂k, k) in zip(F̂, kgrid(length(F̂) ÷ 2)) )


function evaltrig_grid(F̂::AbstractVector, M::Integer)
    N = length(F̂) ÷ 2
    @assert 2 * N == length(F̂)
    @assert M >= N 
    F̂_M = zeros(ComplexF64, 2*M)
    F̂_M[1:N] .= F̂[1:N]
    F̂_M[end-N+1:end] .= F̂[end-N+1:end]
    x = xgrid(M) 
    Fx = real.(ifft(F̂_M) * 2*M)
    return Fx, x
end

# ------------------------ 
# Utilities for dim > 1 

"""
Given a one-dimensional array y, return d d-dimensional arrays 
 y ⊗ 1 ⊗ ... ⊗ 1   (x1-coordinate)
 1 ⊗ y ⊗ 1 ⊗ ...   (x2-coordinate)
... 
 1 ⊗ ... ⊗ 1 ⊗ y   (xd-coordinate)
"""
function tensorgrid(d, x1)
    dims = ntuple(i -> length(x1), d)
    X = reshape(x1 * ones(Bool, length(x1)^(d-1))', dims)
    pdim(i, d) = (dd = collect(1:d); dd[1] = i; dd[i] = 1; tuple(dd...))
    return ntuple(i -> permutedims(X, pdim(i,d)), d)
end

"""
d-dimensional x grid 
"""
xgrid(d, N) = tensorgrid(d, xgrid(N))

"""
d-dimensional k-grid 
"""
kgrid(d, N) = tensorgrid(d, kgrid(N))


"""
construct the coefficients of the trigonometric interpolant
in d dimensions
"""
function triginterp_fft(f::Function, N, d::Integer)
    XX = xgrid(d, N)
    # nodal values at interpolation nodes
    F = f.(XX...) 
    return fft(F) / (2*N)^d
end 

function evaltrig_grid(F̂::AbstractArray{T, 2}, M::Integer) where {T}
    N = size(F̂, 1) ÷ 2;
    @assert size(F̂) == (2*N, 2*N)
    @assert M >= N
    F̂_M = zeros(ComplexF64, (2*M, 2*M)) 
    kk1 = 1:N; kk2 = N+1:2*N; kk3 = 2*M-N+1:2*M
    F̂_M[kk1, kk1] .= F̂[kk1, kk1]
    F̂_M[kk1, kk3] .= F̂[kk1, kk2]
    F̂_M[kk3, kk1] .= F̂[kk2, kk1] 
    F̂_M[kk3, kk3] .= F̂[kk2, kk2]
    x = xgrid(M) 
    Fx = real.(ifft(F̂_M) * (2*M)^2)
    return Fx, x
end

function evaltrig_grid(F̂::AbstractArray{T, 3}, M::Integer) where {T}
    N = size(F̂, 1) ÷ 2;
    @assert size(F̂) == (2*N, 2*N, 2*N)
    @assert M >= N
    F̂_M = zeros(ComplexF64, (2*M, 2*M, 2*M))
    kk1 = 1:N; kk2 = N+1:2*N; kk3 = 2*M-N+1:2*M
    F̂_M[kk1, kk1, kk1] .= F̂[kk1, kk1, kk1]
    F̂_M[kk1, kk1, kk3] .= F̂[kk1, kk1, kk2]
    F̂_M[kk1, kk3, kk1] .= F̂[kk1, kk2, kk1]
    F̂_M[kk1, kk3, kk3] .= F̂[kk1, kk2, kk2]
    F̂_M[kk3, kk1, kk1] .= F̂[kk2, kk1, kk1]
    F̂_M[kk3, kk1, kk3] .= F̂[kk2, kk1, kk2]
    F̂_M[kk3, kk3, kk1] .= F̂[kk2, kk2, kk1]
    F̂_M[kk3, kk3, kk3] .= F̂[kk2, kk2, kk2]
    x = xgrid(M) 
    Fx = real.(ifft(F̂_M) * (2*M)^3)
    return Fx, x
end


=#

nothing

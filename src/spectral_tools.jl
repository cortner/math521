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


"""
Evaluate the trigonometric polynomial with coefficients `F̂` on 
an equispaced grid with `2M` gridpoints. 
"""
function evaltrig_grid(F̂, M::Integer; convert = real)
    N = length(F̂) ÷ 2
    @assert 2 * N == length(F̂)
    @assert M >= N 
    F̂_M = zeros(ComplexF64, 2*M)
    F̂_M[1:N] .= F̂[1:N]
    F̂_M[end-N+1:end] .= F̂[end-N+1:end]
    x = xgrid(M) 
    Fx = real.(ifft(F̂_M) * (2*M))
    return Fx, x
end

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
returns X1, X2, ...
"""
xgrid(d, N) = tensorgrid(d, xgrid(N))

"""
d-dimensional k-grid 
returns K1, K2, ...
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

"""
Evaluate a multivariate trigonometric polynomial on a 
refined grid. 
"""
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

;
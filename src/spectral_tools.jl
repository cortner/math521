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
;
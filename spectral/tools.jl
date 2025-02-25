
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
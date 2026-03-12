
function  aitken(x::AbstractVector)
    n = length(x)
    d = x[2:n] - x[1:n-1]     # length n-1
    dd = d[2:n-1] - d[1:n-2]  # length n-2
    y = x[2:n-1] - d[1:n-2] .* d[2:n-1] ./ dd
    return y[end]
end


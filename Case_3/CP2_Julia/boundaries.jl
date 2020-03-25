function boundaries(X::AbstractArray)
    row, col = size(X)
    X_new  = zeros(row+2, col +2)
    for i in 1:row
        for j in 1:col
            X_new[i+1,j+1] = X[i,j]
        end
    end
    return X_new
end

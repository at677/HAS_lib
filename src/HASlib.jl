module HASlib
function hexgrid(n::Int64)
    a = [1.0 -0.5;
         0.0 sqrt(3.0)/2.0]

    max_index = 2*n+1

    nelements = 1+6*n+3*(n)*(n-1)
    g = zeros(Float64,nelements,2)
    gi = zeros(Int32,nelements,2)
    index = 1
    for i = 1:max_index
        maxj = i+n    
        if maxj > max_index
            maxj = max_index
        end
        minj = i-n    
        if minj < 1
            minj = 1
        end
        for j = minj:maxj
            t = a*[i-n-1 ; j-n-1]
            g[index,:]=t
            gi[index,:]=[i-n-1; j-n-1]
            index+=1
        end
    end
    return g,gi
end
end

module HASlib
export hexgrid, h2m, cc
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

function h2m(m)
    hb = 1.054571800e-34
    hb^2/(2*m)
end

function cc(Ei,theta,m,maxgrid=5,st_p_wl=100)

    Ei = Ei*1.6e-19

    ki = sqrt(Ei/h2m(m))
    Ki = ki*sin(theta*180/π)

    g,gi = hexgrid(maxgrid)
    g = g.*1e10

    kz2 = -((g[:,1].+Ki).^2 .+ g[:,2].^2).+ ki^2

    index = sortperm(kz2)
    sorted_kz2 = kz2[index]
    
    #find lowest open channel
    min_index_open = argmax(sorted_kz2 .>= 0.0)
    println(sorted_kz2[min_index_open])
    println(min_index_open)
    println(sorted_kz2[min_index_open-1])

    e_max = maximum(kz2)
    step = 2*π / (st_p_wl*sqrt(e_max))

    #coupling



end

end

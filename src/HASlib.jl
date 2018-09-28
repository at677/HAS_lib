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

function W(z)
    I*center_pot(z,D,kappa)/h2m(m) - M_K
end

function cc(Ei,theta,m,maxgrid=5,st_p_wl=100,n_open=10,n_close=10)

    Ei = Ei*1.6e-19

    ki = sqrt(Ei/h2m(m))
    Ki = ki*sin(theta*180/π)

    g,gi = hexgrid(maxgrid)
    g = g.*1e10

    
    kz2 = -((g[:,1].+Ki).^2 .+ g[:,2].^2).+ ki^2

    index = sortperm(kz2)
    sorted_kz2 = kz2[index]
    
    #find lowest open channel
    min_index_open = findfirst(x -> x >= 0, sorted_kz2)
    println(sorted_kz2[min_index_open],"->",min_index_open,"\n",
            sorted_kz2[min_index_open-1],"->",min_index_open-1)


    open_index = min_index_open:length(sorted_kz2)
    if n_open >= 0 & length(sorted_kz2) <= min_index_open+n_open
        open_index = min_index_open:min_index_open+n_open
    end

    close_index = 1:min_index_open-1
    if n_close >= 0 &  min_index_open-1-n_close <= 1
        close_index = min_index_open-1-n_close:min_index_open-1
    end

    println("open: ",open_index," close:",close_index)
    
    open_kz2 = sorted_kz2[open_index]
    open_g = g[index,:][open_index,:]
    open_gi = gi[index,:][open_index,:]

    close_kz2 = sorted_kz2[close_index]
    close_g = g[index,:][close_index,:]
    close_gi = gi[index,:][close_index,:]

    display(hcat(open_kz2,open_gi,open_g))

    # calculate coupling

    n_total = length(open_kz2)+length(close_kz2)

    M_0 = I
    M_K = Array{Float64,2}(undef,n_total,n_total)
    M_E = Array{Float64,2}(undef,n_total,n_total)


    e_max = maximum(kz2)
    step = 2*π / (st_p_wl*sqrt(e_max))

    #coupling



end

end

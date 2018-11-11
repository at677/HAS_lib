module HASlib
export hexgrid, h2m, cc, fox_goodwin, CCSim

using LinearAlgebra

struct CCSim
    n_total
    n_open
    n_closed
    kz2
    g
    gi
end

function Base.show(result::CCSim)
    print("total channels: ",result.n_total)
    print(" (open: ",result.n_open)
    print(", closed: ",result.n_closed,")")
    print(", kz2: ",result.kz2,")")
end

"""
Calculate a hexagonal grid with a radius n and a lattice constant of 1
angle between a1 and a2 = 120 deg

# Arguments
- `n:Int64`
- `a₁=1`: base vector length 1
- `a₂=1`: base vector length 2

pattern (n = 3):
```
  + + +
 + + + +
+ + + + +
 + + + +
  + + +
```

# returns (g,gi)
 - g .. 2d array of lattice point coordinates
 - gi ... index of lattice point
"""
function hexgrid(n::Int64, a₁ = 1,a₂ = 1)
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


function fox_goodwin_step!(w_p1,w_0,w_m1,r,z,i)
        aa = 2*I+10*w_0
        bb = dot(I-w_m1,r)
        cc = I-w_p1
        r = (aa.-bb)/cc
        w_m1 = w_0
        w_0 = w_p1
        r
end

function fox_goodwin(z,w)
    h = z[2]-z[1]
    c1 = h^2/12

    w_m1 = c1*w(z[1])
    w_0 = c1*w(z[2])

    nr = size(w_0)[1]

    r  = fill!(similar(w_0),0)
    w_p1 = similar(w_0)

    for i in 2:length(z)-1
        w_p1 = c1*w(z[i+1])
        fox_goodwin_step!(w_p1,w_0,w_m1,r,z,i)
    end
    r
end


function cc(Ei,theta,m,z,maxgrid=5,st_p_wl=100,n_open=10,n_close=10)

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
        open_index = min_index_open:min_index_open+n_open-1
    end

    close_index = 1:min_index_open-1
    if n_close >= 0 &  min_index_open-1-n_close <= 1
        close_index = min_index_open-n_close:min_index_open-1
    end

    println("open: ",open_index," close:",close_index)
    
    open_kz2 = sorted_kz2[open_index]
    open_g = g[index,:][open_index,:]
    open_gi = gi[index,:][open_index,:]

    close_kz2 = sorted_kz2[close_index]
    close_g = g[index,:][close_index,:]
    close_gi = gi[index,:][close_index,:]

    # calculate coupling

    n_total = length(open_kz2)+length(close_kz2)
    println(n_total)

    M_0 = I
    M_K = Array{Float64,2}(undef,n_total,n_total)
    M_E = Array{Float64,2}(undef,n_total,n_total)

    e_max = maximum(kz2)
    step = 2*π / (st_p_wl*sqrt(e_max))

    #coupling
    w(z) = M_0.*z + M_K

    r = fox_goodwin(z,w)


    kz2 = sorted_kz2[min_index_open-n_close:min_index_open+n_open-1]

    S0_m1 = zeros(size(r)[1])
    S0_00 = zeros(size(r)[1])
    CE_m1 = zeros(size(r)[1])
    CE_00 = zeros(size(r)[1])

    for (i,erg) in enumerate(kz2)
        if erg > 0
            kz = sqrt(erg)
            S0_m1[i] = sin(kz*z[end-1])/sqrt(kz)
            S0_00[i] = sin(kz*z[end])/sqrt(kz)
            CE_m1[i] = cos(kz*z[end-1])/sqrt(kz)
            CE_00[i] = cos(kz*z[end])/sqrt(kz)
        else
            kz = sqrt(abs(erg))
            CE_m1[i] = exp(-kz*z[end-1])/sqrt(kz)
            CE_00[i] = exp(-kz*z[end])/sqrt(kz)
        end
    end

    S0_m1 = Diagonal(S0_m1)
    S0_00 = Diagonal(S0_00)
    CE_m1 = Diagonal(CE_m1)
    CE_00 = Diagonal(CE_00)

    SS = S0_m1 .- dot(r,S0_00)
    CC = dot(r,CE_00) .- CE_m1

    KK = CC\SS

    K_open = KK[n_close+1:end,n_close+1:end]
    S = broadcast(abs2,(I-K_open*(0+1im))\(I+K_open*(0+1im)))

    spec = (open_gi[:,1] .== 0) .& (open_gi[:,2] .== 0)

    hcat(S[spec,:]',open_gi,open_kz2)

    result = CCSim(n_total,
                   n_open,
                   n_close,
                   vcat(close_kz2, open_kz2),
                   vcat(close_g, open_g),
                   vcat(close_gi, open_gi))
end

end

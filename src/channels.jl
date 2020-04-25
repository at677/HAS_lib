using Parameters
using Printf
using LinearAlgebra

export calc_radius_factor, sort_and_cut, get_channels

"""
Object containing the scattering channels and the asympotic energy
"""
struct ScatteringChannels
    "reciprical lattice points"
    g
    "reciprical lattice index"
    gi
    "kz2 after scattering"
    kz2::Array{Real}
end

function Base.show(ch::ScatteringChannels)
    @unpack g,gi,kz2 = ch
    println("g","\t","index","\t","kz²")
    for i in 1:length(g)
        @printf("%8.2f %8.2f %3d %3d %8.2f\n", g[i][1],g[i][2],gi[i][1],gi[i][1],kz2[i])
        end
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



"""
# Arguments
* lattice Lattice object which is used for the reciprical lattice vectors
* ki wavevector of incident beam
* theta angle of incident beam and the surface normal in radians
* n range of generated channels
"""
function get_channels(b1,b2,ki, theta, n::T = 10) where (T <: Int)
    g = []
    gi = []
    kz2= []
    for i in -n:n, j in -n:n
        push!(gi,[i,j])
        cur_g = b1.*i+b2.*j
        push!(g,cur_g)
        tmp1 = ki.*[sin(theta),0].+cur_g
        push!(kz2,ki^2-tmp1'tmp1)
    end
    ScatteringChannels(g,gi,kz2)
end

function filter_channels_radius(ch::ScatteringChannels,
                                max_norm)
    @unpack kz2,g,gi = ch
    new_g = []
    new_gi = []
    new_kz²= []
    for (i,cur_g) in enumerate(g)
        if norm(cur_g) <= max_norm
            push!(new_g,cur_g)
            push!(new_gi,gi[i])
            push!(new_kz²,kz2[i])
        end
    end
    ScatteringChannels(new_g,new_gi,new_kz²)
end

"""
Calculate radius of a circle fitting in 4 cells
cells are drawn by g = i*b1+j+b2 for i and j in -1:1
"""
function calc_radius_factor(lattice::Generic2DLattice)
    @unpack b1,b2 = lattice
    d1 = [b1[2],-b1[1]]./norm(b1)
    d2 = [b2[2],-b2[1]]./norm(b2)
    r1 = abs(d1'b2)
    r2 = abs(d2'b1)
    min(r1,r2)
end

"""
sort by kz² and cut away channels
# Arguments
* ch channels
* n_open number of open channels to keep
* n_closed number of closed channels to keep

"""
function sort_and_cut(ch::ScatteringChannels, n_open::T, n_closed::U) where {T <: Int, U <:Int}

    @unpack g,gi,kz2 = ch

    index = sortperm(collect(zip(kz2,getindex.(gi,1),getindex.(gi,2))),rev=true)
    sorted_kz2 = kz2[index]
    
    #find lowest open channel
    first_closed = findfirst(x -> x < 0, sorted_kz2)

    open_index = 1:first_closed-1
    if n_open >= 0 && n_open < length(open_index)
        open_index = (first_closed-n_open:first_closed-1)
    end

    close_index = first_closed:length(kz2)
    if n_closed >= 0 &&  n_closed < length(close_index)
        close_index = first_closed:first_closed+n_closed-1
    end

    selected_index = vcat(open_index,close_index)

    ScatteringChannels(g[index,:][selected_index],
                       gi[index,:][selected_index],
                       sorted_kz2[selected_index])

end


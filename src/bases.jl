using Parameters
using Printf
using LinearAlgebra

export Generic2DLattice, calc_radius_factor, sort_and_cut, get_channels

"""
Object containing the scattering channels and the asympotic energy
"""
struct ScatteringChannels
    "reciprical lattice points"
    g
    "reciprical lattice index"
    gi
    "kz² after scattering"
    kz²
end

function Base.show(ch::ScatteringChannels)
    @unpack g,gi,kz² = ch
    println("g","\t","index","\t","kz²")
    for i in 1:length(g)
        @printf("%8.2f %8.2f %3d %3d %8.2f\n", g[i][1],g[i][2],gi[i][1],gi[i][1],kz²[i])
        end
end

abstract type Abstract2DLattice end 

"""
Structure containing the real and the reciprical basis vectors of a surface
"""
struct Generic2DLattice <: Abstract2DLattice
    a1
    a2
    b1
    b2
end

"""
# Arguments
* a1 first lattice vector
* a2 second lattice vector

# Example
` Generic2DLattice([1,0],[1,0]) `
"""
function Generic2DLattice(a1, a2)
    area = a2[2]*a1[1]-a2[1]*a1[2]
    b1 = 2π * [ a2[2],-a2[1]] / area
    b2 = 2π * [-a1[2], a1[1]] / area
    Generic2DLattice(a1,a2,b1,b2)
end

"""
# Arguments
* lattice Lattice object which is used for the reciprical lattice vectors
* ki wavevector of incident beam
* theta angle of incident beam and the surface normal in radians
* n range of generated channels
"""
function get_channels(lattice::Generic2DLattice,
                      ki, theta, n::T = 10) where (T <: Int)
    @unpack b1,b2 = lattice
    g = []
    gi = []
    kz²= []
    for i in -n:n, j in -n:n
        push!(g,[i,j])
        cur_g = b1*i+b2*j
        push!(gi,cur_g)
        tmp1 = ki.*[sin(theta),cos(theta)].+cur_g
        push!(kz²,ki^2-tmp1'tmp1)
    end
    ScatteringChannels(g,gi,kz²)
end

function filter_channels_radius(ch::ScatteringChannels,
                                max_norm)
    @unpack kz²,g,gi = ch
    new_g = []
    new_gi = []
    new_kz²= []
    for (i,cur_g) in enumerate(g)
        if norm(cur_g) <= max_norm
            push!(new_g,cur_g)
            push!(new_gi,gi[i])
            push!(new_kz²,kz²[i])
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

    @unpack g,gi,kz² = ch

    index = sortperm(kz²)
    sorted_kz2 = kz²[index]
    
    #find lowest open channel
    min_index_open = findfirst(x -> x >= 0, sorted_kz2)

    open_index = min_index_open:length(sorted_kz2)
    if n_open >= 0 && length(sorted_kz2) >= min_index_open+n_open - 1
        open_index = min_index_open:min_index_open+n_open-1
    end

    close_index = 1:min_index_open-1
    if n_closed >= 0 &&  min_index_open-1-n_closed >= 0
        close_index = min_index_open-n_closed:min_index_open-1
    end

    selected_index = vcat(close_index,open_index)

    ScatteringChannels(g[index,:][selected_index],
                       gi[index,:][selected_index],
                       sorted_kz2[selected_index])

end


using LinearAlgebra

export Generic2DLattice

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

struct HexGrid <: Abstract2DLattice
    a
    l::Generic2DLattice
end

HexGrid(a) = HexGrid(a,Generic2DLattice([a,0],[cos(pi/180*120),sin(pi/180*120)].*a))

function rotate_reciprocal(lattice::Generic2DLattice,direction)
	unitx = [lattice.b1 lattice.b2] * direction
	ang = -atan(unitx[2],unitx[1])
	R = [cos(ang) -sin(ang)
	     sin(ang) cos(ang)]
	R*lattice.b1,R*lattice.b2
end
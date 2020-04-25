using Test
using PyCall

haslib = pyimport("haslib")

mel = 1.6021766208000002e-22
A = 1e-10
a = 4.3*A
Ei = 4.196*mel
theta = 45.75

D = 5*mel
h = 0.1
xi = 0.38/A

crystal=haslib.Crystal(a, a, 120)
crystal.calc_reciprocals()
grid = HASlib.HexGrid(a)

@test grid.l.a1 ≈ crystal.a1[1:2]
@test grid.l.a2 ≈ crystal.a2[1:2]
@test grid.l.b1 ≈ crystal.b1_orig[1:2]
@test grid.l.b2 ≈ crystal.b2_orig[1:2]

function test_rotation(dir)
	crystal=haslib.Crystal(a, a, 120)
	grid = HASlib.HexGrid(a)
	crystal.rotate_crystal(dir...)
	b1,b2 = HASlib.rotate_reciprocal(grid.l,dir)
	(b1 ≈ crystal.b1[1:2]) && (b2 ≈ crystal.b2[1:2])
end


@test test_rotation([2,0])
@test test_rotation([0,2])
@test test_rotation([2,1])
@test test_rotation([2,-1])

b1,b2 = HASlib.rotate_reciprocal(grid.l,[2,-1])
mb1,mb2 = HASlib.rotate_reciprocal(grid.l,[-2,1])
@test b1 ≈ -mb1
@test b2 ≈ -mb2


b1,b2 = HASlib.rotate_reciprocal(grid.l,[1,0])
ERG = haslib.icc_2d_2.iCC_2D_2(h,1.0,a,D,xi,Ei,theta,b1,b2,floq=0,Gope=-1,Gclo=50,maxrad=20,sta_z=-1*A,end_z=6*A,st_p_wl=150)

m_He3 = 5.0082343984338715e-27
ki = sqrt(Ei/HASlib.h2m(m_He3))

channels = HASlib.get_channels(b1,b2,ki,theta*pi/180,5)
channels = HASlib.sort_and_cut(channels,-1,50)
result = cc(channels,a,h,D,xi,theta,m_He3,(-1*A,6*A),st_p_wl=150)

@test sort(result[:,1]) ≈ sort(ERG[:,7])

using BenchmarkTools
display(@benchmark cc(channels,a,h,D,xi,theta,m_He3,(-1*A,6*A),st_p_wl=150))
display(@benchmark haslib.CC_2D_old(h,a,D,xi,Ei,theta,b1,b2,G_ope=-1,G_clo=50,maxrad=20,sta_z=-1*A,end_z=6*A,st_p_wl=150))
display(@benchmark haslib.icc_2d_2.iCC_2D_2(h,1,a,D,xi,Ei,theta,b1,b2,Gope=-1,Gclo=50,floq=0,maxrad=20,sta_z=-1*A,end_z=6*A,st_p_wl=150))

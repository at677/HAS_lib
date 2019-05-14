
mel = 1.6021766208000002e-22
A = 1e-10

a = 4.3*A
Ei = 4.196*mel
theta = 45.75
D = 5*mel
h = 0.1
xi = 0.38/A

m_He3 = 5.0082343984338715e-27

ki = sqrt(Ei/HASlib.h2m(m_He3))
grid = HASlib.HexGrid(a)
b1,b2 = HASlib.rotate_reciprocal(grid.l,[1,0])
channels = HASlib.get_channels(b1,b2,ki,theta*pi/180,5)
channels = HASlib.sort_and_cut(channels,5,0)

@test_nowarn cc(channels,a,h,D,xi,theta,m_He3,(-1*A,6*A))
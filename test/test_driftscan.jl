
param = DriftScan(
m = m_He3,
theta = 45.75,
D = 5*mel,
h = 0.1,
xi = 0.38/A,
a = 4.3*A,
energies = range(4*mel,5*mel,length=10),
z = (-1*A,6*A)
)

@test_nowarn driftscan(param)

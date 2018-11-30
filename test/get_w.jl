using HASlib

a = 1
k = 1

hg = HASlib.HexGrid(a)
chan = HASlib.get_channels(hg.l, k, 30, 1)
corr = HASlib.CorrugatedMorse(1,1/7,0.001,hg)
w = HASlib.get_w(0,chan,k,corr)



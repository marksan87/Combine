#!/usr/bin/env python
from ROOT import *

_file = TFile("higgsCombineTest.MultiDimFit.mH120.root","read")


_tree = _file.Get("limit")

print _tree.GetEntries()



_MT = []
_deltaNLL  = []

for event in _tree:
    #if event.MT > 1725 and event.MT < 1745:
    #    _MT.append(event.MT / 10.0)
    #    _deltaNLL.append(2*event.deltaNLL)
    _MT.append(event.MT / 10.0)
    _deltaNLL.append(2*event.deltaNLL)



from array import array

MT = array('d',_MT)
deltaNLL = array('d',_deltaNLL)

obs = "p_{T}(ll)"
graph = TGraph(len(MT),MT,deltaNLL)
graph.SetTitle("%s Likelihood Fit" % obs)
graph.GetXaxis().SetTitle("m_{t} [GeV]")
graph.GetYaxis().SetTitle("-2 * NLL")
c = TCanvas("c","c",1200,800)

#graph.GetXaxis().SetRangeUser(171, 176)

graph.SetMarkerStyle(8)
graph.Draw("ap")

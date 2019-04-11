#!/usr/bin/env python
from ROOT import *
from array import array

f = TFile.Open("bin1_cut220_mtTemplatesForCH.root")
h = f.Get("rec_ptll/ttactual1725").Clone()
h.SetDirectory(0)
f.Close()

nq = 2 
xq = array('d', [0.] * (nq))  # Input range from [0 ... 1]
yq = array('d', [0.] * (nq))  # Output quantiles
for i in range(nq): xq[i] = float(i+1)/float(nq)
h.GetQuantiles(nq, yq, xq)

bins = [0] + [i for i in yq]
rebinned = h.Rebin(nq, "rebinned", array('d',bins))
#rebinned.Scale(1., "width")
#rebinned.Draw("hist")
rebinned.Draw("e1")

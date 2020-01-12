#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
import sys
from array import array

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="tt_bin10_morph2D.root")

args = parser.parse_args()

f = TFile.Open(args.inF)
g = f.Get("tt_ptll_graph2D")

g.Draw()
g.GetXaxis().SetTitle("m_{t} [GeV]")
g.GetXaxis().SetTitleOffset(1.5)
g.GetYaxis().SetTitle("p_{T}(ll) [GeV]")
g.GetYaxis().SetTitleOffset(1.5)

g.GetXaxis().SetRangeUser(166, 179)
g.GetYaxis().SetRangeUser(-0.5, 200.5)

g.Draw("err")

for p in xrange(g.GetN()):
    x = Double()
    y = Double()
    z = Double()
    g.GetPoint(p, x, y, z)
    print "(%.1f, %.1f, %.1f)" % (x,y,z)


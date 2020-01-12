#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="higgsCombine_bin9_grid.MultiDimFit.mH125.root")
parser.add_argument("-b", "--bin", type=int, default=2)
args = parser.parse_args()

f = TFile.Open(args.inF)
t = f.Get("limit")

t.Draw("2*deltaNLL:MT/10.", "2*deltaNLL < 1.2")


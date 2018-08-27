#!/usr/bin/env python
from ROOT import *
from array import array
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-o", dest="outF", default="interp.root", help="Output root file")
args = parser.parse_args()

inputs = ["linear", "cspline", "polynomial", "akima"]
masses = range(1665, 1795, 10)
morph = {}
for interp in inputs:
    morph[interp] = {}
    f = TFile.Open("debug_%s.root" % interp, "read")
    for m in masses:
        morph[interp][m] = f.Get("emu_signal_morph/morph_point_%d" % m)
        morph[interp][m].SetDirectory(0)

    f.Close()

binG = {}
nbins = morph[inputs[0]][masses[0]].GetNbinsX()
for interp in inputs:
    binG[interp] = {}
    for b in xrange(1, nbins+1):
        binG[interp][b] = TGraphErrors(len(masses), array('d', masses), array('d', [morph[interp][m].GetBinContent(b) for m in masses]), array('d', [0.] * len(masses)), array('d', [morph[interp][m].GetBinError(b) for m in masses]))
        binG[interp][b].SetName("bin_%d_%s" % (b, interp))
        binG[interp][b].SetTitle("Bin %d  %s" % (b, interp))
        binG[interp][b].GetXaxis().SetTitle("m_{t}")
        binG[interp][b].GetYaxis().SetTitle("Bin %d Contents" % b)

outF = TFile.Open(args.outF, "recreate")
for interp in inputs:
    for b in xrange(1, nbins+1):
        binG[interp][b].Write()

outF.Close()

#!/usr/bin/env python
from ROOT import *
from array import array
import os
import sys
from argparse import ArgumentParser
from pprint import pprint

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input templates")
parser.add_argument("--obs", default="rec_ptll", help="reco lvl + observable")
parser.add_argument("-o", "--outF", default="binErrors.root", help="output root file")
args = parser.parse_args()

f = TFile.Open(args.inF, 'read')
masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]

actualH = {}
rate = {}
for m in masses:
    actualH[m] = f.Get("%s/ttactual%d" % (args.obs, m)).Clone()
    actualH[m].SetDirectory(0)
    rate[m] = actualH[m].Integral()
    actualH[m].Scale(1./actualH[m].Integral())

f.Close()


g = {}
g2D = {}
centers = {}
contents = {}
errors = {}


for m in masses:
    centers[m] = []
    contents[m] = []
    errors[m] = []

    nBins = actualH[m].GetNbinsX()
    for b in xrange(1,nBins+1):
        centers[m].append(actualH[m].GetBinCenter(b))
        contents[m].append(actualH[m].GetBinContent(b))
        errors[m].append(actualH[m].GetBinError(b))

    g[m] = TGraph(nBins, array('d',contents[m]), array('d',errors[m]))
    g[m].SetName("errors_mt%d" % m)
    g[m].SetTitle("Bin Error vs Bin Content  m_{t} = %.1f" % (m/10.))
    g[m].GetXaxis().SetTitle("Bin content")
    g[m].GetXaxis().SetTitleOffset(1.1)
    g[m].GetYaxis().SetTitle("Bin error")
    g[m].GetYaxis().SetTitleOffset(1.1)
    g[m].SetMarkerStyle(22)

    
    g2D[m] = TGraph2D(nBins, array('d', centers[m]), array('d', contents[m]), array('d', errors[m]))
    g2D[m].SetName("g2D_errors_mt%d" % m)
    g2D[m].GetXaxis().SetTitle("p_{T}(ll)")
    g2D[m].GetYaxis().SetTitle("Bin Content")
    g2D[m].GetZaxis().SetTitle("Bin Error")


f = TFile.Open(args.outF, "recreate")
for m in masses:
    g[m].Write()
    g2D[m].Write()

f.Close()


# Fit nominal mass graph
g[1725].Draw("ap")
fit = TF1("sqrt", "[0] * x ^ (0.5)")
g[1725].Fit("sqrt")
print "\n"

# Get parameter [0]
const = fit.GetParameter(0)
print "\tp0 =", const
print "Neff**-0.5 =", actualH[1725].GetEffectiveEntries() ** -0.5

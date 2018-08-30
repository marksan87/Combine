#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", dest="inF", default="templates.root", help="Input debug root file from CombineHarvester with morphed templates")
parser.add_argument("-o", dest="outF", default="morph2D.root", help="Output root file")
args = parser.parse_args()
if args.outF.find(".root") < 0:
    args.outF += ".root"

f = TFile.Open(args.inF)
deltaM = 1
masses = range(1665, 1786, deltaM)
#masses = range(1665, 1795, 10)
morph = {}

for m in masses:
    morph[m] = f.Get("rec_ptll/tt%d" % m)
    #morph[m] = f.Get("rec_ptll/tW%d" % m)
    morph[m].SetDirectory(0)

#h = TH2D("morph_points", "Morphed Points", len(masses), masses[0] / 10.0, masses[-1] / 10.0, 50, 30, 180)
h = TH2D("morph_points", "Morphed p_{T}(ll) Templates", len(masses), (masses[0] - 0.5*deltaM) / 10.0, (masses[-1] + 0.5*deltaM)/10.0, morph[masses[0]].GetNbinsX(), morph[masses[0]].GetXaxis().GetXmin(), morph[masses[0]].GetXaxis().GetXmax())
h.GetXaxis().SetTitle("m_{t} [GeV]")
h.GetYaxis().SetTitle("p_{T}(ll) [GeV]")
h.GetZaxis().SetTitle("Events")

h.GetXaxis().SetTitleOffset(1.6)
h.GetYaxis().SetTitleOffset(1.8)
h.GetZaxis().SetTitleOffset(1.5)

for m in masses:
    for b in xrange(1, morph[m].GetNbinsX() + 1):
        h.Fill(m / 10.0, morph[m].GetBinCenter(b), morph[m].GetBinContent(b))

h.Draw("lego2z")
outF = TFile.Open(args.outF, "recreate")
h.Write()
for m in masses:
    morph[m].Write()

outF.Close()

print "Results saved in %s" % args.outF

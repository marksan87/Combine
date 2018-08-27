#!/usr/bin/env python
from ROOT import *
from array import array
from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument("-i", dest="inDir", default="RootFiles", help="Input template directory")
parser.add_argument("-obs", dest="obs", default="pt_ll", help="Observable")
parser.add_argument("-f", dest="func", default="pol3", help="interpolation function")
parser.add_argument("-d", dest="deltaM", default=1, help="mass spacing (in 0.1 GeV)", type=int)
parser.add_argument("--plots", action="store_true", default=False, help="Create per-bin and difference plots")
parser.add_argument("-o", dest="outF", default="morphed_templates.root", help="Output root file")
args = parser.parse_args()

if args.inDir[-1] == "/": 
    args.inDir = args.inDir[:-1]

masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
if args.deltaM < 1e-10:
    print "Mass spacing too small! Defaulting to 0.1 GeV"
    args.deltaM = 1
    
morphed_masses = range(masses[0], masses[-1]+args.deltaM, args.deltaM) 

actual = {}

for m in masses:
    f = TFile.Open("%s/mc_TT_mt%d.root" % (args.inDir, m), "read")
    actual[m] = f.Get(args.obs).Clone("signal%d" % m)
    actual[m].SetDirectory(0)
    f.Close()

#f = TFile.Open(args.inF, "read")
#for m in masses:
#    actual[m] = f.Get("emu/signal%d" % m)
#    actual[m].SetDirectory(0)

outF = TFile.Open(args.outF, "recreate")


rates = {}
# Normalize
for m in masses:
    rates[m] = actual[m].Integral()
    actual[m].Scale(1.0/actual[m].Integral())


binG = {}
nbins = actual[masses[0]].GetNbinsX()
xmin = actual[masses[0]].GetXaxis().GetXmin()
xmax = actual[masses[0]].GetXaxis().GetXmax()

# Morphed templates
morph = {}
for m in morphed_masses:
    morph[m] = TH1F("morph_mt%d" % m, "Morph Point  m_{t} = %.1f" % (float(m)/10.), nbins, xmin, xmax)

diff = {}  # Difference between actual and morphed templates

binG = {}
binMorphG = {}
fit = {}
fitFunc = {}



for b in xrange(1, nbins+1):
    binG[b] = TGraphErrors(len(masses), array('d', [m/10.0 for m in masses]), array('d', [actual[m].GetBinContent(b) for m in masses]), array('d', [0.] * len(masses)), array('d', [actual[m].GetBinError(b) for m in masses]))
    binG[b].SetName("bin_%d" % b)
    binG[b].SetTitle("Bin %d" % b)
    binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
    binG[b].GetYaxis().SetTitle("Entries")
    binG[b].GetYaxis().SetTitleOffset(1.3)
    binG[b].Write()
    
    fit[b] = binG[b].Fit(args.func, "S")
    fitFunc[b] = binG[b].GetFunction(args.func)
    errors = array('d', [0.] * len(morphed_masses))
    fit[b].GetConfidenceIntervals(len(morphed_masses), 1, 1, array('d', [m/10.0 for m in morphed_masses]), errors, 2./3., False)
    for i,m in enumerate(morphed_masses):
        if b == 1:
            print "m = %.1f\tfit = %f" % (m/10.0, fitFunc[b].Eval(m/10.0))
        morph[m].SetBinContent(b, fitFunc[b].Eval(m/10.0))
        morph[m].SetBinError(b, errors[i]) 

for b in xrange(1, nbins+1):
    binMorphG[b] = TGraphErrors(len(morphed_masses), array('d', [m/10.0 for m in morphed_masses]), array('d', [morph[m].GetBinContent(b) for m in morphed_masses]), array('d', [0.] * len(morphed_masses)), array('d', [morph[m].GetBinError(b) for m in morphed_masses]))
    binMorphG[b].SetName("morphed_bin_%d" % b)
    binMorphG[b].SetTitle("Morphed Bin %d" % b)
    binMorphG[b].GetXaxis().SetTitle("m_{t} [GeV]")
    binMorphG[b].GetYaxis().SetTitle("Entries")
    binMorphG[b].GetYaxis().SetTitleOffset(1.3)
    binMorphG[b].Write()


for m in masses:
    diff[m] = actual[m].Clone("diff_mt%d" % m)
    diff[m].SetDirectory(0)
    diff[m].Add(morph[m], -1)
    diff[m].Divide(actual[m])
    #diff[m].Scale(rates[m])
    diff[m].SetTitle("(actual - morphed) / actual   m_{t} = %.1f" % (m/10.))
    diff[m].GetYaxis().SetTitle("(actual - morphed) / actual")
    diff[m].GetXaxis().SetTitleOffset(1.2)

for m in morphed_masses:
    morph[m].Scale(rates[1725])


for m in masses:
    actual[m].Scale(rates[m])
    actual[m].Write()
    diff[m].Write()
for m in morphed_masses:
    morph[m].Write()
outF.Close()

if args.plots:
    gROOT.SetBatch(True)
    outDir = "diff"
    os.system("mkdir -p %s" % outDir)
    c = TCanvas("foo","bar", 2000, 1200)
    line = TF1("line", "0", -99999, 99999)
    gStyle.SetOptStat(0)
    for m in masses:
        diff[m].GetYaxis().SetRangeUser(-0.3, 0.3)
        diff[m].Draw()
        c.Update()
        tline = TLine(c.GetUxmin(), 0, c.GetUxmax(), 0)
        tline.Draw("SAME")
        c.SaveAs("%s/diff_mt%d.png" % (outDir, m))

    os.system("mkdir -p bins")

    for b in xrange(1, nbins+1):
        l = TLegend(0.425, 0.8, 0.575, 0.9)
        binMorphG[b].SetFillColor(0)
        binG[b].SetFillColor(0)
        binG[b].SetLineWidth(2)
        binMorphG[b].Draw("ALP")
        binG[b].SetLineColor(kBlue)
        binG[b].Draw("LP SAME")
        fit = binG[b].Fit(args.func, "S")
        l.AddEntry(binMorphG[b], "Morphed")
        l.AddEntry(binG[b], "Actual")
        l.SetFillStyle(0)
        #l.AddEntry(fit.Get(), "Fit")
        l.Draw("SAME")
        c.SaveAs("bins/morph_bin_%d.png" % b)

print "Results saved in %s" % args.outF

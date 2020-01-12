#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
import sys
from pprint import pprint
from glob import glob

systematics=[\
    "fit",
    "BkgNorm",
    "Lumi",
    "BTagSF",
    "CRGluon",
    "CRQCD",
    "CRerdON",
    "EleIDEff",
    "EleRecoEff",
    "EleScale",
    "EleSmear",
    "JEC",
    "JER",
    "MuIDEff",
    "MuIsoEff",
    "MuScale",
    "MuTrackEff",
    "Pdf",
    "Q2",
    "TrigEff",
    "UE",
    "fsr",
    "hdamp",
    "isr",
    "pileup",
    "toppt",
    "DS",
    ]


gROOT.SetBatch(True)
gStyle.SetStatFormat("3.3f")

parser = ArgumentParser()
parser.add_argument("-i", "--inDir", default="results_pseudodataToys/condor_4403756", help="input file or directory containing toy root files")
parser.add_argument("-o", "--outDir", default="pseudodataToys_plots", help="output directory")
args = parser.parse_args()

inDir = args.inDir
outDir = args.outDir
if inDir[-1] == "/": inDir = inDir[:-1]
if outDir[-1] == "/": outDir = outDir[:-1]

os.system("mkdir -p %s" % outDir)
#if inDir.find(".root") > -1:
#    # Single file
#    toyfiles = [inDir]
#else:
#    # Directory
#    toyfiles = glob("%s/higgsCombine*mH125*.root" % inDir)

hists = {}

print "Loading toy data from %s" % inDir
for syst in systematics:
    toyDir = "%s/%s_toyFits" % (inDir,syst)
    toyF = glob("%s/*%s*.root" % (toyDir,syst if syst != "fit" else "initialFit"))[0]
    print "Syst %s from file: %s" % (syst,toyF)
    f = TFile.Open(toyF, "read")
    t = f.Get("limit")
    t.Draw("MT/10>>%sH" % syst, "quantileExpected == -1")
    exec("%sH.SetDirectory(0)" % syst)
    exec("%sH.SetTitle('%s Toys')" % (syst,syst))
    exec("hists[syst] = %sH" % syst)

    f.Close()

print "Found %d toy files" % len(systematics)

c = TCanvas("c","c",1200,800)
for syst,h in hists.items():
    h.Draw("hist")
    h.GetXaxis().SetTitle("m_{t} [GeV]")
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitle("Entries")
    h.GetYaxis().SetTitleOffset(1.3)
    
    c.SaveAs("%s/%s_toys.png" % (outDir,syst))


f = TFile.Open("%s/hists.root" % outDir, "recreate") 
for syst,h in hists.items():
    h.Write()

f.Close()

print "Histograms saved to %s" % outDir

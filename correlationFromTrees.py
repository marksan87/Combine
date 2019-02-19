#!/usr/bin/env python
import ROOT
from ROOT import TFile,TCanvas,gROOT,gStyle,TH2D
from argparse import ArgumentParser
from array import array
import os
from time import sleep
from pprint import pprint
import pickle

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "ptneg":"p_{T}(l^{-})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

gStyle.SetOptStat(0)

def makeHist2D(tree, config, xobs, yobs):
    name = "%s_vs_%s" % (yobs[4:], xobs[4:])
    hist = TH2D(name, "%s  vs  %s" % (config[yobs]["title"],config[xobs]["title"]), config[xobs]["nbins"], config[xobs]["min"], config[xobs]["max"], config[yobs]["nbins"], config[yobs]["min"], config[yobs]["max"])
    hist.GetXaxis().SetTitle(config[xobs]["title"])
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitle(config[yobs]["title"])
    hist.GetYaxis().SetTitleOffset(1.4)
    tree.Draw("%s:%s >> %s" % (yobs,xobs,name), "weight")
    hist.SetDirectory(0)
    return hist


def plotCorrelations(inDir, configF, outDir): 
    # Load config file for distributions and bin min,max,count
    config = {}
    with open(configF, "r") as f:
        for i,line in enumerate(f):
            #print "line %d: len = %d" %(i,len(line.strip()))
            if (i == 0 and line[0] == "#") or len(line.strip()) == 0: continue
            l = line.split()
            config[l[0]] = {"nbins":int(l[1]), "min":float(l[2]), "max":float(l[3]), "title":obsTitle[l[0][4:]]}

    print "\nFound observables from %s:\n" % args.config
    for obs,vals in config.items():
        print "%s\t%d bins from %.1f-%.1f GeV  (bin width = %.1f GeV)" % (obs if len(obs) >= 8 else obs+"\t", vals["nbins"], vals["min"], vals["max"], (vals["max"] - vals["min"]) / float(vals["nbins"])) 
    print "\n"

    ttmasses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
    print "Processing ttrees from %s/mc_TT_mt1725.root" % inDir

    f = TFile.Open("%s/mc_TT_mt1725.root" % inDir, "read")
    tree = f.Get("goodEvents")

    histos = []
    obsList = ["rec_ptll", "rec_ptpos", "rec_ptneg", "rec_Mll", "rec_Epos", "rec_ptp_ptm", "rec_Ep_Em"]
    for xobs in obsList:
        for yobs in obsList[1:]:
            histos.append(makeHist2D(tree=tree, config=config, xobs=xobs, yobs=yobs))
        obsList = obsList[1:]
    f.Close()

    os.system("mkdir -p %s" % outDir)
    c = TCanvas("c","c",1200,1200)
    c.Clear()
    c.SetRightMargin(0.13)
    f = TFile.Open("%s/plots.root" % outDir, "recreate")
    for h in histos:
        h.Draw("colz")
        c.SaveAs("%s/%s.png" % (outDir,h.GetName()))
        h.Write()

    f.Close()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", dest="inDir", default="plots2018", help="Input directory containing ttrees")
    parser.add_argument("-o", dest="outDir", default="correlations", help="Output directory")
    parser.add_argument("-c", dest="config", default="config.txt", help="Config file specifying observables and bin ranges/sizes")
    args = parser.parse_args()

    if args.inDir[-1] == "/": args.inDir = args.inDir[:-1]
    if args.outDir[-1] == "/": args.outDir = args.outDir[:-1]
   
    plotCorrelations(inDir = args.inDir, configF = args.config, outDir = args.outDir)



#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
from glob import glob
from pprint import pprint
from array import array

gROOT.SetBatch(True)

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="", help="CH datacard")
parser.add_argument("-o", "--outDir", default="")#"grid_scans")
parser.add_argument("--minmt", type=float, default=171.0)
parser.add_argument("--maxmt", type=float, default=174.0)
parser.add_argument("--log", default="combine.log", help ="log file with output of combine commands")
parser.add_argument("--noFreeze", action="store_true", default=False, help="float other nuisances when scanning one np")
parser.add_argument("-p", "--points", default=100)
parser.add_argument("-c", "--cut", type=float, default=1.3, help="2*deltaNLL to cut on")
parser.add_argument("-d", "--data_obs", action="store_true", default=False, help="use data_obs instead of asimov dataset")
args = parser.parse_args()

minmt = int(10*args.minmt)
maxmt = int(10*args.maxmt)

freezePars = not args.noFreeze  # Controls whether to freeze other nps when scanning one np
if freezePars:
    print "Other nuisances fixed to best-fit values when scanning a np"
else:
    print "All other nuisances floated while scanning a np"


useAsimov = not args.data_obs
if useAsimov:
    print "Using Asimov dataset"
else:
    print "Using data_obs"

if args.outDir == "":
    args.outDir = args.inF.replace(".txt", "").replace(".root","")
    args.outDir = args.outDir.replace("_cardMorph", "_gridScan")

outDir = args.outDir
cutNLL = args.cut

print "2*deltaNLL cut at", cutNLL

command = "mkdir -p %s; " % outDir

cardRoot = args.inF.replace(".txt",".root")
if args.inF.find(".txt") >= 0:
    # Create root workspace
    command += "text2workspace.py %s; " % args.inF
    
command += "cp %s %s/; rm -f %s/%s" % (cardRoot, outDir, outDir, args.log)
os.system(command)


f = TFile.Open("%s/%s" % (outDir, cardRoot))
w = f.Get("w")
nset = w.set("nuisances")
nuisances = []
niter = nset.createIterator()

for i in range(len(nset)):
    nuisances.append(niter.Next().getPlotLabel())

#nuisances = ['bin2']
print "Found np's:", nuisances
f.Close()
#pprint(nuisances)

# saveSpecifiedNuis causes seg fault
#options = "--saveSpecifiedNuis all --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"
#options = "--setRobustFitTolerance 0.2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"
options = "--floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"

# Run impact grid scans
for n in nuisances:
    freezeArgs = ""
    if freezePars:
        # Comma separated list of all additional nps
        #freezeArgs = " --freezeParameters "
        freezeArgs = " --freezeParameters "
        for freezeMe in nuisances:
            if freezeMe != n:
                freezeArgs += "%s," % freezeMe

        # Remove trailing comma
        freezeArgs = freezeArgs[:-1]

    combineCmd = "combine -M MultiDimFit -n _%s_grid --algo grid --points %d%s --redefineSignalPOIs MT,r -P %s %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=1725 -d %s |& tee -a %s" % (n, args.points, freezeArgs, n, options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), cardRoot, args.log)
    #combineCmd = "combine -M MultiDimFit -n _%s_grid --algo grid --points %d --redefineSignalPOIs MT,r,%s -P %s %s --setParameterRanges MT=%d,%d:r=0,2 %s -m 125 --setParameters MT=1725 -d %s |& tee -a %s" % (n, args.points, n, n, options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), cardRoot, args.log)
    os.system("pushd %s; echo '%s'; %s" % (outDir, combineCmd, combineCmd))
    #os.system("pushd %s; combine -M MultiDimFit -n _%s_grid --algo grid --points %d --redefineSignalPOIs MT,r -P %s %s --setParameterRanges MT=%d,%d:r=0,2 %s -m 125 --setParameters MT=1725 -d %s |& tee -a %s" % (outDir, n, args.points, n, options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), cardRoot , args.log) )

# NLL vs parameter
paramG = {}
paramG_cut = {}

# NLL vs MT for each parameter
mtG = {}
mtG_cut = {}

c = TCanvas("c","C",1200,800)

# Get graphs
for n in nuisances:
    f = TFile.Open("%s/higgsCombine_%s_grid.MultiDimFit.mH125.root" % (outDir, n) )
    t = f.Get("limit")
    
    # All entries
    nll = []
    mt = []
    par = []

    # Entries passing 2*deltaNLL cut
    nll_cut = []
    mt_cut = []
    par_cut = []

    for i in range(1, t.GetEntriesFast()+1):
        t.GetEntry(i)
        nll.append(2*t.deltaNLL)
        mt.append(t.MT / 10.)
        exec("par.append(t.%s)" % n)
        
        if abs(nll[-1]) < cutNLL:
            nll_cut.append(2*t.deltaNLL)
            mt_cut.append(t.MT / 10.)
            exec("par_cut.append(t.%s)" % n)

    paramG[n] = TGraph(len(nll), array('d',par), array('d',nll))
    paramG[n].SetName("nll_vs_%s" % n)
    paramG[n].SetTitle("NLL vs %s" % n)
    paramG[n].GetXaxis().SetTitle(n)
    paramG[n].GetYaxis().SetTitle("2*#Delta{NLL}")
    paramG[n].SetMarkerStyle(22)

    paramG_cut[n] = TGraph(len(nll_cut), array('d',par_cut), array('d',nll_cut))
    paramG_cut[n].SetName("nll_vs_%s" % n)
    paramG_cut[n].SetTitle("NLL vs %s" % n)
    paramG_cut[n].GetXaxis().SetTitle(n)
    paramG_cut[n].GetYaxis().SetTitle("2*#Delta{NLL}")
    paramG_cut[n].SetMarkerStyle(22)
    
    mtG[n] = TGraph(len(nll), array('d', mt), array('d',nll))
    mtG[n].SetName("%s_nll_vs_mt" % n)
    mtG[n].SetTitle("%s  NLL vs m_{t}" % n)
    mtG[n].GetXaxis().SetTitle("m_{t} [GeV]")
    mtG[n].GetYaxis().SetTitle("2*#Delta{NLL}")
    mtG[n].SetMarkerStyle(22)
    
    mtG_cut[n] = TGraph(len(nll_cut), array('d', mt_cut), array('d',nll_cut))
    mtG_cut[n].SetName("cut_%s_nll_vs_mt" % n)
    mtG_cut[n].SetTitle("%s  NLL vs m_{t}" % n)
    mtG_cut[n].GetXaxis().SetTitle("m_{t} [GeV]")
    mtG_cut[n].GetYaxis().SetTitle("2*#Delta{NLL}")
    mtG_cut[n].SetMarkerStyle(22)

    paramG[n].Draw("alp")
    c.SaveAs("%s/nll_vs_%s.png" % (outDir,n))

    mtG[n].Draw("alp")
    c.SaveAs("%s/%s_nll_vs_mt.png" % (outDir,n))
    
    paramG_cut[n].Draw("alp")
    c.SaveAs("%s/cut_nll_vs_%s.png" % (outDir,n))

    mtG_cut[n].Draw("alp")
    c.SaveAs("%s/cut_%s_nll_vs_mt.png" % (outDir,n))

# Write output file
f = TFile.Open("%s/%s.root" %(outDir, outDir), "recreate")
for n in nuisances:
    paramG[n].Write()
    paramG_cut[n].Write()
    mtG[n].Write()
    mtG_cut[n].Write()
f.Close()

gROOT.SetBatch(False)
#c = TCanvas("c2","c2",1200,800)



#t = f.Get("limit")

#t.Draw("2*deltaNLL:MT/10.", "2*deltaNLL < 1.2")


#os.system("combine --help")

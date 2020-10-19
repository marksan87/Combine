#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
import sys
from glob import glob
from pprint import pprint
from array import array

gROOT.SetBatch(True)

oneSidedSysts = ["CRerdON", "CRGluon", "CRQCD", "amcanlo", "madgraph", "herwigpp", "toppt", "DS"]

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="", help="CH datacard")
parser.add_argument("-o", "--outDir", default="")#"grid_scans")
parser.add_argument("-m", "--mt", type=float, default=172.5, help="initial mt for generating asimov dataset")
parser.add_argument("--minmt", type=float, default=166.5)
parser.add_argument("--maxmt", type=float, default=178.5)
parser.add_argument("-f", "--format", "--formats", dest="formats", nargs="+", default=["png"], help="output file formats")
parser.add_argument("--log", default="combine.log", help ="log file with output of combine commands")
parser.add_argument("--noFreeze", action="store_true", default=False, help="float other nuisances when scanning one np")
parser.add_argument("--robustHesse", action="store_true", default=False, help="use option --robustHesse 1")
parser.add_argument("--nosysts", action="store_true", default=False, help="don't include systematic likelihood scans")
parser.add_argument("-p", "--points", type=int, default=200)
parser.add_argument("-c", "--cut", type=float, default=1.3, help="2*deltaNLL to cut on")
parser.add_argument("--noconst", dest="noconst", action="store_true", default=False, help="absolute likelihood")
parser.add_argument("-d", "--data_obs", action="store_true", default=False, help="use data_obs instead of asimov dataset")
parser.add_argument("-v", "--verbosity", type=int, default=0, help="verbosity level")
parser.add_argument("-q", "--quiet", action="store_true", default=False, help="no terminal output (will still log with specified verbosity level")
args = parser.parse_args()

useSysts = not args.nosysts
useRobustHesse = args.robustHesse
if useRobustHesse:
    print "Using --robustHesse 1"
formats = [fmt.strip(".") for fmt in args.formats]
print "Saving plots as:", formats
freezePars = not args.noFreeze  # Controls whether to freeze other nps when scanning one np
if not args.quiet:
    if freezePars:
        print "Other nuisances fixed to best-fit values when scanning a np"
    else:
        print "All other nuisances floated while scanning a np"


useAsimov = not args.data_obs
if not args.quiet:
    if useAsimov:
        print "Using Asimov dataset"
    else:
        print "Using data_obs"

outDir = args.outDir
inF = args.inF

if args.outDir == "":
    outDir = os.path.basename(inF).replace(".txt", "").replace(".root","")
    outDir = outDir.replace("_cardMorph", "_gridScans")

os.system("mkdir -p %s" % outDir)
if args.outDir != "":
    # Copy card data into output directory
    os.system("cp %s %s/%s" % (inF,outDir,os.path.basename(inF)))
    #inF = os.path.basename(inF)

cutNLL = args.cut

if not args.quiet: print "2*deltaNLL cut at", cutNLL

command = ""
cardRoot = inF.replace(".txt",".root")
if inF.find(".txt") >= 0:
    # Create root workspace
    command += "text2workspace.py %s; " % inF
    
command += "cp %s %s/; rm -f %s/%s" % (cardRoot, outDir, outDir, args.log)
os.system(command)

cardRoot = os.path.basename(cardRoot)

loggingCmd = " |& tee -a %s" % args.log if not args.quiet else " &>> %s" % args.log
silenceOutput = " &>> %s" % args.log if args.quiet else ""

f = TFile.Open("%s/%s" % (outDir, os.path.basename(cardRoot)))


w = f.Get("w")
MTvar = w.var("MT")
precision = len(str(int(MTvar.getVal()))) - 3      # Decimal places of precision

print "precision:", precision

mt = int((10**precision)*args.mt)
minmt = int((10**precision)*args.minmt)
maxmt = int((10**precision)*args.maxmt)



nset = w.set("nuisances")
nuisances = []
if useSysts:
    try:
        niter = nset.createIterator()
        for i in range(len(nset)):
            nuisances.append(niter.Next().getPlotLabel())
    except ReferenceError:
        # No nuisances found!
        niter = None
    if not args.quiet: print "Found np's:", nuisances
f.Close()
#pprint(nuisances)

# saveSpecifiedNuis causes seg fault
#options = "--saveSpecifiedNuis all --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"
#options = "--setRobustFitTolerance 0.2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"

options = "--floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1 --cminDefaultMinimizerStrategy=1 --cminFallbackAlgo Minuit2,0"
if args.noconst:
    options += " --saveNLL --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
#options = "--saveNLL --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1 --cminDefaultMinimizerStrategy=1 --cminFallbackAlgo Minuit2,0"
if useRobustHesse:
    options += " --robustHesse 1"
#options = "--floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1 --robustHesse 1"


# Initial fit
combineCmd = "combine -M MultiDimFit -n _bestfit --algo singles --redefineSignalPOIs MT,r %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d -d %s -v %d --saveWorkspace %s" % (options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, cardRoot, args.verbosity, loggingCmd)
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

# Run impact grid scans
for n in ["MT"]+nuisances:
    freezeArgs = ""
    if freezePars and n != "MT" and len(nuisances) > 1:
        
        # Comma separated list of all additional nps
        #freezeArgs = " --freezeParameters "
        freezeArgs = " --freezeParameters "
        for freezeMe in nuisances:
            if freezeMe != n:
                freezeArgs += "%s," % freezeMe

        # Remove trailing comma
        freezeArgs = freezeArgs[:-1]

    combineCmd = "combine -M MultiDimFit -n _%s_grid --algo grid --points %d%s --redefineSignalPOIs MT,r -P %s %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d higgsCombine_bestfit.MultiDimFit.mH125.root --snapshotName MultiDimFit -v %d %s" % (n, args.points, freezeArgs, n, options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, args.verbosity, loggingCmd)
    os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

# NLL vs parameter
paramG = {}
paramG_cut = {}

# NLL vs MT for each parameter
mtG = {}
mtG_cut = {}

c = TCanvas("c","C",1200,800)

# Get graphs
for n in ["MT"]+nuisances:
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
        exec("pval = t.%s" % n)
        if n in oneSidedSysts and pval < 0: continue
            #exec("if t.%s < 0: continue" % n)
        
        if args.noconst:
            nll.append(2*(t.deltaNLL + t.nll + t.nll0))
        else:
            nll.append(2*t.deltaNLL)
        mt.append(t.MT / 10.)
        exec("par.append(t.%s)" % n)
        
        if abs(nll[-1]) < cutNLL:
            nll_cut.append(nll[-1])
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

    if n != "MT":
        paramG[n].Draw("alp")
        for fmt in formats:
            c.SaveAs("%s/nll_vs_%s.%s" % (outDir,n,fmt))
        
        paramG_cut[n].Draw("alp")
        for fmt in formats:
            c.SaveAs("%s/cut_nll_vs_%s.%s" % (outDir,n,fmt))

    mtG[n].Draw("alp")
    for fmt in formats:
        c.SaveAs("%s/%s_nll_vs_mt.%s" % (outDir,n,fmt))
    
    mtG_cut[n].Draw("alp")
    for fmt in formats:
        c.SaveAs("%s/cut_%s_nll_vs_mt.%s" % (outDir,n,fmt))

# Write output file
f = TFile.Open("%s/%s.root" %(outDir, outDir), "recreate")
for n in ["MT"]+nuisances:
    paramG[n].Write()
    paramG_cut[n].Write()
    mtG[n].Write()
    mtG_cut[n].Write()
f.Close()

gROOT.SetBatch(False)

print "\nOutput saved to %s\n" % outDir
#c = TCanvas("c2","c2",1200,800)



#t = f.Get("limit")

#t.Draw("2*deltaNLL:MT/10.", "2*deltaNLL < 1.2")


#os.system("combine --help")

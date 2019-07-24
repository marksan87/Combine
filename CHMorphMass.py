#!/usr/bin/env python
from __future__ import division
import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT
from argparse import ArgumentParser
import sys
from ROOT.TMath import Ceil,Log10
from pprint import pprint

systematics = {\
    "Lumi":"lnN", \
    "BkgNorm":"lnN", \
    "EleIDEff":"shape", \
    "EleRecoEff":"shape", \
    "EleScale":"shape", \
    "EleSmear":"shape", \
    "MuIDEff":"shape", \
    "MuIsoEff":"shape", \
    "MuTrackEff":"shape", \
    "MuScale":"shape", \
    "TrigEff":"shape", \
    "BTagSF":"shape", \
    "JEC":"shape", \
    "JER":"shape", \
    "pileup":"shape", \
    "Pdf":"shape", \
    "isr":"shape", \
    "fsr":"shape", \
    "toppt":"shape", \
    "Q2":"shape", \
    "hdamp":"shape", \
    "UE":"shape", \
    "CRerdON":"shape", \
    "CRGluon":"shape", \
    "CRQCD":"shape", \
    "DS":"shape", \
#    "amcanlo":"shape", \
#    "madgraph":"shape", \
#    "herwigpp":"shape", \
    }

experimentalSysts = ["pileup", "Lumi", "BTagSF", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER",'BkgNorm' ]
theorySysts = ["toppt", "Q2", "isr", "fsr", "Pdf", "hdamp", "UE", "CRerdON", "CRGluon", "CRQCD", "amcanlo", "madgraph", "herwigpp", "DS" ]

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input template root file")
parser.add_argument("-o", "--outF", default="", help="append name to output datacard and rootfile")
parser.add_argument("--minmt", type=float, default=166.5, help="min MT for morph mass range")
parser.add_argument("--maxmt", type=float, default=178.5, help="max MT for morph mass range")
parser.add_argument("--noBinStats", action="store_true", default=False, help="don't include binwise statistics for signal")
parser.add_argument("--deltaMT", type=float, default=0.1, help="MT increments for morphing range")
parser.add_argument("--precision", type=int, default=1, help="decimal places of precision for morphed templates")
parser.add_argument("--splineInterp", default="CSPLINE", help="CH roofit morphing interpolation mode")
parser.add_argument("--systs", default="", nargs="+", choices=(["None","none"] + systematics.keys()), help="ONLY plot these systematics (or 'none' for no systematics)") 
parser.add_argument("--obs", default="ptll", help="observable")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reco level")
parser.add_argument("-v", "--verbosity", type=int, default=3, help="verbosity level (0+)")
args = parser.parse_args()


addBinStats = not args.noBinStats

if args.outF != "":
    outF = args.outF + "_"
elif args.inF.find("mtTemplatesForCH.root") >= 0:
    outF = args.inF.replace("mtTemplatesForCH.root","")
else:
    outF = ""

if args.precision < 0:
    print "Cannot have %d decimal places! Defaulting to 1" % args.precision
    args.precision = 1

elif args.deltaMT < 10**(-args.precision):
    # Increase precision to leading decimal value in deltaMT
    args.precision = int(Ceil(Log10(1./args.deltaMT)))
    print "Using precision: %d" % args.precision

if args.systs == "":
    print "All systematics included"
elif "none" in args.systs or "None" in args.systs:
    systematics = {}
    print "No systematics included"
else:
    systsToRemove = systematics.keys()
    for _s in args.systs:
        systsToRemove.remove(_s)
    for _s in systsToRemove:
        systematics.pop(_s, None)
    print "Using systematics: ", systematics.keys()


precision = args.precision
decimalScaling = 10**precision


outCardName = "%scardMorph.txt" % (outF)
outRootName = "%soutputfileMorph.root" % (outF)

#outCardName = "%scardMorph.txt" % ("%s_" % args.outF if args.outF != "" else "")
#outRootName = "%soutputfileMorph.root" % ("%s_" % args.outF if args.outF != "" else "")
cb = ch.CombineHarvester()
cb.SetVerbosity(args.verbosity)

"""
masses = [  \
          '1665',
          '1695',
          '1715',
          '1725',
          '1735',
          '1755',
          '1785' \
          ]
"""
#deltaM = 1
#masses = [str(1665 + i * deltaM) for i in xrange(121)]
#masses = [str(1695 + i * deltaM) for i in xrange(61)]
masses = range(int(args.minmt*decimalScaling), int((args.maxmt+args.deltaMT)*decimalScaling), int(args.deltaMT*decimalScaling))
masses = [str(m) for m in masses]


exec('minstr = "%.' + str(precision) + 'f" % (args.minmt)')
exec('maxstr = "%.' + str(precision) + 'f" % (args.maxmt)')
    
print "\nUsing %d morphed templates in range %s to %s at %s GeV increments\n" % (len(masses), minstr, maxstr, str(round(args.deltaMT,precision)))

ttOnlySysts = ['toppt','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

#signals = ['tt','tW']
signals = ['tttW']

if len(signals) == 1 and signals[0] == "tt":
    print "TTbar only signal. Removing DS systematic"
    theorySysts.remove("DS")
backgrounds = ['DY','TTV','Diboson','ST_bkgd'] #,'WJets']

if args.obs == "diff":
    cats = [(0,"rec_ptll_M0_E0"), \
            (1,"rec_ptll_M0_E1"), \
            (2,"rec_ptll_M0_E2"), \
            (3,"rec_ptll_M1_E0"), \
            (4,"rec_ptll_M1_E1"), \
            (5,"rec_ptll_M1_E2"), \
            (6,"rec_ptll_M2_E0"), \
            (7,"rec_ptll_M2_E1"), \
            (8,"rec_ptll_M2_E2"), \
    ]
else:
    cats = [(0,"%s_%s" %(args.reco,args.obs)), \
        ]

if addBinStats and args.obs != "diff" and 'tttW' in signals:
    f = ROOT.TFile.Open(args.inF, "read")
    nbins = f.Get("%s_%s/tttW%s" % (args.reco, args.obs,masses[0])).GetNbinsX()
    print "Using binwise statistical uncertainty nps in %d bins" % nbins
    for _b in xrange(1, nbins+1):
        systematics["bin%d"%_b] = "shape"
    f.Close()

cb.AddObservations(['*'],['tt'],['13TeV'],['mtop'],cats)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['DY'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['TTV'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['Diboson'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['ST_bkgd'],  cats,False)
#cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['WJets'],   cats,False)


for sig in signals:
    cb.AddProcesses(masses,  ['tt'],['13TeV'],['mtop'],[sig],cats,True)


# Add systematics
if len(systematics) > 0:
    for syst,systype in systematics.items():
        if syst in ttOnlySysts and 'tt' in signals: 
            cb.cp().process(['tt']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))
        elif syst in tWOnlySysts and 'tW' in signals:
            cb.cp().process(['tW']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))
        elif syst == "BkgNorm":
            cb.cp().process(['DY']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.76923,1.30)) )
            cb.cp().process(['tW']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.94787,1.055)) )
        elif syst == "Lumi":
            cb.cp().process(signals+backgrounds).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.97561,1.025)) )
        else:
            cb.cp().signals().AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))


# Extract templates from input root file
cb.cp().backgrounds().ExtractShapes(args.inF,"$BIN/$PROCESS","$BIN/$PROCESS")
cb.cp().signals().ExtractShapes(args.inF,"$BIN/$PROCESS$MASS","$BIN/$PROCESS$MASS_$SYSTEMATIC")
#if len(systematics) > 0:
#    cb.cp().signals().ExtractShapes(args.inF,"$BIN/$PROCESS$MASS","$BIN/$PROCESS$MASS_$SYSTEMATIC")
#else:
#    print "Extracting shapes without systematics"
#    cb.cp().signals().ExtractShapes(args.inF,"$BIN/$PROCESS$MASS","$BIN/$PROCESS$MASS")


# Bin-by-bin uncertainties
#bbb = ch.BinByBinFactory()
#bbb.SetAddThreshold(-0.1).SetMergeThreshold(0.01).SetFixNorm(False).SetVerbosity(2)
##bbb.SetAddThreshold(0.1).SetMergeThreshold(0.5).SetFixNorm(False)
#bbb.MergeBinErrors(cb.cp().backgrounds())
#bbb.AddBinByBin(cb.cp().backgrounds(), cb)
#bbb.MergeBinErrors(cb.cp().signals())
#bbb.AddBinByBin(cb.cp().signals(), cb)
#ch.SetStandardBinNames(cb)

# Create workspace
w = ROOT.RooWorkspace('morph', 'morph')

mT = ROOT.RooRealVar("MT","",172.5*decimalScaling,165*decimalScaling,180*decimalScaling)
mT.setConstant(True)


debug = ROOT.TFile('debug.root', 'RECREATE')

# Create morphing splines
print "BuildRooMorphing"
for b in cb.cp().channel(['mtop']).bin_set():   # rec_ptll, rec_ptpos
    for sig in signals: # tt, tW
        morphing.BuildRooMorphing(w, cb, b, sig, mT,force_template_limit=False,file=debug)
        #morphing.BuildRooMorphing(w, cb, b, sig, mT,force_template_limit=False,file=debug, interp_mode="CSPLINE")

cb.AddWorkspace(w,True)

print "Extracting pdfs"
# Extract morphed templates
cb.cp().signals().ExtractPdfs(cb, 'morph', '$BIN_$PROCESS_morph', '')

print "Done extracting pdfs"
# Set up groups
if len(systematics) > 0:
    cb.SetGroup('all', theorySysts + experimentalSysts)
    cb.SetGroup('theory', theorySysts)
    cb.SetGroup('exp', experimentalSysts)

# Write to datacard
cb.PrintAll()
#cb.WriteDatacard(args.outF, "outputfileMorph.root")
f = ROOT.TFile.Open(outRootName, "recreate")
cb.WriteDatacard(outCardName, f)
f.Close()
#cb.WriteDatacard(outCardName, outRootName)

print "Datacard saved to %s" % outCardName

#!/usr/bin/env python
import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT
from argparse import ArgumentParser
import sys


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

#    "amcanlo":"shape", \
#    "madgraph":"shape", \
#    "herwigpp":"shape", \
    }

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input template root file")
parser.add_argument("-o", "--outF", default="", help="append name to output datacard and rootfile")
parser.add_argument("--splineInterp", default="CSPLINE", help="CH roofit morphing interpolation mode")
parser.add_argument("--systs", default="", nargs="+", choices=systematics.keys(), help="ONLY plot these systematics") 
parser.add_argument("--obs", default="ptll", help="observable")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reco level")
args = parser.parse_args()

outCardName = "%scardMorph.txt" % ("%s_" % args.outF if args.outF != "" else "")
outRootName = "%soutputfileMorph.root" % ("%s_" % args.outF if args.outF != "" else "")

cb = ch.CombineHarvester()
cb.SetVerbosity(3)

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
deltaM = 1
masses = [str(1665 + i * deltaM) for i in xrange(121)]
#masses = [str(1695 + i * deltaM) for i in xrange(61)]

ttOnlySysts = ['toppt','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

signals = ['tt','tW']
backgrounds = ['DY','TTV','Diboson','ST_bkgd','WJets']
cats = [(0,"%s_%s" %(args.reco,args.obs)), \
#        (1,'rec_ptpos'), \
#        (2,'rec_ptp_ptm'), \
#        (3,'rec_Ep_Em'), \
#        (4,'rec_Epos'), \
#        (5,'rec_Mll'), \
        ]
cb.AddObservations(['*'],['tt'],['13TeV'],['mtop'],cats)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['DY'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['TTV'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['Diboson'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['ST_bkgd'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['WJets'],   cats,False)
cb.AddProcesses(masses,  ['tt'],['13TeV'],['mtop'],['tt'],cats,True)
cb.AddProcesses(masses,  ['tt'],['13TeV'],['mtop'],['tW'],cats,True)


# Add systematics
for syst,systype in systematics.items():
    if syst in ttOnlySysts: 
        cb.cp().process(['tt']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))
    elif syst in tWOnlySysts:
        cb.cp().process(['tW']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))
    elif syst == "BkgNorm":
        cb.cp().process(['DY']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.70,1.30)) )
        cb.cp().process(['tW']).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.945,1.055)) )
    elif syst == "Lumi":
        cb.cp().process(signals+backgrounds).AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], (0.975,1.025)) )
    else:
        cb.cp().signals().AddSyst(cb,syst,systype,ch.SystMap('era')(['13TeV'], 1.))


# Extract templates from input root file
cb.cp().backgrounds().ExtractShapes(args.inF,"$BIN/$PROCESS","$BIN/$PROCESS")
cb.cp().signals().ExtractShapes(args.inF,"$BIN/$PROCESS$MASS","$BIN/$PROCESS$MASS_$SYSTEMATIC")


# Bin-by-bin uncertainties
#bbb = ch.BinByBinFactory()
#bbb.SetAddThreshold(0.1).SetMergeThreshold(0.5).SetFixNorm(True)
#bbb.MergeBinErrors(cb.cp().backgrounds())
#bbb.AddBinByBin(cb.cp().backgrounds(), cb)
#bbb.MergeBinErrors(cb.cp().signals())
#bbb.AddBinByBin(cb.cp().signals(), cb)


# Create workspace
w = ROOT.RooWorkspace('morph', 'morph')

mT = ROOT.RooRealVar("MT","",1730,1650,1800)
mT.setConstant(True)


debug = ROOT.TFile('debug.root', 'RECREATE')

# Create morphing splines
for b in cb.cp().channel(['mtop']).bin_set():   # rec_ptll, rec_ptpos
    for sig in signals: # tt, tW
        morphing.BuildRooMorphing(w, cb, b, sig, mT,force_template_limit=False,file=debug)
        #morphing.BuildRooMorphing(w, cb, b, sig, mT,force_template_limit=False,file=debug, interp_mode="CSPLINE")

cb.AddWorkspace(w,True)

# Extract morphed templates
cb.cp().signals().ExtractPdfs(cb, 'morph', '$BIN_$PROCESS_morph', '')


# Write to datacard
cb.PrintAll()
#cb.WriteDatacard(args.outF, "outputfileMorph.root")
cb.WriteDatacard(outCardName, outRootName)

print "Datacard saved to %s" % outCardName

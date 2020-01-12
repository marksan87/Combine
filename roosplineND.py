#!/usr/bin/env python
from ROOT import *
import os
import sys
from argparse import ArgumentParser
from array import array

# Load dynamic library with RooSplineND class
gSystem.Load("libHiggsAnalysisCombinedLimit.so")

observables = ["ptll","Mll","Epos","Eneg","ptpos","ptneg","Ep_Em","ptp_ptm"]
obsTitle = {\
        "ptll":"p_{T}(ll)",
        "ptpos":"p_{T}(l^{+})",
        "ptneg":"p_{T}(l^{-})",
        "Epos":"E(l^{+})",
        "Eneg":"E(l^{-})",
        "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})",
        "Ep_Em":"E(l^{+}) + E(l^{-})",
        "Mll":"M(ll)"
        }
parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="dataobs_ptll_bin10_mtTemplatesForCH.root", help="input root file")
parser.add_argument("--obs", default="ptll", choices=observables, help="observable")
parser.add_argument("-e", "--eps", type=float, default=1., help="tunable parameter for radial basis function")
parser.add_argument("--reco", default="rec", choices=["rec","gen"])
parser.add_argument("-o", "--outDir", default="RooSplineND_output", help="output directory")
args = parser.parse_args()

inF = args.inF
if inF == "":
    print "Missing input file!"
    sys.exit()
outDir = args.outDir
if outDir[-1] == "/": outDir = outdir[:-1]

obs = args.obs
reco = args.reco
outDir = args.outDir
recoObs = "%s_%s" % (reco,obs)


h = {}
masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
#masses = [1695, 1715, 1725, 1735, 1755]
#masses = [1715, 1725, 1735]

f = TFile.Open(inF, "read")
for m in masses:
    h[m] = f.Get("%s/ttactual%d" % (recoObs,m))
    h[m].SetDirectory(0)
f.Close()


tree = TTree("templates", "input template histogram data: mt vs obserable")
nEntries = len(masses) * h[1725].GetNbinsX()

binCenter = array('f', [0.])
binContent = array('f', [0.])
binError = array('f', [0.])     # TODO: Use errors!
mt = array('f', [0.])

tree.Branch("%s" % recoObs, binCenter, "%s/F" % recoObs)
tree.Branch("Events", binContent, "Events/F")
tree.Branch("Error", binError, "Error/F")
tree.Branch("mt", mt, "mt/F")

actualG = TGraph2DErrors()
actualG.SetName("actual")
actualG.SetTitle("%s %s  actual ; %s %s [GeV] ; m_{t} [GeV] ; Events" % (reco,obsTitle[obs],reco,obsTitle[obs]))

# Fill TTree
for m in masses:
    for b in range(1,h[m].GetNbinsX()+1):
        mt[0] = m/10.
        binCenter[0] = h[m].GetBinCenter(b)
        binContent[0] = h[m].GetBinContent(b)
        binError[0] = h[m].GetBinContent(b)

        point = actualG.GetN()
        actualG.SetPoint(point, binCenter[0], mt[0], binContent[0]) 
        actualG.SetPointError(point, 0., 0., binError[0]) 

        tree.Fill()

tree.Print()

obsVal = RooRealVar(recoObs, recoObs, 1.0, -10., 210.)
#mtVal  = RooRealVar("mt","mt", 0.05, 166.5, 178.5)
mtVal  = RooRealVar("mt","mt", 1, 165, 180)
#mtVal  = RooRealVar("mt","mt", 0.1, 165, 180)
argList = RooArgList(obsVal,mtVal)

# Interpolation width, tunable parameter
#eps = 2.5    
eps = 1.5
spline = RooSplineND("spline", "spline", argList, tree, "Events", eps, True, "", "TPS")

morphG = TGraph2D()
morphG.SetName("morph")
morphG.SetTitle("%s %s  morphed ; %s %s [GeV] ; m_{t} [GeV] ; Events" % (reco,obsTitle[obs],reco,obsTitle[obs]))

morphedMasses = [m/10. for m in range(min(masses), max(masses)+1)]
#morphedObsVals = [(o+0.5) for o in range(30,171)]
morphedObsVals = [(o+0.5) for o in range(0,201)]

for mt in morphedMasses: 
    for o in morphedObsVals:
        
        obsVal.setVal(o)
        mtVal.setVal(mt)

        point = morphG.GetN()
        morphG.SetPoint(point, o, mt, spline.getVal())

#actualG.Draw("err")
#morphG.Draw("Pcol same")
#morphG.Draw("tri2")


os.system("mkdir -p %s" % outDir)
f = TFile.Open("%s/%s.root" % (outDir, os.path.basename(outDir)), "recreate")

actualG.Write()
morphG.Write()



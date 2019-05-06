#!/usr/bin/env python
import sys
from ROOT import TH2D, TCanvas, TFile 
from argparse import ArgumentParser
from array import array

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}
signalTitle = {"tt":"t#bar{t}", "tW":"tW"}

parser = ArgumentParser()
parser.add_argument("-i", dest="inF", default="mtTemplatesForCH.root", help="Input debug root file from CombineHarvester with morphed templates")
#parser.add_argument("--sig", default="tt", choices=["tt","tW"], help="signal sample")
parser.add_argument("-s", "--sys", default="", help="systematic")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reconstruction level")
parser.add_argument("--obs", dest="obs", default="ptll", help="Kinematic observable")
parser.add_argument("-o", dest="outF", default="morph2D.root", help="Output root file")
args = parser.parse_args()

if args.obs not in obsTitle.keys():
    print "Invalid observable:", obs
    print "Choose from:"
    print obsTitle.keys()
    sys.exit()

if args.outF.find(".root") < 0:
    args.outF += ".root"

f = TFile.Open(args.inF)
deltaM = 1
masses = range(1665, 1786, deltaM)
#masses = range(1665, 1795, 10)
#signal = ["tt","tW"]
massBins = [m/10. for m in masses]

signal = ["tt","tW"]
#systematics = ["nominal", "pileupUp", "pileupDown", "Q2Up", "Q2Down", "PdfUp", "PdfDown"]
systematics = ["nominal"]
recoObs = "%s_%s" % (args.reco, args.obs)
h = {}
for s in signal:
    h[s] = {}
    for syst in systematics:
        morph = {}
        for m in masses:
            #morph[m] = f.Get("rec_%s/%s%d%s" % (args.obs,s,m,"" if syst=="nominal" else "_" + syst))
            morph[m] = f.Get("%s/%s%d%s" % (recoObs,s,m,"" if syst=="nominal" else "_" + syst))
            morph[m].SetDirectory(0)

        obsBins = morph[masses[0]].GetXaxis().GetXbins()
        obsBins = [p for p in obsBins]
        if len(obsBins) < 1:
            # Uniform binning
            obsBins = []
            for _b in xrange(1, morph[masses[0]].GetNbinsX()+2):
                obsBins.append(morph[masses[0]].GetBinLowEdge(_b))

        h[s][syst] = TH2D("%s%s_morph_points" % (s,"" if syst=="nominal" else "_" + syst), "%s Morphed %s%s Templates" % (signalTitle[s], obsTitle[args.obs],"" if syst == "nominal" else "_" + syst), len(massBins)-1, array('d',massBins), morph[masses[0]].GetNbinsX(), array('d', obsBins) )
        h[s][syst].GetXaxis().SetTitle("m_{t} [GeV]")
        h[s][syst].GetYaxis().SetTitle("%s [GeV]" % obsTitle[args.obs])
        h[s][syst].GetZaxis().SetTitle("Events")

        h[s][syst].GetXaxis().SetTitleOffset(1.6)
        h[s][syst].GetYaxis().SetTitleOffset(1.8)
        h[s][syst].GetZaxis().SetTitleOffset(1.5)

        for m in masses:
            for b in xrange(1, morph[m].GetNbinsX() + 1):
                h[s][syst].Fill(m / 10.0, morph[m].GetBinCenter(b), morph[m].GetBinContent(b))


outF = TFile.Open(args.outF, "recreate")
for s in signal:
    for syst in systematics:
        h[s][syst].Write()
#for m in masses:
#    morph[m].Write()

outF.Close()

print "Results saved in %s" % args.outF

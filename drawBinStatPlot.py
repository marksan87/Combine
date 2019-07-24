#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser

obsTitle = {"ptll":"p_{T}(ll)",
            "ptpos":"p_{T}(l^{+})",
            "Epos":"E(l^{+})",
            "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})",
            "Ep_Em":"E(l^{+}) + E(l^{-})",
            "Mll":"M(ll)",
            "leadLepPt":"Leading Lepton p_{T}",
            "leadJetPt":"Leading Jet p_{T}",
            "ptll_M0_E0":"p_{T}(ll)  M(ll) < 80  E(l^{+}) < 60",
            "ptll_M0_E1":"p_{T}(ll)  M(ll) < 80  60 <= E(l^{+}) < 100",
            "ptll_M0_E2":"p_{T}(ll)  M(ll) < 80  100 <= E(l^{+}) < 600",
            "ptll_M1_E0":"p_{T}(ll)  80 <= M(ll) < 130  E(l^{+}) < 60",
            "ptll_M1_E1":"p_{T}(ll)  80 <= M(ll) < 130  60 <= E(l^{+}) < 100",
            "ptll_M1_E2":"p_{T}(ll)  80 <= M(ll) < 130  100 <= E(l^{+}) < 600",
            "ptll_M2_E0":"p_{T}(ll)  130 <= M(ll) < 600  E(l^{+}) < 60",
            "ptll_M2_E1":"p_{T}(ll)  130 <= M(ll) < 600  60 <= E(l^{+}) < 100",
            "ptll_M2_E2":"p_{T}(ll)  130 <= M(ll) < 600  100 <= E(l^{+}) < 600",
            }

gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument("-i", "--inF", help="input template file with bin stats")
parser.add_argument("-b", "--bin", type=int, default=3, help="bin np")
parser.add_argument("--obs", default="ptll")
parser.add_argument("--reco", default="rec")
parser.add_argument("-m", "--mass", type=int, default=1725, help="10*mt")
parser.add_argument("-o", "--outF", default="", help="output plot file")
args = parser.parse_args()

recoObs = "%s_%s" % (args.reco, args.obs)
if args.outF == "":
    args.outF = "%s_MC_stats_bin_%d.png" % (recoObs,args.bin)

print "Plotting %s  mt = %.1f  bin %d" % (recoObs, args.mass/10., args.bin)

f = TFile.Open(args.inF, "read")

nom = f.Get("%s/tttW%d" % (recoObs, args.mass))
nom.SetDirectory(0)

up = f.Get("%s/tttW%d_bin%dUp" % (recoObs, args.mass, args.bin))
up.SetDirectory(0)
dn = f.Get("%s/tttW%d_bin%dDown" % (recoObs, args.mass, args.bin))
dn.SetDirectory(0)

nom.SetLineWidth(2)
up.SetLineWidth(2)
dn.SetLineWidth(2)

nom.SetLineColor(kBlack)
up.SetLineColor(kRed)
dn.SetLineColor(kBlue)

#nom.Scale(1., "width")
#up.Scale(1., "width")
#dn.Scale(1., "width")

l = TLegend(0.12, 0.7, 0.38, 0.88)
l.SetBorderSize(0)

l.AddEntry(up, "up")
l.AddEntry(nom, "nom")
l.AddEntry(dn, "down")

c = TCanvas("c","c",1200,800)

up.SetTitle("%s %s  m_{t} = %.1f GeV  Bin %d" % (args.reco, obsTitle[args.obs], args.mass/10., args.bin))
up.GetXaxis().SetTitle("%s %s [GeV]" % (args.reco, obsTitle[args.obs]))
up.GetXaxis().SetTitleOffset(1.2)
up.GetYaxis().SetTitle("Entries")
up.GetYaxis().SetTitleOffset(1.4)
up.Draw("hist")
dn.Draw("hist same")
nom.Draw("hist same")
l.Draw("same")

c.SaveAs(args.outF)

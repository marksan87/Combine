#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import sys
gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="diff_mtTemplatesForCH.root")
parser.add_argument("-o", "--outF", default="dndist.png")
args = parser.parse_args()

_diffHists = ["rec_ptll_M0_E0", "rec_ptll_M0_E1", "rec_ptll_M0_E2", "rec_ptll_M1_E0", "rec_ptll_M1_E1", "rec_ptll_M1_E2", "rec_ptll_M2_E0", "rec_ptll_M2_E1", "rec_ptll_M2_E2"]
#_diffHists = ["rec_ptll_M0_E0", "rec_ptll_M0_E1", "rec_ptll_M0_E2"]

_diffText = {
    "rec_ptll_M0_E0":"{0 < Mll < 80}{0 < Epos < 60}",
    "rec_ptll_M0_E1":"{0 < Mll < 80}{60 < Epos < 100}",
    "rec_ptll_M0_E2":"{0 < Mll < 80}{100 < Epos < 600}",
    "rec_ptll_M1_E0":"{80 < Mll < 130}{0 < Epos < 60}",
    "rec_ptll_M1_E1":"{80 < Mll < 130}{60 < Epos < 100}",
    "rec_ptll_M1_E2":"{80 < Mll < 130}{100 < Epos < 600}",
    "rec_ptll_M2_E0":"{130 < Mll < 600}{0 < Epos < 60}",
    "rec_ptll_M2_E1":"{130 < Mll < 600}{60 < Epos < 100}",
    "rec_ptll_M2_E2":"{130 < Mll < 600}{100 < Epos < 600}",
}

masses = [1695,1725,1755]
colors = {1695:kBlue, 1725:kBlack, 1755:kRed}
h = {}
ratioUp = {}
ratioDown = {}
l = {}
nomRatioWithErr = {}
f = TFile.Open(args.inF, 'read')
print "Loading templates from %s..." % args.inF,
for d in _diffHists:
    h[d] = {}
    for m in masses:
        h[d][m] = f.Get("%s/tttW%d" % (d,m))
        h[d][m].SetDirectory(0)
        h[d][m].SetLineColor(colors[m])
        h[d][m].SetLineWidth(2)
        h[d][m].SetTitle("")

    nomNoErr = h[d][1725].Clone("nomNoErr_%s" % d)
    for b in xrange(1,nomNoErr.GetNbinsX()+1):
        nomNoErr.SetBinError(b,0.)
    ratioUp[d] = h[d][1755].Clone("ratioUp_%s" % d)
    ratioDown[d] = h[d][1695].Clone("ratioDown_%s" % d)
    ratioUp[d].SetDirectory(0)
    ratioUp[d].SetTitle("")
    ratioDown[d].SetDirectory(0)
    ratioDown[d].SetTitle("")
    #ratioUp[d].Divide(h[d][1725])
    #ratioDown[d].Divide(h[d][1725])
    ratioUp[d].Divide(nomNoErr)
    ratioDown[d].Divide(nomNoErr)
    
    nomRatioWithErr[d] = h[d][1725].Clone("nomErr_%s" % d)
    nomRatioWithErr[d].SetDirectory(0)
    nomRatioWithErr[d].Divide(nomNoErr)
    
    #l[d] = TLine(ratioUp[d].GetXaxis().GetBinLowEdge(1), 1., ratioUp[d].GetXaxis().GetBinUpEdge(ratioUp[d].GetNbinsX()), 1.)
    #l[d].SetLineWidth(2)

f.Close()
print "done"
sys.stdout.flush()
#h["rec_ptll_M0_E0"][1725].Draw("hist")


c = TCanvas("c","c",3600,800)

c.Divide(len(_diffHists),1)

l = TLine(ratioUp[_diffHists[0]].GetXaxis().GetBinLowEdge(1), 1., ratioUp[_diffHists[0]].GetXaxis().GetBinUpEdge(ratioUp[_diffHists[0]].GetNbinsX()), 1.)
l.SetLineWidth(2)

t = TLatex()
#t.SetTextSize(0.1)
t.SetTextAlign(12)
for i,d in enumerate(_diffHists):
    c.cd(i+1)
    #c.Draw()
    pad1 = TPad("pad1_%d"%i, "pad1_%d"%i, 0., 0.25, 1., 1.)
    pad2 = TPad("pad2_%d"%i, "pad2_%d"%i, 0., 0., 1., 0.25)
    pad1.Draw()
    pad2.Draw()
    pad1.cd()
    for j,m in enumerate(masses):
        if j == 0:
            h[d][m].Draw("hist")
            h[d][m].GetXaxis().SetTitle("p_{T}(ll) [GeV]")
            h[d][m].GetXaxis().SetTitleSize(0.06)
            h[d][m].GetXaxis().SetTitleOffset(0.7)
            #h[d][m].GetYaxis().SetTitle("Entries")
            #h[d][m].GetYaxis().SetTitleOffset(1.6)
            h[d][m].GetYaxis().SetTitleSize(0.06)
            h[d][m].GetYaxis().SetTitleOffset(1.0)
        else:
            h[d][m].Draw("hist same")

    #t.DrawLatexNDC(0.5, 0.5, "#splitline{0 < M < 10}{3 < E < 10}")
    t.DrawLatexNDC(0.5, 0.5, "#splitline%s" % _diffText[d])

    pad2.cd()
#    ratioUp = h[d][1755].Clone("%s_ratioUp" % d)
#    ratioUp.Divide(h[d][1725])
#    ratioDown = h[d][1695].Clone("%s_ratioDown" % d)
#    ratioDown.Divide(h[d][1725])
#
#    maxRatioY = 1.0
#    maxRatioY = max(maxRatioY, max(ratioUp.GetMaximum())))
#    minRatioY = 1.0
#    minRatioY = min(minRatioY, min(ratioUp.GetMinimum())))
    minRatioY = 0.88
    maxRatioY = 1.12
    
    ratioUp[d].Draw("hist e1")
    ratioUp[d].GetYaxis().SetRangeUser(minRatioY, maxRatioY)
    ratioUp[d].GetYaxis().SetTitle("ratio wrt nominal")
    ratioUp[d].GetYaxis().SetTitleOffset(0.5)
    ratioUp[d].GetYaxis().SetTitleSize(0.1)
    ratioDown[d].Draw("hist e1 same")
    #l.Draw("same")
    nomRatioWithErr[d].Draw("hist e1 same")

c.SaveAs(args.outF)

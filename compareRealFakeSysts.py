#!/usr/bin/env python
from ROOT import *
import os
import sys
from array import array
from pprint import pprint
from argparse import ArgumentParser

gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument("-r", "--real", default="bin10_cut220_mtTemplatesForCH.root", help="template file containing real systematics")
parser.add_argument("-f", "--fake", default="fakeJEC_bin10_cut220_mtTemplatesForCH.root", help="template file containing fake systematics")
parser.add_argument("-a", "--useActual", action="store_true", default=False, help="plot actual templates instead of morphed")
parser.add_argument("-s", "--syst", default="JEC", help="systematic to plot")
parser.add_argument("-o", "--outDir", default="fakeRealSyst", help="output directory")

args = parser.parse_args()
if args.outDir[-1] == "/": args.outDir = args.outDir[:-1]
os.system("mkdir -p %s" % args.outDir)

recoObs = "rec_ptll"
signal = "tt" + ("actual" if args.useActual else "")
gROOT.SetBatch(True)
masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]

c = TCanvas("c", "c", 800, 1200)
pad1 = TPad("pad1", "pad1", 0., 0.5, 1., 1.)
pad2 = TPad("pad2", "pad2", 0., 0.25, 1., 0.5)
pad3 = TPad("pad3", "pad3", 0., 0., 1., 0.25)

canvasCompare = TCanvas("c2", "c2", 800,1200)
compPad1 = TPad("compPad1", "compPad1", 0., 0.5, 1., 1.)
compPad2 = TPad("compPad2", "compPad2", 0., 0., 1., 0.5)

realF = TFile.Open(args.real, "read")
fakeF = TFile.Open(args.fake, "read")

for m in masses:
    print "Now on mass %d" % m
    nom = realF.Get("%s/ttactual%d" % (recoObs,m)).Clone()
    nom.SetDirectory(0)
    nom.SetLineWidth(2)

    realUp = realF.Get("%s/%s%d_%sUp" % (recoObs, signal, m, args.syst)).Clone()
    realDn = realF.Get("%s/%s%d_%sDown" % (recoObs, signal, m, args.syst)).Clone()

    realUp.SetDirectory(0)
    realUp.SetLineColor(kRed)
    realUp.SetFillColor(kWhite)
    realUp.SetLineWidth(2)
    realDn.SetDirectory(0)
    realDn.SetLineColor(kBlue)
    realDn.SetFillColor(kWhite)
    realDn.SetLineWidth(2)

    resRealUp = realUp.Clone()
    resRealUp.SetDirectory(0)
    resRealUp.Add(nom, -1)
    ratioRealUp = realUp.Clone()
    ratioRealUp.SetDirectory(0)
    ratioRealUp.Divide(nom)

    resRealDn = realDn.Clone()
    resRealDn.SetDirectory(0)
    resRealDn.Add(nom, -1)
    ratioRealDn = realDn.Clone()
    ratioRealDn.SetDirectory(0)
    ratioRealDn.Divide(nom)


    maxRealRes = -100
    maxRealRatio = -100
    minRealRes = 100
    minRealRatio = 100

    maxRealRes = max(maxRealRes, max(resRealUp.GetMaximum(), resRealDn.GetMaximum()))
    minRealRes = min(minRealRes, min(resRealUp.GetMinimum(), resRealDn.GetMinimum()))
    maxRealRatio = max(maxRealRatio, max(ratioRealUp.GetMaximum(), ratioRealDn.GetMaximum()))
    minRealRatio = min(minRealRatio, min(ratioRealUp.GetMinimum(), ratioRealDn.GetMinimum()))

    paddingRealRes = 0.05 * abs(maxRealRes - minRealRes)
    paddingRealRatio = 0.05 * abs(maxRealRatio - minRealRatio)

    fakeUp = fakeF.Get("%s/%s%d_%sUp" % (recoObs, signal, m, args.syst)).Clone()
    fakeDn = fakeF.Get("%s/%s%d_%sDown" % (recoObs, signal, m, args.syst)).Clone()

    fakeUp.SetDirectory(0)
    fakeUp.SetLineColor(kRed)
    fakeUp.SetFillColor(kWhite)
    fakeUp.SetLineWidth(2)
    fakeDn.SetDirectory(0)
    fakeDn.SetLineColor(kBlue)
    fakeDn.SetFillColor(kWhite)
    fakeDn.SetLineWidth(2)

    resFakeUp = fakeUp.Clone()
    resFakeUp.SetDirectory(0)
    resFakeUp.Add(nom, -1)
    ratioFakeUp = fakeUp.Clone()
    ratioFakeUp.SetDirectory(0)
    ratioFakeUp.Divide(nom)

    resFakeDn = fakeDn.Clone()
    resFakeDn.SetDirectory(0)
    resFakeDn.Add(nom, -1)
    ratioFakeDn = fakeDn.Clone()
    ratioFakeDn.SetDirectory(0)
    ratioFakeDn.Divide(nom)

    resLine = TLine(nom.GetXaxis().GetBinLowEdge(1), 0., nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 0.)
    ratioLine = TLine(nom.GetXaxis().GetBinLowEdge(1), 1., nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 1.)

    resLine.SetLineWidth(2)
    ratioLine.SetLineWidth(2)

    maxFakeRes = -100
    maxFakeRatio = -100
    minFakeRes = 100
    minFakeRatio = 100

    maxFakeRes = max(maxFakeRes, max(resFakeUp.GetMaximum(), resFakeDn.GetMaximum()))
    minFakeRes = min(minFakeRes, min(resFakeUp.GetMinimum(), resFakeDn.GetMinimum()))
    maxFakeRatio = max(maxFakeRatio, max(ratioFakeUp.GetMaximum(), ratioFakeDn.GetMaximum()))
    minFakeRatio = min(minFakeRatio, min(ratioFakeUp.GetMinimum(), ratioFakeDn.GetMinimum()))

    paddingFakeRes = 0.05 * abs(maxFakeRes - minFakeRes)
    paddingFakeRatio = 0.05 * abs(maxFakeRatio - minFakeRatio)


    maxReal = max(realUp.GetMaximum(), max(nom.GetMaximum(), realDn.GetMaximum()))
    maxFake = max(fakeUp.GetMaximum(), max(nom.GetMaximum(), fakeDn.GetMaximum()))

    l = TLegend(0.7, 0.7, 0.88, 0.88)
    l.SetBorderSize(0)
    l.AddEntry(realUp, "%s Up" % args.syst)
    l.AddEntry(nom, "nominal")
    l.AddEntry(realDn, "%s Dn" % args.syst)

    c.cd()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### Real syst ###
    pad1.cd()
    realUp.SetTitle("Real %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    realUp.GetYaxis().SetRangeUser(0, 1.05*maxReal) 
    realUp.Draw("hist")
    nom.Draw("hist same")
    realDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resRealUp.GetYaxis().SetRangeUser(minRealRes - paddingRealRes, maxRealRes + paddingRealRes)
    resRealUp.SetTitle("Residual")
    resRealUp.Draw("hist")
    resRealDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioRealUp.GetYaxis().SetRangeUser(minRealRatio - paddingRealRatio, maxRealRatio + paddingRealRatio)
    ratioRealUp.SetTitle("Ratio")
    ratioRealUp.Draw("hist")
    ratioRealDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/real%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    c.cd()
    c.ResetDrawn()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### Fake syst ###
    pad1.cd()
    fakeUp.SetTitle("Fake %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    fakeUp.GetYaxis().SetRangeUser(0, 1.05*maxFake)
    fakeUp.Draw("hist")
    nom.Draw("hist same")
    fakeDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resFakeUp.GetYaxis().SetRangeUser(minFakeRes - paddingFakeRes, maxFakeRes + paddingFakeRes)
    resFakeUp.SetTitle("Residual")
    resFakeUp.Draw("hist")
    resFakeDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioFakeUp.GetYaxis().SetRangeUser(minFakeRatio - paddingFakeRatio, maxFakeRatio + paddingFakeRatio)
    ratioFakeUp.SetTitle("Ratio")
    ratioFakeUp.Draw("hist")
    ratioFakeDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/fake%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    c.cd()
    c.ResetDrawn()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### Comparison ###
    pad1.cd()

    fakeUp.SetLineColor(kOrange)
    fakeDn.SetLineColor(kCyan)

    resFakeUp.SetLineColor(kOrange)
    ratioFakeUp.SetLineColor(kOrange)
    resFakeDn.SetLineColor(kCyan)
    ratioFakeDn.SetLineColor(kCyan)


    l = TLegend(0.7, 0.6, 0.88, 0.88)
    l.SetBorderSize(0)
    l.AddEntry(realUp, "Real %s Up" % args.syst)
    l.AddEntry(fakeUp, "Fake %s Up" % args.syst) 
    l.AddEntry(nom, "nominal")
    l.AddEntry(realDn, "Real %s Dn" % args.syst)
    l.AddEntry(fakeDn, "Fake %s Dn" % args.syst)

    realUp.SetTitle("Real and Fake %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    realUp.Draw("hist")
    fakeUp.Draw("hist same")
    nom.Draw("hist same")
    realDn.Draw("hist same")
    fakeDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resRealUp.GetYaxis().SetRangeUser(min(minRealRes, minFakeRes) - max(paddingRealRes, paddingFakeRes), max(maxRealRes, maxFakeRes) + max(paddingRealRes, paddingFakeRes))
    resRealUp.Draw("hist")
    resFakeUp.Draw("hist same")
    resRealDn.Draw("hist same")
    resFakeDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioRealUp.GetYaxis().SetRangeUser(min(minRealRatio, minFakeRatio) - max(paddingRealRatio, paddingFakeRatio), max(maxRealRatio, maxFakeRatio) + max(paddingRealRatio, paddingFakeRatio))
    ratioRealUp.Draw("hist")
    ratioFakeUp.Draw("hist same")
    ratioRealDn.Draw("hist same")
    ratioFakeDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/comparison_real_fake_%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))


    canvasCompare.cd()
    canvasCompare.ResetDrawn()
    canvasCompare.Draw()
    compPad1.Draw()
    compPad2.Draw()

    compPad1.cd()
    ratioFakeRealUp = fakeUp.Clone("ratioFakeRealUp")
    ratioFakeRealUp.SetDirectory(0)
    ratioFakeRealUp.Divide(realUp)
    ratioFakeRealUp.GetYaxis().SetRangeUser(ratioFakeRealUp.GetBinContent(ratioFakeRealUp.GetMinimumBin())*0.95, ratioFakeRealUp.GetBinContent(ratioFakeRealUp.GetMaximumBin())*1.05)
    ratioFakeRealUp.SetTitle("Fake / Real %s Up  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    ratioFakeRealUp.SetLineColor(kRed)

    ratioFakeRealDn = fakeDn.Clone("ratioFakeRealDn")
    ratioFakeRealDn.SetDirectory(0)
    ratioFakeRealDn.Divide(realDn)
    ratioFakeRealDn.SetTitle("Fake / Real %s Down  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    ratioFakeRealDn.SetLineColor(kBlue)

    compPad1.cd()
    ratioFakeRealUp.Draw("hist")
    ratioLine.Draw("same")

    compPad2.cd()
    ratioFakeRealDn.Draw("hist")
    ratioLine.Draw("same")
    canvasCompare.SaveAs("%s/ratio_fake_real_%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    
realF.Close()
fakeF.Close()

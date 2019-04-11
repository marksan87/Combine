#!/usr/bin/env python
from ROOT import *
import os
import sys
from array import array
from pprint import pprint
from argparse import ArgumentParser

gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument("--old", default="newMorphRates_bin10_cut220_mtTemplatesForCH.root", help="old template file")
parser.add_argument("--new", default="secondtry_morphRates_bin10_cut220_mtTemplatesForCH.root", help="new template file")
parser.add_argument("-a", "--useActual", action="store_true", default=False, help="plot actual templates instead of morphed")
parser.add_argument("-s", "--syst", default="EleScale", help="systematic to plot")
parser.add_argument("-o", "--outDir", default="compOldNewHists", help="output directory")

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

oldF = TFile.Open(args.old, "read")
newF = TFile.Open(args.new, "read")

for m in masses:
    print "Now on mass %d" % m
    nom = oldF.Get("%s/ttactual%d" % (recoObs,m)).Clone()
    nom.SetDirectory(0)
    nom.SetLineWidth(2)

    oldUp = oldF.Get("%s/%s%d_%sUp" % (recoObs, signal, m, args.syst)).Clone()
    oldDn = oldF.Get("%s/%s%d_%sDown" % (recoObs, signal, m, args.syst)).Clone()

    oldUp.SetDirectory(0)
    oldUp.SetLineColor(kRed)
    oldUp.SetFillColor(kWhite)
    oldUp.SetLineWidth(2)
    oldDn.SetDirectory(0)
    oldDn.SetLineColor(kBlue)
    oldDn.SetFillColor(kWhite)
    oldDn.SetLineWidth(2)

    resOldUp = oldUp.Clone()
    resOldUp.SetDirectory(0)
    resOldUp.Add(nom, -1)
    ratioOldUp = oldUp.Clone()
    ratioOldUp.SetDirectory(0)
    ratioOldUp.Divide(nom)

    resOldDn = oldDn.Clone()
    resOldDn.SetDirectory(0)
    resOldDn.Add(nom, -1)
    ratioOldDn = oldDn.Clone()
    ratioOldDn.SetDirectory(0)
    ratioOldDn.Divide(nom)


    maxOldRes = -100
    maxOldRatio = -100
    minOldRes = 100
    minOldRatio = 100

    maxOldRes = max(maxOldRes, max(resOldUp.GetMaximum(), resOldDn.GetMaximum()))
    minOldRes = min(minOldRes, min(resOldUp.GetMinimum(), resOldDn.GetMinimum()))
    maxOldRatio = max(maxOldRatio, max(ratioOldUp.GetMaximum(), ratioOldDn.GetMaximum()))
    minOldRatio = min(minOldRatio, min(ratioOldUp.GetMinimum(), ratioOldDn.GetMinimum()))

    paddingOldRes = 0.05 * abs(maxOldRes - minOldRes)
    paddingOldRatio = 0.05 * abs(maxOldRatio - minOldRatio)

    newUp = newF.Get("%s/%s%d_%sUp" % (recoObs, signal, m, args.syst)).Clone()
    newDn = newF.Get("%s/%s%d_%sDown" % (recoObs, signal, m, args.syst)).Clone()

    newUp.SetDirectory(0)
    newUp.SetLineColor(kRed)
    newUp.SetFillColor(kWhite)
    newUp.SetLineWidth(2)
    newDn.SetDirectory(0)
    newDn.SetLineColor(kBlue)
    newDn.SetFillColor(kWhite)
    newDn.SetLineWidth(2)

    resNewUp = newUp.Clone()
    resNewUp.SetDirectory(0)
    resNewUp.Add(nom, -1)
    ratioNewUp = newUp.Clone()
    ratioNewUp.SetDirectory(0)
    ratioNewUp.Divide(nom)

    resNewDn = newDn.Clone()
    resNewDn.SetDirectory(0)
    resNewDn.Add(nom, -1)
    ratioNewDn = newDn.Clone()
    ratioNewDn.SetDirectory(0)
    ratioNewDn.Divide(nom)

    resLine = TLine(nom.GetXaxis().GetBinLowEdge(1), 0., nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 0.)
    ratioLine = TLine(nom.GetXaxis().GetBinLowEdge(1), 1., nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 1.)

    resLine.SetLineWidth(2)
    ratioLine.SetLineWidth(2)

    maxNewRes = -100
    maxNewRatio = -100
    minNewRes = 100
    minNewRatio = 100

    maxNewRes = max(maxNewRes, max(resNewUp.GetMaximum(), resNewDn.GetMaximum()))
    minNewRes = min(minNewRes, min(resNewUp.GetMinimum(), resNewDn.GetMinimum()))
    maxNewRatio = max(maxNewRatio, max(ratioNewUp.GetMaximum(), ratioNewDn.GetMaximum()))
    minNewRatio = min(minNewRatio, min(ratioNewUp.GetMinimum(), ratioNewDn.GetMinimum()))

    paddingNewRes = 0.05 * abs(maxNewRes - minNewRes)
    paddingNewRatio = 0.05 * abs(maxNewRatio - minNewRatio)


    maxOld = max(oldUp.GetMaximum(), max(nom.GetMaximum(), oldDn.GetMaximum()))
    maxNew = max(newUp.GetMaximum(), max(nom.GetMaximum(), newDn.GetMaximum()))

    l = TLegend(0.7, 0.7, 0.88, 0.88)
    l.SetBorderSize(0)
    l.AddEntry(oldUp, "%s Up" % args.syst)
    l.AddEntry(nom, "nominal")
    l.AddEntry(oldDn, "%s Dn" % args.syst)

    c.cd()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### Old syst ###
    pad1.cd()
    oldUp.SetTitle("Old %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    oldUp.GetYaxis().SetRangeUser(0, 1.05*maxOld) 
    oldUp.Draw("hist")
    nom.Draw("hist same")
    oldDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resOldUp.GetYaxis().SetRangeUser(minOldRes - paddingOldRes, maxOldRes + paddingOldRes)
    resOldUp.SetTitle("Residual")
    resOldUp.Draw("hist")
    resOldDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioOldUp.GetYaxis().SetRangeUser(minOldRatio - paddingOldRatio, maxOldRatio + paddingOldRatio)
    ratioOldUp.SetTitle("Ratio")
    ratioOldUp.Draw("hist")
    ratioOldDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/old%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    c.cd()
    c.ResetDrawn()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### New syst ###
    pad1.cd()
    newUp.SetTitle("New %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    newUp.GetYaxis().SetRangeUser(0, 1.05*maxNew)
    newUp.Draw("hist")
    nom.Draw("hist same")
    newDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resNewUp.GetYaxis().SetRangeUser(minNewRes - paddingNewRes, maxNewRes + paddingNewRes)
    resNewUp.SetTitle("Residual")
    resNewUp.Draw("hist")
    resNewDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioNewUp.GetYaxis().SetRangeUser(minNewRatio - paddingNewRatio, maxNewRatio + paddingNewRatio)
    ratioNewUp.SetTitle("Ratio")
    ratioNewUp.Draw("hist")
    ratioNewDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/new%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    c.cd()
    c.ResetDrawn()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


### Comparison ###
    pad1.cd()

    newUp.SetLineColor(kOrange)
    newDn.SetLineColor(kCyan)

    resNewUp.SetLineColor(kOrange)
    ratioNewUp.SetLineColor(kOrange)
    resNewDn.SetLineColor(kCyan)
    ratioNewDn.SetLineColor(kCyan)


    l = TLegend(0.7, 0.6, 0.88, 0.88)
    l.SetBorderSize(0)
    l.AddEntry(oldUp, "Old %s Up" % args.syst)
    l.AddEntry(newUp, "New %s Up" % args.syst) 
    l.AddEntry(nom, "nominal")
    l.AddEntry(oldDn, "Old %s Dn" % args.syst)
    l.AddEntry(newDn, "New %s Dn" % args.syst)

    oldUp.SetTitle("Old and New %s  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    oldUp.Draw("hist")
    newUp.Draw("hist same")
    nom.Draw("hist same")
    oldDn.Draw("hist same")
    newDn.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resOldUp.GetYaxis().SetRangeUser(min(minOldRes, minNewRes) - max(paddingOldRes, paddingNewRes), max(maxOldRes, maxNewRes) + max(paddingOldRes, paddingNewRes))
    resOldUp.Draw("hist")
    resNewUp.Draw("hist same")
    resOldDn.Draw("hist same")
    resNewDn.Draw("hist same")
    resLine.Draw("same")

    pad3.cd()
    ratioOldUp.GetYaxis().SetRangeUser(min(minOldRatio, minNewRatio) - max(paddingOldRatio, paddingNewRatio), max(maxOldRatio, maxNewRatio) + max(paddingOldRatio, paddingNewRatio))
    ratioOldUp.Draw("hist")
    ratioNewUp.Draw("hist same")
    ratioOldDn.Draw("hist same")
    ratioNewDn.Draw("hist same")
    ratioLine.Draw("same")

    c.SaveAs("%s/comparison_old_new_%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))


    canvasCompare.cd()
    canvasCompare.ResetDrawn()
    canvasCompare.Draw()
    compPad1.Draw()
    compPad2.Draw()

    compPad1.cd()
    ratioNewOldUp = newUp.Clone("ratioNewOldUp")
    ratioNewOldUp.SetDirectory(0)
    ratioNewOldUp.Divide(oldUp)
    ratioNewOldUp.GetYaxis().SetRangeUser(ratioNewOldUp.GetBinContent(ratioNewOldUp.GetMinimumBin())*0.95, ratioNewOldUp.GetBinContent(ratioNewOldUp.GetMaximumBin())*1.05)
    ratioNewOldUp.SetTitle("New / Old %s Up  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    ratioNewOldUp.SetLineColor(kRed)

    ratioNewOldDn = newDn.Clone("ratioNewOldDn")
    ratioNewOldDn.SetDirectory(0)
    ratioNewOldDn.Divide(oldDn)
    ratioNewOldDn.SetTitle("New / Old %s Down  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    ratioNewOldDn.SetLineColor(kBlue)

    compPad1.cd()
    ratioNewOldUp.Draw("hist")
    ratioLine.Draw("same")

    compPad2.cd()
    ratioNewOldDn.Draw("hist")
    ratioLine.Draw("same")
    canvasCompare.SaveAs("%s/ratio_new_old_%s_mt%d_%s.png" % (args.outDir, args.syst, m, "actual" if args.useActual else "morphed"))

    
oldF.Close()
newF.Close()

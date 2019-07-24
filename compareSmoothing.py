#!/usr/bin/env python
from ROOT import *
import os
import sys
from array import array
from pprint import pprint
from argparse import ArgumentParser

gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument("--nosmooth", default="nosmooth_perbinMorph_ptll_q10_mtTemplatesForCH.root", help="old template file")
parser.add_argument("--smooth", default="smooth_perbinMorph_ptll_q10_mtTemplatesForCH.root", help="new template file")
parser.add_argument("-a", "--useActual", action="store_true", default=False, help="plot actual templates instead of morphed")
parser.add_argument("--noerrors", action="store_true", default=False, help="don't draw error bars")
parser.add_argument("-s", "--syst", default="", help="systematic to plot")
parser.add_argument("-o", "--outDir", default="compSmoothing", help="output directory")

args = parser.parse_args()
if args.outDir[-1] == "/": args.outDir = args.outDir[:-1]
os.system("mkdir -p %s" % args.outDir)

recoObs = "rec_ptll"
signal = "tt"

gROOT.SetBatch(True)
masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]

c = TCanvas("c", "c", 800, 1200)
pad1 = TPad("pad1", "pad1", 0., 0.5, 1., 1.)
pad2 = TPad("pad2", "pad2", 0., 0.25, 1., 0.5)
pad3 = TPad("pad3", "pad3", 0., 0., 1., 0.25)

canvasCompare = TCanvas("c2", "c2", 800,1200)
compPad1 = TPad("compPad1", "compPad1", 0., 0.5, 1., 1.)
compPad2 = TPad("compPad2", "compPad2", 0., 0., 1., 0.5)

nosmoothF = TFile.Open(args.nosmooth, "read")
smoothF = TFile.Open(args.smooth, "read")

for m in masses:
    print "Now on mass %d" % m
    nosmoothNom = nosmoothF.Get("%s/%sactual%d" % (recoObs,signal,m)).Clone()
    nosmoothNom.SetDirectory(0)
    nosmoothNom.SetLineWidth(2)
    nosmoothNom.SetLineColor(kBlack)

    smoothNom = smoothF.Get("%s/%ssmooth%d" % (recoObs,signal,m)).Clone()
    smoothNom.SetDirectory(0)
    smoothNom.SetLineWidth(2)
    
    nosmoothMorphed = nosmoothF.Get("%s/%s%d" % (recoObs, signal, m)).Clone()

    nosmoothMorphed.SetDirectory(0)
    nosmoothMorphed.SetLineColor(kBlue)
    nosmoothMorphed.SetFillColor(kWhite)
    nosmoothMorphed.SetLineWidth(2)

    nosmoothNomNoErrs = nosmoothNom.Clone("noerrs_mt%d" % m)
    for b in xrange(1, nosmoothNomNoErrs.GetNbinsX()+1):
        nosmoothNomNoErrs.SetBinError(b, 0.)

    resNosmoothMorphed = nosmoothMorphed.Clone()
    resNosmoothMorphed.SetDirectory(0)
    resNosmoothMorphed.Add(nosmoothNomNoErrs, -1)
    ratioNosmoothMorphed = nosmoothMorphed.Clone()
    ratioNosmoothMorphed.SetDirectory(0)
    ratioNosmoothMorphed.Divide(nosmoothNomNoErrs)

    resSmooth = smoothNom.Clone("residual_smooth_mt%d" % m)
    resSmooth.SetDirectory(0)
    resSmooth.Add(nosmoothNomNoErrs, -1)
    ratioSmooth = smoothNom.Clone("ratio_smooth_mt%d" % m)
    ratioSmooth.SetDirectory(0)
    ratioSmooth.Divide(nosmoothNomNoErrs)

    

    

    maxOldRes = -100
    maxOldRatio = -100
    minOldRes = 100
    minOldRatio = 100

    maxOldRes = resNosmoothMorphed.GetMaximum()
    minOldRes = resNosmoothMorphed.GetMinimum()
    maxOldRatio = ratioNosmoothMorphed.GetMaximum()
    minOldRatio = ratioNosmoothMorphed.GetMinimum()

    paddingOldRes = 0.05 * abs(maxOldRes - minOldRes)
    paddingOldRatio = 0.05 * abs(maxOldRatio - minOldRatio)

    smoothMorphed = smoothF.Get("%s/%s%d" % (recoObs, signal, m)).Clone()

    smoothMorphed.SetDirectory(0)
    smoothMorphed.SetLineColor(kRed)
    smoothMorphed.SetFillColor(kWhite)
    smoothMorphed.SetLineWidth(2)

    resSmoothMorphed = smoothMorphed.Clone()
    resSmoothMorphed.SetDirectory(0)
    resSmoothMorphed.Add(nosmoothNom, -1)
    ratioSmoothMorphed = smoothMorphed.Clone()
    ratioSmoothMorphed.SetDirectory(0)
    ratioSmoothMorphed.Divide(nosmoothNom)

    resLine = TLine(smoothNom.GetXaxis().GetBinLowEdge(1), 0., smoothNom.GetXaxis().GetBinUpEdge(smoothNom.GetNbinsX()), 0.)
    ratioLine = TLine(smoothNom.GetXaxis().GetBinLowEdge(1), 1., smoothNom.GetXaxis().GetBinUpEdge(smoothNom.GetNbinsX()), 1.)

    resLine.SetLineWidth(2)
    ratioLine.SetLineWidth(2)

    maxNewRes = -100
    maxNewRatio = -100
    minNewRes = 100
    minNewRatio = 100

    maxNewRes = max(resSmoothMorphed.GetMaximum(), resSmooth.GetMaximum())
    minNewRes = min(resSmoothMorphed.GetMinimum(), resSmooth.GetMinimum())
    maxNewRatio = max(ratioSmoothMorphed.GetMaximum(), ratioSmooth.GetMaximum())
    minNewRatio = min(ratioSmoothMorphed.GetMinimum(), ratioSmooth.GetMinimum())

    paddingNewRes = 0.05 * abs(maxNewRes - minNewRes)
    paddingNewRatio = 0.05 * abs(maxNewRatio - minNewRatio)


    maxOld = max(nosmoothMorphed.GetMaximum(), nosmoothNom.GetMaximum())
    maxNew = max(smoothMorphed.GetMaximum(), smoothNom.GetMaximum())

    c.cd()
    c.ResetDrawn()
    c.Draw()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()


    ### Comparison ###
    pad1.cd()

    smoothMorphed.SetLineColor(kRed)

    resSmoothMorphed.SetLineColor(kRed)
    ratioSmoothMorphed.SetLineColor(kRed)
    nosmoothNom.Scale(1., "width")
    smoothNom.Scale(1., "width")

    smoothNom.SetLineColor(kGreen+2)
    resSmooth.SetLineColor(kGreen+2)
    ratioSmooth.SetLineColor(kGreen+2)
    #nosmoothNom.SetLineStyle(7)
    l = TLegend(0.7, 0.6, 0.88, 0.88)
    l.SetBorderSize(0)
    l.AddEntry(nosmoothNom, "Actual")
    l.AddEntry(smoothNom, "Smoothed pre-morphing")
    l.AddEntry(nosmoothMorphed, "Morphed: No smoothing")
    l.AddEntry(smoothMorphed, "Morphed: With smoothing")

    nosmoothMorphed.SetTitle("With and without smoothing  m_{t} = %.1f GeV" % (m/10.))
    nosmoothMorphed.Scale(1., "width")
    nosmoothMorphed.GetXaxis().SetTitle("p_{T}(ll) [GeV]")
    nosmoothMorphed.GetXaxis().SetTitleOffset(1.2)
    nosmoothMorphed.GetYaxis().SetTitle("Entries / GeV")
    nosmoothMorphed.GetYaxis().SetTitleOffset(1.3)
    nosmoothMorphed.Draw("hist")
    smoothMorphed.Scale(1, "width")
    smoothMorphed.Draw("hist same")
    
    nosmoothNom.Draw("hist same")
    smoothNom.Draw("hist same")
    l.Draw("same")

    pad2.cd()
    resNosmoothMorphed.GetYaxis().SetRangeUser(min(minOldRes, minNewRes) - max(paddingOldRes, paddingNewRes), max(maxOldRes, maxNewRes) + max(paddingOldRes, paddingNewRes))
    resNosmoothMorphed.SetTitle("")
    resNosmoothMorphed.GetYaxis().SetTitle("residual")
    resNosmoothMorphed.GetYaxis().SetTitleSize(0.07)
    resNosmoothMorphed.GetYaxis().SetTitleOffset(0.6)
    resNosmoothMorphed.Draw("hist%s" % ("" if args.noerrors else " e1"))
    resSmoothMorphed.Draw("hist%s same" % ("" if args.noerrors else " e1"))
    resSmooth.Draw("hist%s same" % ("" if args.noerrors else " e1"))
    resLine.Draw("same")

    pad3.cd()
    ratioNosmoothMorphed.GetYaxis().SetRangeUser(min(minOldRatio, minNewRatio) - max(paddingOldRatio, paddingNewRatio), max(maxOldRatio, maxNewRatio) + max(paddingOldRatio, paddingNewRatio))
    ratioNosmoothMorphed.SetTitle("")
    ratioNosmoothMorphed.GetYaxis().SetTitle("ratio")
    ratioNosmoothMorphed.GetYaxis().SetTitleSize(0.07)
    ratioNosmoothMorphed.GetYaxis().SetTitleOffset(0.6)
    ratioNosmoothMorphed.Draw("hist%s" % ("" if args.noerrors else " e1"))
    ratioSmoothMorphed.Draw("hist%s same" % ("" if args.noerrors else " e1"))
    ratioSmooth.Draw("hist%s same" % ("" if args.noerrors else " e1"))
    ratioLine.Draw("same")

    c.SaveAs("%s/comparison_smoothing_mt%d.png" % (args.outDir, m))


    canvasCompare.cd()
    canvasCompare.ResetDrawn()
    canvasCompare.Draw()
    compPad1.Draw()
    compPad2.Draw()

    compPad1.cd()
    ratioNewNosmoothMorphed = smoothMorphed.Clone("ratioNewNosmoothMorphed")
    ratioNewNosmoothMorphed.SetDirectory(0)
    ratioNewNosmoothMorphed.Divide(nosmoothMorphed)
    ratioNewNosmoothMorphed.GetYaxis().SetRangeUser(ratioNewNosmoothMorphed.GetBinContent(ratioNewNosmoothMorphed.GetMinimumBin())*0.95, ratioNewNosmoothMorphed.GetBinContent(ratioNewNosmoothMorphed.GetMaximumBin())*1.05)
    ratioNewNosmoothMorphed.SetTitle("New / Old %s Up  m_{t} = %.1f GeV  %s" % (args.syst, m/10., "actual" if args.useActual else "morphed"))
    ratioNewNosmoothMorphed.SetLineColor(kRed)


    compPad1.cd()
    ratioNewNosmoothMorphed.Draw("hist%s" % ("" if args.noerrors else " e1"))
    ratioLine.Draw("same")

    compPad2.cd()
    ratioLine.Draw("same")
    canvasCompare.SaveAs("%s/ratio_smoothing_mt%d.png" % (args.outDir, m))

    
nosmoothF.Close()
smoothF.Close()

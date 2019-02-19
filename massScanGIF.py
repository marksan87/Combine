#!/usr/bin/env python
import ROOT
from ROOT import *
import os
import sys
from argparse import ArgumentParser
from datetime import datetime

gROOT.SetBatch(True)
obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

observables = obsTitle.keys()
ttOnlySysts = ['hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

exclusiveStats = {"tt":ttOnlySysts, "tW":tWOnlySysts}

systematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS']

oneSidedSysts = ["toppt", "CRerdON", "CRGluon", "CRQCD", "DS", "amcanlo", "madgraph", "herwigpp" ]

padRatio = 0.25
padOverlap = 0.15
padGap = 0.08

H = 600;
W = 800;


# references for T, B, L, R                                                                                                             
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.1*W


#############################################
#  Quiet                                    #
#  Usage: Quiets info, warnings, or errors  #
#                                           #
#  Ex: TCanvas.c1.SaveAs("myplot.png")      #
#      Quiet(c1.SaveAs)("myplot.png")       #
#############################################
def Quiet(func, level = ROOT.kInfo + 1):
    def qfunc(*args,**kwargs):
        oldlevel = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = level
        try:
            return func(*args,**kwargs)
        finally:
            ROOT.gErrorIgnoreLevel = oldlevel
    return qfunc



parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input template root file")
parser.add_argument("-s", "--sig", default="tt", choices=["tt","tW"], help="signal to plot")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="rec level")
parser.add_argument("--cut", default=0, type=int, help="upper cut for name only")
parser.add_argument("-b", "--rebin", default=1, type=int, help="rebin in increments of 2 GeV")
parser.add_argument("--obs", default="ptll", choices=observables, help="observable to plot")
parser.add_argument("--useActual", action="store_true", default=False, help="use actual templates instead of morphed templates")
parser.add_argument("--noScaling", action="store_true", default=False, help="don't normalize templates to unity")
parser.add_argument("--syst", default="", help="systematic to plot")
parser.add_argument("-o", "--outF", default="", help="output mass scan gif file")

args = parser.parse_args()

if args.syst != "":
    if args.syst not in systematics:
        print "Invalid systematic %s! Choose from:" % args.syst
        print systematics
        sys.exit()

    else:
        if (args.sig == "tt" and args.syst in tWOnlySysts) or \
           (args.sig == "tW" and args.syst in ttOnlySysts):
            print "%s %s not available!" % (args.sig,args.syst)
            sys.exit()

print "\nCreating gif for: %s %s %s %s %s" % (args.sig, args.reco, args.obs, "%s" % ("nominal" if args.syst == "" else "%s" % args.syst), "actual" if args.useActual else "morph" ) 

if args.outF == "":
    args.outF = "massScan_%s_%s_%s%s_%s%s%s.gif" % (args.sig, args.reco, args.obs, "" if args.syst == "" else ("_%s" % args.syst), "actual" if args.useActual else "morph", "" if args.cut == 0 else "_cut%s"%args.cut, "" if args.rebin == 1 else "_rebin%d"%args.rebin)

c = TCanvas("c","c",800,600)

if args.useActual:
    masses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
else:
    masses = [1665 + i for i in xrange(121)]

hists = {}
up = {}
down = {}

ratioUp = {}
ratioDown = {}
f = TFile.Open(args.inF, 'r')
print "Loading templates from file %s.." % args.inF,



maxY = -1
minRatioY = 10
maxRatioY = -10
for m in masses:
    hists[m] = f.Get("%s_%s/%s%s%d" % (args.reco,args.obs,args.sig,"actual" if args.useActual else "",m)).Clone()
    hists[m].SetDirectory(0)
    hists[m].SetLineColor(kBlack)
    hists[m].SetLineWidth(2)


    if not args.noScaling:
        hists[m].Scale(1./hists[m].Integral())
    if args.rebin > 1:
        hists[m].Rebin(args.rebin)
    maxY = max(hists[m].GetMaximum(), maxY)

    if args.syst == "":
        up = None
        down = None

    else:
        up[m] = f.Get("%s_%s/%s%s%d_%sUp" % (args.reco,args.obs,args.sig,"actual" if args.useActual else "",m,args.syst)).Clone()
        up[m].SetDirectory(0)
        up[m].SetLineColor(kRed)
        up[m].SetLineWidth(2)
        if not args.noScaling:
            up[m].Scale(1./up[m].Integral())
        maxY = max(up[m].GetMaximum(), maxY)
        
        if args.syst not in oneSidedSysts:
            down[m] = f.Get("%s_%s/%s%s%d_%sDown" % (args.reco,args.obs,args.sig,"actual" if args.useActual else "",m,args.syst)).Clone()
            down[m].SetDirectory(0)
            down[m].SetLineColor(kBlue)
            down[m].SetLineWidth(2)
            if not args.noScaling:
                down[m].Scale(1./down[m].Integral())
            maxY = max(down[m].GetMaximum(), maxY)
        
        if args.rebin > 1:
            up[m].Rebin(args.rebin)
            if args.syst not in oneSidedSysts:
                down[m].Rebin(args.rebin)



binW = hists[masses[0]].GetBinCenter(2) - hists[masses[0]].GetBinCenter(1)
        
for m in masses:
    if args.syst == "":
        ratioUp = hists[m].Clone()
        ratioUp.Divide(hists[1725])
        minRatioY = min(minRatioY, ratioUp.GetMinimum() - ratioUp.GetBinError(ratioUp.GetMaximumBin()))
        maxRatioY = max(maxRatioY, ratioUp.GetMaximum() + ratioUp.GetBinError(ratioUp.GetMaximumBin()))
    else: 
        ratioUp = up[m].Clone() 
        ratioUp.Divide(hists[m])
        if args.syst not in oneSidedSysts:
            ratioDown = down[m].Clone()
            ratioDown.Divide(hists[m])
            minRatioY = min(minRatioY, ratioDown.GetMinimum() - ratioDown.GetBinError(ratioDown.GetMaximumBin()))
            maxRatioY = max(maxRatioY, ratioDown.GetMaximum() + ratioDown.GetBinError(ratioDown.GetMaximumBin()))


f.Close()
print "done"

canvasRatio = TCanvas('c1Ratio','c1Ratio',W,H)
canvasRatio.SetFillColor(0)
canvasRatio.SetBorderMode(0)
canvasRatio.SetFrameFillStyle(0)
canvasRatio.SetFrameBorderMode(0)
canvasRatio.SetLeftMargin( L/W )
canvasRatio.SetRightMargin( R/W )
canvasRatio.SetTopMargin( T/H )
canvasRatio.SetBottomMargin( B/H )
canvasRatio.SetTickx(0)
canvasRatio.SetTicky(0)
canvasRatio.Draw()
canvasRatio.cd()

pad1 = TPad("zxc_p1","zxc_p1",0,padRatio-padOverlap,1,1)
pad2 = TPad("qwe_p2","qwe_p2",0,0,1,padRatio+padOverlap)
pad1.SetLeftMargin( L/W )
pad1.SetRightMargin( R/W )
pad1.SetTopMargin( T/H/(1-padRatio+padOverlap) )
pad1.SetBottomMargin( (padOverlap+padGap)/(1-padRatio+padOverlap) )
pad2.SetLeftMargin( L/W )
pad2.SetRightMargin( R/W )
pad2.SetTopMargin( (padOverlap)/(padRatio+padOverlap) )
pad2.SetBottomMargin( B/H/(padRatio+padOverlap) )

pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetFrameBorderMode(0)
pad1.SetTickx(0)
pad1.SetTicky(0)

pad2.SetFillColor(0)
pad2.SetFillStyle(4000)
pad2.SetBorderMode(0)
pad2.SetFrameFillStyle(0)
pad2.SetFrameBorderMode(0)
pad2.SetTickx(0)
pad2.SetTicky(0)



tempdir = ".tmpGif_" + datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
os.mkdir(tempdir)
tempfiles = ""
print "Drawing frames...",
sys.stdout.flush()

#for m,h in sorted(hists.items(), key=lambda kv: kv[1]):
for m in masses:
    l = TLegend(0.75,0.75,0.9,0.9)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    
    canvasRatio.cd()
    canvasRatio.ResetDrawn()
    canvasRatio.Draw()
    canvasRatio.cd()

    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    gStyle.SetOptStat(1111)
    hists[1725].SetLineColor(kBlack)
    hists[m].Draw("hist e9")
    if args.syst == "":
        hists[m].SetLineColor(kMagenta)
    hists[m].SetTitle("%s %s %s%s  %s  m_{t} = %.1f GeV" % ("t#bar{t}" if args.sig == "tt" else args.sig,args.reco, obsTitle[args.obs], "" if args.syst == "" else " %s" % args.syst, "actual" if args.useActual else "morphed", m/10.) )
    #hists[m].GetXaxis().SetTitle("%s %s [GeV]" % (args.reco, obsTitle[args.obs]) )
    hists[m].GetYaxis().SetTitle("Entries / %d GeV" % binW )
    hists[m].GetXaxis().SetTitleOffset(1.1)
    hists[m].GetYaxis().SetTitleOffset(1.3)
    hists[m].GetYaxis().SetRangeUser(0.0, 1.05*maxY)
    gPad.Update()


    if args.syst == "":
        hists[1725].Draw("hist same e9")
        l.AddEntry(hists[m], "m_{t} = %.1f" % (m/10.))
        l.AddEntry(hists[1725], "nominal")
    else:
        l.AddEntry(hists[m], "nominal")


    st = hists[m].FindObject("stats")  # TPaveStats object
    st.SetX1NDC(0.7)
    st.SetX2NDC(0.9)
    st.SetY2NDC(0.9)

    if args.syst != "":
        up[m].Draw("hist same e9")
        if args.syst not in oneSidedSysts:
            l.AddEntry(up[m], "%s Up" % args.syst)
            l.AddEntry(down[m], "%s Down" % args.syst)
            down[m].Draw("hist same e9")
        else:
            l.AddEntry(up[m], "%s" % args.syst)

    if args.syst == "" and m != 1725: hists[m].Draw("hist same e9")
    l.Draw("same")

    pad2.cd()
    gStyle.SetOptStat(0)

    if args.syst == "":
        # Mass ratio
        ratioUp = hists[m].Clone()
        ratioUp.Divide(hists[1725])

        # Ratio plot lables
        ratioUp.SetTitle("")
        ratioUp.GetXaxis().SetTitleSize(0.05)
        #print "ratioUp.GetXaxis().GetTitleSize() = %.2f" % ratioUp.GetXaxis().GetTitleSize()
        #ratioUp.GetXaxis().SetTitle(obsTitle[histName[4:]] + " [GeV]")
        ratioUp.GetXaxis().SetTitle("%s %s [GeV]" % (args.reco,obsTitle[args.obs]) )
        ratioUp.GetXaxis().SetTitleOffset(1.1)

        #ratioUp.GetYaxis().SetRangeUser(0.75,1.25)
        ratioUp.GetXaxis().SetTitleSize(0.2)
        ratioUp.GetYaxis().SetRangeUser(0.95*minRatioY,1.05*maxRatioY)
        ratioUp.GetYaxis().SetTitle("ratio wrt nominal")
        
        ratioUp.Draw("hist e9")
        line = TLine(ratioUp.GetXaxis().GetBinLowEdge(1), 1.0, ratioUp.GetXaxis().GetBinUpEdge(ratioUp.GetNbinsX()), 1.0)
        line.SetLineWidth(2)
        line.Draw("same")
    else:
        ratioUp = up[m].Clone()
        ratioUp.Divide(hists[m])

        # Ratio plot lables
        ratioUp.SetTitle("")
        ratioUp.GetXaxis().SetTitleSize(0.06)
        ratioUp.GetXaxis().SetTitle("%s %s [GeV]" % (args.reco,obsTitle[args.obs]) )
        ratioUp.GetXaxis().SetTitleOffset(1.1)
        #ratioUp.GetYaxis().SetRangeUser(0.95*minRatioY,1.05*maxRatioY)
        ratioUp.GetYaxis().SetRangeUser(0.95,1.05)
#        if not args.noScaling:
#            ratioUp.GetYaxis().SetRangeUser(0.8,1.2)
#        else: 
#            if args.syst == "Q2":
#                ratioUp.GetYaxis().SetRangeUser(0.5,1.5)
#            elif args.syst in ["EleIDEff", "EleRecoEff", "PU", "MuIDEff", "MuIsoEff", "MuTrackEff"]:
#                ratioUp.GetYaxis().SetRangeUser(0.95,1.05)
#            else:
#                ratioUp.GetYaxis().SetRangeUser(0.8,1.2)

        ratioUp.GetYaxis().SetTitleSize(0.06)
        ratioUp.GetYaxis().SetTitleOffset(0.5)
        ratioUp.GetYaxis().SetTitle("ratio wrt nominal")
        #ratioUp.GetYaxis().SetTitleOffset(1.1)
        #ratioUp.Draw("hist")
           
        one = hists[m].Clone()
        one.Divide(hists[m])
        if args.syst not in oneSidedSysts:
            ratioDown = down[m].Clone()
            ratioDown.Divide(hists[m])
         #   ratioDown[m].Draw("hist same")
        
       # one.Draw("hist same")
        ratioUp.Draw("hist e9")
        line = TLine(ratioUp.GetXaxis().GetBinLowEdge(1), 1.0, ratioUp.GetXaxis().GetBinUpEdge(ratioUp.GetNbinsX()), 1.0)
        line.SetLineWidth(2)
        line.Draw("same")
        #one.Draw("hist same e9")
        if args.syst not in oneSidedSysts:
            ratioDown.Draw("hist same e9")

    tempfiles += "tmp_%s.png\n" % masses.index(m) 
    Quiet(canvasRatio.SaveAs)(tempdir + "/tmp_%s.png" % masses.index(m))
print "done!"



print "Creating gif %s..." % args.outF,
sys.stdout.flush()

command = "convert -loop 0 -delay %d @file.tx -delay 150 tmp_%d.png %s" % (150 if args.useActual else 20, len(masses)-1, args.outF)
os.system("cd %s" % tempdir + ";echo '%s' > file.tx" % tempfiles + ";" + command + ";mv %s ../; cd .." % args.outF)
print "done!"
os.system("rm -rf %s" % tempdir)


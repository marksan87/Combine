#!/usr/bin/env python
import ROOT
from ROOT import *
import os
import sys
from argparse import ArgumentParser
from datetime import datetime
from ROOT.TMath import Ceil,Log10

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

observables = obsTitle.keys()
ttOnlySysts = ['hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

exclusiveStats = {"tt":ttOnlySysts, "tW":tWOnlySysts}

systematics = ["pileup", "Lumi", "BTagSF", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS']

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

gStyle.SetOptStat(0)

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
parser.add_argument("-a", "--useActual", action="store_true", default=False, help="use actual templates instead of morphed templates")
parser.add_argument("--norm", action="store_true", default=False, help="don't normalize templates to unity")
parser.add_argument("--deltaMT", type=float, default=0.1, help="mass scan delta")
parser.add_argument("--precision", type=int, default=1, help="decimals of precision")
parser.add_argument("--minmt", type=float, default=166.5, help="min mt of scan range")
parser.add_argument("--maxmt", type=float, default=178.5, help="max mt of scan range")
parser.add_argument("-c", "--color", default="kRed", help="color of mt distribution")
parser.add_argument("-p", "--percentDiff", action="store_true", default=False, help="use percent difference instead of ratio")
parser.add_argument("-d", "--saveDebugFrames", action="store_true", default=False, help="save the frame temp directory")
parser.add_argument("--syst", default="", help="systematic to plot")
parser.add_argument("-o", "--outF", default="", help="output mass scan gif file")


try:
    args = parser.parse_args()
except:
    # Print usage
    parser.print_help()
    sys.exit()

exec("color=%s" % (args.color))

normalize = args.norm

if args.precision < 0:
    print "Cannot have %d decimal places! Defaulting to 1" % args.precision
    args.precision = 1

elif args.deltaMT < 10**(-args.precision):
    # Increase precision to leading decimal value in deltaMT
    args.precision = int(Ceil(Log10(1./args.deltaMT)))
#    print "Using precision: %d" % args.precision

precision = args.precision
decimalScaling = 10**precision

morph_masses = range(int(args.minmt*decimalScaling), int((args.maxmt+args.deltaMT)*decimalScaling), int(args.deltaMT*decimalScaling))


gROOT.SetBatch(True)
#if not args.saveDebugFrames:
#    gROOT.SetBatch(True)
#else:
#    gROOT.SetBatch(False)

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
    #masses = [1665, 1695, 1715, nominalMt, 1735, 1755, 1785]
    masses = [166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5]
    masses = [int(m*decimalScaling) for m in masses]
else:
    #masses = [1665 + i for i in xrange(121)]
    masses = morph_masses 

nominalMt = int(172.5*decimalScaling)
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


    if normalize: 
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
        if normalize: 
            up[m].Scale(1./up[m].Integral())
        maxY = max(up[m].GetMaximum(), maxY)
        
        if args.syst not in oneSidedSysts:
            down[m] = f.Get("%s_%s/%s%s%d_%sDown" % (args.reco,args.obs,args.sig,"actual" if args.useActual else "",m,args.syst)).Clone()
            down[m].SetDirectory(0)
            down[m].SetLineColor(kBlue)
            down[m].SetLineWidth(2)
            if normalize: 
                down[m].Scale(1./down[m].Integral())
            maxY = max(down[m].GetMaximum(), maxY)
        
        if args.rebin > 1:
            up[m].Rebin(args.rebin)
            if args.syst not in oneSidedSysts:
                down[m].Rebin(args.rebin)


print "done"
useVariableBinning = None
variableBins = None

# Check for variable binning
hist = hists[masses[0]]
if hist.GetNbinsX() <= 1:
    useVariableBinning = False
else:
    useVariableBinning = False
    binW = hist.GetBinCenter(2) - hist.GetBinCenter(1)

    # Test for variable binning
    for _bin in xrange(2,hist.GetNbinsX()):
        nextBinW = hist.GetBinCenter(_bin+1) - hist.GetBinCenter(_bin)
        if abs(binW - nextBinW) > 1e-4:
            # Variable binning detected
            print "Variable binning found!"
            useVariableBinning = True
            # Don't rebin if variable binning:
            rebin=1
            break

if useVariableBinning:
    _bins = hist.GetXaxis().GetXbins()
    variableBins = [_bins[_b] for _b in xrange(len(_bins))]
#    print "Variable binning detected:"
#    print variableBins

    for m in masses:
        hists[m].Scale(1., "width")
        if args.syst != "":
            up[m].Scale(1., "width")
            if m in down:
                down[m].Scale(1., "width")
else:
    binW = hists[masses[0]].GetBinCenter(2) - hists[masses[0]].GetBinCenter(1)
        
for m in masses:
    if args.syst == "":
        ratioUp[m] = hists[m].Clone()
        ratioUp[m].SetDirectory(0)
        if args.percentDiff and m != nominalMt:
            ratioUp[m].Add(hists[nominalMt], -1)
        ratioUp[m].Divide(hists[nominalMt])
        
        #minRatioY = min(minRatioY, ratioUp[m].GetMinimum() - ratioUp[m].GetBinError(ratioUp[m].GetMaximumBin()))
        #maxRatioY = max(maxRatioY, ratioUp[m].GetMaximum() + ratioUp[m].GetBinError(ratioUp[m].GetMaximumBin()))
        if m != nominalMt:
            minRatioY = min(minRatioY, ratioUp[m].GetMinimum())
            maxRatioY = max(maxRatioY, ratioUp[m].GetMaximum())
    else: 
        ratioUp[m] = up[m].Clone() 
        ratioUp[m].SetDirectory(0)
        if args.percentDiff and m != nominalMt:
            ratioUp[m].Add(hists[nominalMt], -1)
        ratioUp[m].Divide(hists[m])
        if args.percentDiff:
            if m != nominalMt:
                minRatioY = min(minRatioY, ratioUp[m].GetMinimum())
                maxRatioY = max(maxRatioY, ratioUp[m].GetMaximum())
        else:
            minRatioY = min(minRatioY, ratioUp[m].GetMinimum())
            maxRatioY = max(maxRatioY, ratioUp[m].GetMaximum())
        
        if args.syst not in oneSidedSysts:
            ratioDown[m] = down[m].Clone()
            ratioDown[m].SetDirectory(0)
            if args.percentDiff and m != nominalMt:
                ratioDown[m].Add(hists[nominalMt], -1)
            ratioDown[m].Divide(hists[m])
            #minRatioY = min(minRatioY, ratioDown[m].GetMinimum() - ratioDown[m].GetBinError(ratioDown[m].GetMaximumBin()))
            #maxRatioY = max(maxRatioY, ratioDown[m].GetMaximum() + ratioDown[m].GetBinError(ratioDown[m].GetMaximumBin()))
            if args.percentDiff:
                if m != nominalMt:
                    minRatioY = min(minRatioY, ratioDown[m].GetMinimum())
                    maxRatioY = max(maxRatioY, ratioDown[m].GetMaximum())
            else:
                minRatioY = min(minRatioY, ratioDown[m].GetMinimum())
                maxRatioY = max(maxRatioY, ratioDown[m].GetMaximum())


f.Close()

# Ratio pad Y axis range
ratioYpadding = 0.05 * (maxRatioY - minRatioY)
minRatioY -= ratioYpadding
maxRatioY += ratioYpadding


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
    l = TLegend(0.75,0.75,0.88,0.88)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    
    canvasRatio.cd()
    canvasRatio.ResetDrawn()
    canvasRatio.Draw()
    canvasRatio.cd()

    pad1.Draw()
    pad2.Draw()

    pad1.cd()

#    gStyle.SetOptStat(1111)
    hists[nominalMt].SetLineColor(kBlack)
    hists[m].Draw("hist e9")
    if args.syst == "":# and m != nominalMt:
        hists[m].SetLineColor(color)
    
    exec('massStr = "%.' + str(precision) + 'f" % (float(m)/decimalScaling)')
    hists[m].SetTitle("%s %s %s%s  %s  m_{t} = %s GeV" % ("t#bar{t}" if args.sig == "tt" else args.sig,args.reco, obsTitle[args.obs], "" if args.syst == "" else " %s" % args.syst, "actual" if args.useActual else "morphed", massStr) )
    #hists[m].GetXaxis().SetTitle("%s %s [GeV]" % (args.reco, obsTitle[args.obs]) )
    hists[m].GetYaxis().SetTitle("Normalized Entries / %s" % ("%.0f GeV" % binW if not useVariableBinning else "Bin") )
    hists[m].GetXaxis().SetTitleOffset(1.1)
    hists[m].GetYaxis().SetTitleOffset(1.3)
    hists[m].GetYaxis().SetRangeUser(0.0, 1.05*maxY)
    #gPad.Update()
    pad1.Update()
    

    if args.syst == "":
        l.AddEntry(hists[nominalMt], "nominal")
        if m != nominalMt:
            hists[nominalMt].Draw("hist same e9")
            l.AddEntry(hists[m], "m_{t} = %s" % massStr)
    else:
        l.AddEntry(hists[m], "nominal")


#    st = hists[m].FindObject("stats")  # TPaveStats object
#    st.SetX1NDC(0.7)
#    st.SetX2NDC(0.9)
#    st.SetY2NDC(0.9)

    if args.syst != "":
        up[m].Draw("hist same e9")
        if args.syst not in oneSidedSysts:
            l.AddEntry(up[m], "%s Up" % args.syst)
            l.AddEntry(down[m], "%s Down" % args.syst)
            down[m].Draw("hist same e9")
        else:
            l.AddEntry(up[m], "%s" % args.syst)

    if args.syst == "" and m != nominalMt: hists[m].Draw("hist same e9")
    l.Draw("same")

    pad2.cd()
#    gStyle.SetOptStat(0)

    if args.syst == "":
        # Mass ratio
    #    ratioUp[m] = hists[m].Clone()
    #    ratioUp[m].Divide(hists[nominalMt])

        # Ratio plot lables
        ratioUp[m].SetLineColor(color)
        ratioUp[m].SetTitle("")
        ratioUp[m].GetXaxis().SetTitleSize(0.1)
        #print "ratioUp[m].GetXaxis().GetTitleSize() = %.2f" % ratioUp[m].GetXaxis().GetTitleSize()
        #ratioUp[m].GetXaxis().SetTitle(obsTitle[histName[4:]] + " [GeV]")
        ratioUp[m].GetXaxis().SetTitle("%s %s [GeV]" % (args.reco,obsTitle[args.obs]) )
        ratioUp[m].GetXaxis().SetTitleOffset(1.1)

        #ratioUp[m].GetYaxis().SetRangeUser(0.75,1.25)
        ratioUp[m].GetYaxis().SetTitleSize(0.05)
        ratioUp[m].GetYaxis().SetTitleOffset(0.8)
        ratioUp[m].GetYaxis().SetRangeUser(minRatioY,maxRatioY)
        ratioUp[m].GetYaxis().SetTitle("%s wrt nominal" % ("percent diff" if args.percentDiff else "ratio") )
        
        if not (args.percentDiff and m == nominalMt):
            ratioUp[m].Draw("hist e9")
#        if args.percentDiff:
#            if m != nominalMt:
#                ratioUp[m].Draw("hist e9")
#        else:
#            ratioUp[m].Draw("hist e9")
        lineY = 0. if args.percentDiff else 1.
        line = TLine(ratioUp[m].GetXaxis().GetBinLowEdge(1), lineY, ratioUp[m].GetXaxis().GetBinUpEdge(ratioUp[m].GetNbinsX()), lineY)
        line.SetLineWidth(2)
        line.Draw("same")
    else:
        ratioUp[m] = up[m].Clone()
        ratioUp[m].Divide(hists[m])

        # Ratio plot lables
        ratioUp[m].SetTitle("")
        ratioUp[m].GetXaxis().SetTitleSize(0.03)
        ratioUp[m].GetXaxis().SetTitle("%s %s [GeV]" % (args.reco,obsTitle[args.obs]) )
        ratioUp[m].GetXaxis().SetTitleOffset(1.1)
        ratioUp[m].GetYaxis().SetRangeUser(0.95,1.05)
#        if not args.noScaling:
#            ratioUp[m].GetYaxis().SetRangeUser(0.8,1.2)
#        else: 
#            if args.syst == "Q2":
#                ratioUp[m].GetYaxis().SetRangeUser(0.5,1.5)
#            elif args.syst in ["EleIDEff", "EleRecoEff", "PU", "MuIDEff", "MuIsoEff", "MuTrackEff"]:
#                ratioUp[m].GetYaxis().SetRangeUser(0.95,1.05)
#            else:
#                ratioUp[m].GetYaxis().SetRangeUser(0.8,1.2)

        ratioUp[m].GetYaxis().SetTitleSize(0.1)
        ratioUp[m].GetYaxis().SetTitleOffset(0.5)
        ratioUp[m].GetYaxis().SetTitle("ratio wrt nominal")
        #ratioUp[m].GetYaxis().SetTitleOffset(1.1)
        #ratioUp[m].Draw("hist")
           
        one = hists[m].Clone()
        one.Divide(hists[m])
        #if args.syst not in oneSidedSysts:
            #ratioDown[m] = down[m].Clone()
           # ratioDown[m].Divide(hists[m])
         #   ratioDown[m].Draw("hist same")
        
       # one.Draw("hist same")
        if args.percentDiff and m != nominalMt:
            ratioUp[m].Draw("hist e9")
        
        lineY = 0. if args.percentDiff else 1.
        line = TLine(ratioUp[m].GetXaxis().GetBinLowEdge(1), lineY, ratioUp[m].GetXaxis().GetBinUpEdge(ratioUp[m].GetNbinsX()), lineY)
        line.SetLineWidth(2)
        line.Draw("same")
        #one.Draw("hist same e9")
        if args.syst not in oneSidedSysts:
            ratioDown[m].Draw("hist same e9")

    tempfiles += "tmp_%s.png\n" % masses.index(m) 
    Quiet(canvasRatio.SaveAs)(tempdir + "/tmp_%s.png" % masses.index(m))
print "done!"



print "Creating gif %s ..." % args.outF,
sys.stdout.flush()

#command = "convert -loop 0 -delay %d @file.tx -delay 200 tmp_%d.png %s" % (150 if args.useActual else 15, len(masses)-1, args.outF)
command = "convert -loop 0 -delay %d @file.tx -delay 200 tmp_%d.png %s" % (150 if args.useActual else (15 if precision == 1 else 2), len(masses)-1, args.outF)
os.system("cd %s" % tempdir + ";echo '%s' > file.tx" % tempfiles + ";" + command + ";mv %s ../; cd .." % args.outF)
print "done!\n"

if not args.saveDebugFrames:
    os.system("rm -rf %s" % tempdir)
else:
    print "Debug frames located in %s" % tempdir

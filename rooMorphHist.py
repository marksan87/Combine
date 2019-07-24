#!/usr/bin/env python
#from ROOT import kRed,kOrange,kYellow,kGreen,kCyan,kBlue,kViolet,RooWorkspace,RooDataHist,TFile,TH1D,TH1F,RooArgList,TVectorD,RooFit
from ROOT import *
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument("-i", "--inF", help="Input CH templates file")
parser.add_argument("--nosmooth", action="store_true", default=False, help="morph original templates instead of smoothed")
parser.add_argument("-o", "--outF", default="", help="Output root file")
args = parser.parse_args()

if args.outF == "":
    if args.inF.find("mtTemplatesForCH.root") >= 0:
        outF = args.inF.replace("mtTemplatesForCH.root", "RooMomentMorph.root")
    else:
        outF = "result_rooMomentMorph.root"
else:
    outF = args.outF

recoObs = "rec_ptll"

deltaM = 1
morph_masses = range(1665,1786, deltaM)

w = RooWorkspace("w",False)
colors = {166.5:kRed, 169.5:kOrange, 171.5:kYellow, 172.5:kGreen, 173.5:kCyan, 175.5:kBlue, 178.5:kViolet}
masses = [166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5]
#ptllMinMax = RooRealVar("ptllMinMax", "ptllMinMax", 20,220)
#getattr(w,'import')(ptllMinMax)
w.factory('ptllMinMax[20,220]')
ptllMinMax = w.var('ptllMinMax')

ptllFrame = ptllMinMax.frame()

ralptll = RooArgList(ptllMinMax)

f = TFile.Open(args.inF, "read")

hist = {}
dh = {}
pdf = {}
for m in masses:
    m10 = "%d" % (m*10)
    hist[m] = f.Get("%s/tt%s%s" % (recoObs, "actual" if args.nosmooth else "smooth",m10))
    hist[m].SetDirectory(0)

    #hist[m].Scale(1.0 / hist[m].Integral())

    #dh[m] = RooDataHist('dh%s' % m10, 'dh mt = %.1f' % m, RooArgList(ptllMinMax), RooFit.Import(h[m]))
    dh[m] = RooDataHist('dh%s' % m10, 'dh mt = %.1f' % m, RooArgList(ptllMinMax), hist[m])
    #dh[m].plotOn(ptllFrame, RooFit.LineColor(colors[m]))
    pdf[m] = RooHistPdf('pdf%s' % m10, 'pdf mt = %.1f' % m, RooArgSet(ptllMinMax), dh[m])


otherMorph = {}
otherMeans = []
otherMeanErrs = [] 
for m in morph_masses:
    otherMorph[m] = f.Get("%s/tt%d" % (recoObs, m))
    otherMorph[m].SetDirectory(0)
    otherMeans.append(otherMorph[m].GetMean())
    otherMeanErrs.append(otherMorph[m].GetMeanError())

f.Close()

#ptllFrame.Draw()
pdfs = RooArgList()
paramVec = TVectorD(len(masses))

for i,m in enumerate(masses):
    m10 = "%d" % (m*10)
    exec("getattr(w,'import')(pdf[%.1f])"%m)
    
    p = w.pdf('pdf%s'%m10)
    print "m = %.1f\t" % m
    #print pdf
    #pdf.plotOn(ptllFrame)
    pdfs.add(p)
    paramVec[i] = m


morphingMethod = "RooMomentMorph.Linear"
exec("setting = %s" % morphingMethod)
#setting = RooMomentMorph.NonLinear
#setting = RooMomentMorph.Linear
#setting = RooMomentMorph.NonLinearLinFractions
#setting = RooMomentMorph.NonLinearPosFractions
#setting = RooMomentMorph.SineLinear
mt = w.factory("mt[166.5,178.5]")
morph = RooMomentMorph('morph','morph',mt,ralptll,pdfs,paramVec,setting)
#morph.useHorizontalMorphing(False)
#morph.useHorizontalMorphing(False)
getattr(w,'import')(morph)
nbins = 100
morphed = {}
means = []
meanerrs = []
for m in morph_masses: 
    morphMt = m/10.
    mt.setVal(morphMt)
    
    morphed[m] = morph.createHistogram("morphed_tt%d" % m, ptllMinMax, RooFit.Binning(nbins))
    morphed[m].SetDirectory(0)
    morphed[m].SetTitle("p_{T}(ll) m_{t} = %.1f" % morphMt)
    means.append(morphed[m].GetMean())
    meanerrs.append(morphed[m].GetMeanError())
#for mass in sorted(morphed.keys()):
#    morphed[mass].SetLineColor(int(mass*10)-1565)
#    if mass == 166.5:
#        morphed[mass].Draw("hist")
#    else:
#        morphed[mass].Draw("hist same")


h = TH2D("morph_points", "Morphed p_{T}(ll) Templates", len(morph_masses), (morph_masses[0] - 0.5*deltaM) / 10.0, (morph_masses[-1] + 0.5*deltaM)/10.0, morphed[morph_masses[0]].GetNbinsX(), morphed[morph_masses[0]].GetXaxis().GetXmin(), morphed[morph_masses[0]].GetXaxis().GetXmax())
h.GetXaxis().SetTitle("m_{t} [GeV]")
h.GetYaxis().SetTitle("p_{T}(ll) [GeV]")
h.GetZaxis().SetTitle("Events")

h.GetXaxis().SetTitleOffset(1.6)
h.GetYaxis().SetTitleOffset(1.8)
h.GetZaxis().SetTitleOffset(1.5)

for m in morph_masses:
    for b in xrange(1, morphed[m].GetNbinsX() + 1):
        h.Fill(m / 10.0, morphed[m].GetBinCenter(b), morphed[m].GetBinContent(b))

#h.Draw("lego2z")

from array import array
g = TGraph(len(means), array('d',[m/10. for m in morph_masses]), array('d',means))
#g = TGraphErrors(len(means), array('d',[m/10. for m in morph_masses]), array('d',means), array('d',[0.]*len(means)), array('d',meanerrs))
#g = TGraphErrors(len(otherMeans), array('d',[m/10. for m in morph_masses]), array('d',otherMeans), array('d',[0.]*len(otherMeans)), array('d',otherMeanErrs))


# Set morphing method

g.SetName("morphingMomentGraph")
g.SetTitle("%s %s" % (recoObs, morphingMethod))
#g.SetTitle("RooMomentMorph.SineLinear")
#g.SetTitle("RooMomentMorph.NonLinearLinFractions")
#g.SetTitle("RooMomentMorph.NonLinearPosFractions")
#g.SetTitle("Per-bin morphing")
g.GetXaxis().SetTitle("m_{t} [GeV]")
g.GetYaxis().SetTitle("Mean")


outputFile = TFile.Open(outF, "recreate")
g.Write()
h.Write()
for m in morph_masses:
    morphed[m].Write()

outputFile.Close()
print "Output saved to %s" % outF
g.Fit("pol1")
g.Draw("alp")


#mt.setVal(168.2)
#createHistogram=getattr(RooMomentMorph,'createHistogram')
#hist = morph.createHistogram("foo",ptllMinMax,RooFit.Binning(100))

#hist.Draw("HIST")

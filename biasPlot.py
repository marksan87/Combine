#!/usr/bin/env python
import ROOT
from ROOT import *
from glob import glob
import sys
from pprint import pprint

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

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


def biasPlot(vals, obs, masses):
    gr = TGraphAsymmErrors()
    gStyle.SetMarkerColor(kBlue)
    gStyle.SetLineColor(kBlack)
    for i, m in enumerate(masses):
        mass = str(m)
        np = gr.GetN()
        gr.SetPoint(np, m/10., vals[i][0]/10.)
        gr.SetPointError(np, 0., 0., vals[i][1]/10., vals[i][1]/10.)
        # Bias
        #gr.SetPoint(np, m/10., (vals[i][0]-m)/10.)
        #gr.SetPointError(np, 0., 0., vals[i][1]/(10.*vals[i][2]), vals[i][1]/(10.*vals[i][2]))
        #gr.GetYaxis().SetRangeUser(-0.1, 0.1)
        gr.SetMarkerStyle(22)
        gr.Fit('pol1')


    #gr.SetTitle("%s pseudoexperiments" % obsTitle[obs])
    gr.SetTitle(obsTitle[obs])
    outF = TFile.Open("%s_bias.root"%obs, "RECREATE")

    c1 = TCanvas('c1','c1', 1200, 800)
    gStyle.SetTitleFontSize(0.04)
    c1.SetLeftMargin(0.15)
    c1.SetRightMargin(0.25)
    c1.SetBottomMargin(0.25)
    pad1=TPad('p1','p1',0.0,0.0,1.0,0.97)
    pad1.SetTopMargin(0.1)
    pad1.Draw()
    

    f = gr.GetFunction("pol1")
    p0 = f.GetParameter(0)
    p0_err = f.GetParError(0)
    p1 = f.GetParameter(1)
    p1_err = f.GetParError(1)
    fitline = "t#bar{t}   offset %.f #pm %.f   slope %.3f #pm %.2f" % (p0, p0_err, p1, p1_err)


    pad1.cd()
    gr.GetXaxis().SetTitle("Expected m_{t} [GeV]")
    #gr.GetYaxis().SetTitle("Bias [GeV]")
    gr.GetYaxis().SetTitle("Measured m_{t} [GeV]")
    gr.GetXaxis().SetTitleOffset(1.3)
    gr.GetYaxis().SetTitleOffset(1.3)
    gr.Draw("AP")
    txt=TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(43)
    txt.SetTextSize(18)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary} %3.1f fb^{-1} (13 TeV)' % (35.9) )
    txt.DrawLatex(0.55,0.92, fitline)
    #txt.DrawLatex(0.55,0.88, fitline)
    c1.SaveAs("%s_bias.png"%obs)
   
    gr.Write("%s_bias"%obs)
    outF.Close()



obs = "Ep_Em"
masses = [1675, 1685, 1695, 1705, 1715, 1725, 1735, 1745, 1755, 1765, 1775]
trees = {}
mtVals = []

for m in masses:
    gROOT.SetBatch(True)
    fn = glob("toys/%s/higgsCombineTest*mH%d.*.root"%(obs,m))[0]
    print "mt = %d  File found: %s"%(m, fn),
    f = TFile.Open(fn, "READ")
    tree = f.Get("limit")
    entries = tree.GetEntriesFast()
    print "\tEntries:", entries 
    Quiet(tree.Draw)("MT>>MT%d"%m)
    h = gDirectory.Get("MT%d"%m)
    mtVals.append([h.GetMean(), h.GetRMS(), entries**0.5])
    f.Close()

gROOT.SetBatch(False)

biasPlot(mtVals, obs, masses)
        



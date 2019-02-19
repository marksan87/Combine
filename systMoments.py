#!/usr/bin/env python
import ROOT
from ROOT import *
import os
import sys
from array import array
from pprint import pprint
from argparse import ArgumentParser
import gzip
import pickle
from pprint import pprint
import json

gStyle.SetOptStat(0)
obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

try:
    breakpoint
except NameError:
    from pdb import set_trace as breakpoint


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
            ROOT


def calcMoments(h):
    sums = {0:0., 1:0., 2:0., 3:0., 4:0., 6:0., 8:0.}
    kvals = sorted(sums.keys())
    
    for b in xrange(1,h.GetNbinsX()+1):
        for k in kvals: 
            sums[k] += h.GetBinContent(b) * h.GetBinCenter(b)**k

    for k in kvals[1:]:
        sums[k] /= sums[0]

    neff = h.GetEffectiveEntries()
    m1 = sums[1]
    m1err = ((sums[2] - sums[1]**2)/neff)**0.5
    m2 = sums[2]
    m2err = ((sums[4] - sums[2]**2)/neff)**0.5
    m3 = sums[3]
    m3err = ((sums[6] - sums[3]**2)/neff)**0.5
    m4 = sums[4]
    m4err = ((sums[8] - sums[4]**2)/neff)**0.5
    

    return m1,m1err,m2,m2err,m3,m3err,m4,m4err


def massFromMoment(fit, moment, momentError = None):
    offset = fit.GetParameter(0)
    slope = fit.GetParameter(1)

    if momentError is not None:
        return (moment - offset) / slope, momentError / slope
    else:
        return (moment - offset) / slope

class Calibrator:
    slope = 1.
    slopeError = 0.
    offset = 0.
    offsetError = 0.
    masses = []
    moments = []
    momentErrors = []
    calibG = None

    def __init__(self, masses = [], moments = [], momentErrors = [], debugPrint=""): 
        self.masses = masses
        self.moments = moments
        self.momentErrors = momentErrors
        self.calculate(debugPrint)

    def calculate(self, debugPrint=""):
        # Determine calibration coefficients
        if debugPrint != "": print debugPrint
        
        momentG = TGraphErrors(len(self.masses), array('d', [m/10. for m in self.masses]), array('d', self.moments), array('d', [0.]*len(self.masses)), array('d', self.momentErrors))
        momentG.Fit("pol1", "Q")
        momentFit = momentG.GetFunction("pol1") 

        momentOffset = momentFit.GetParameter(0)
        momentSlope = momentFit.GetParameter(1)

#        print "moment: slope = %.2f\toffset = %.2f" % (momentSlope, momentOffset)

        extractedMasses = []
        extractedMassErrors = []
        for i,m10 in enumerate(self.masses):
            m = m10/10.
            moment = self.moments[i]  
            momentError = self.momentErrors[i]
#            print "mt = %.1f\tmoment = %.2f\terror = %.2f" % (m,moment,momentError)
            extractedMasses.append((moment - momentOffset) / momentSlope)
            extractedMassErrors.append(momentError / momentSlope)

        print "extractedMasses =", extractedMasses
        print "extractedMassErrors =", extractedMassErrors
#        sys.exit()

        self.calibG = TGraphErrors(len(self.masses), array('d', [m/10. for m in self.masses]), array('d', extractedMasses), array('d', [0.]*len(self.masses)), array('d', extractedMassErrors))
        self.calibG.GetXaxis().SetTitle("m_{t}^{actual} [GeV]")
        self.calibG.GetXaxis().SetTitleOffset(1.1)
        self.calibG.GetYaxis().SetTitle("m_{t}^{extracted} [GeV]")
        self.calibG.GetYaxis().SetTitleOffset(1.1)

        self.calibG.Fit("pol1", "Q")
        calibFit = self.calibG.GetFunction("pol1")

        self.offset = calibFit.GetParameter(0)
        self.offsetError = calibFit.GetParError(0)
        self.slope = calibFit.GetParameter(1)
        self.slopeError = calibFit.GetParError(1)


    def calibrate(self, mass):
        return (mass - self.offset) / self.slope

    def calibrateError(self, error):
        return error / self.slope


separateSystSamples = ['isr','fsr','DS','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
oneSidedSysts = ["toppt", "CRerdON", "CRGluon", "CRQCD", "DS", "amcanlo", "madgraph", "herwigpp" ]
#systematics = ["Q2", "pileup"]
#oneSidedSysts = ["toppt"]

ttOnlySysts = ['toppt','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

observables = ["rec_ptll"]
#signal = ["tt", "tW"]
signal = ["tt"]
tWactual = [1695,1725,1755]
mtactual = [1665,1695,1715,1725,1735,1755,1785]
mtmorphed = [1665 + i for i in xrange(121)]
masses = {}
masses["tt"] = {"actual":mtactual, "morph":mtmorphed}
masses["tW"] = {"actual":tWactual, "morph":mtmorphed}


def main():
    systematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS']
    c = TCanvas("c1","c1",1200,800)
    
    parser = ArgumentParser()
    parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input root file")
    parser.add_argument("-o", "--outDir", default="momentPlots", help="output directory")
    parser.add_argument("-f", "--interp", default="", help="interpolation method used (for label only!)")
    parser.add_argument("-j", "--json", default="impacts.json", help="output json file for moment impacts")
    parser.add_argument("--nocalibration", action="store_true", default=False, help="don't apply calibration to extracted masses")
    parser.add_argument("--nogensysts", action="store_true", default=False, help="don't include alternate MC generator systs")
    parser.add_argument("--noplots", action="store_true", default=False, help="only create json file, don't make any plots")
    parser.add_argument("--rebin", default=1, help="rebin value")
    parser.add_argument("--debug", action="store_true", default=False, help="break out after loading templates")
    parser.add_argument("--test", action="store_true", default=False, help="only draw one plot for each type")
    args = parser.parse_args()

    if args.outDir[-1] == "/": args.outDir = args.outDir[:-1]
    os.system("mkdir -p %s/histplots" % args.outDir)
    os.system("mkdir -p %s/calibration" % args.outDir)

    if args.test:
        # Only make one systematic plot
        print "Running in test mode"
        #systematics = ["Q2"]
        systematics = ["toppt"]

    else:
        gROOT.SetBatch(True)

    if args.nogensysts:
        print "Not including alternate MC generators: herwigpp, madgraph, amcanlo"
        for syst in ["herwigpp","madgraph","amcanlo"]:
            systematics.remove(syst)

    f = TFile.Open(args.inF, "read")

    print "Loading templates from %s" % args.inF

    global h,moments,g,uncalibrated_mt,calibrated_mt,fitLines
    h = {}
    moments = {}
    for s in signal:
        h[s] = {"actual":{},"morph":{}}
        moments[s] = {"actual":{},"morph":{}}
        for obs in observables:
            h[s]['actual'][obs] = {"nominal":{}}
            h[s]['morph'][obs] = {"nominal":{}}
            #moments[obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[]} }
            moments[s]['actual'][obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }
            moments[s]['morph'][obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }

            for lvl in ['morph','actual']:
                for m in masses[s][lvl]: 
                    h[s][lvl][obs]["nominal"][m] = f.Get("%s/%s%s%d" % (obs,s,"actual" if lvl == "actual" else "",m)).Clone()
                    h[s][lvl][obs]["nominal"][m].SetDirectory(0)
                    
                    m1,m1err,m2,m2err,m3,m3err,m4,m4err = calcMoments(h[s][lvl][obs]["nominal"][m])
                    moments[s][lvl][obs]["nominal"]["m1"].append(m1)
                    moments[s][lvl][obs]["nominal"]["m1err"].append(m1err)
                    moments[s][lvl][obs]["nominal"]["m2"].append(m2)
                    moments[s][lvl][obs]["nominal"]["m2err"].append(m2err)
                    moments[s][lvl][obs]["nominal"]["m3"].append(m3)
                    moments[s][lvl][obs]["nominal"]["m3err"].append(m3err)
                    moments[s][lvl][obs]["nominal"]["m4"].append(m4)
                    moments[s][lvl][obs]["nominal"]["m4err"].append(m4err)



                    for syst in systematics:
                        if lvl == "actual" and syst in separateSystSamples and m != 1725: continue
                       
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        if syst not in h[s][lvl][obs].keys():
                            h[s][lvl][obs][syst] = {"Up":{}}
                            moments[s][lvl][obs][syst] = {"Up":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }
                            if syst not in oneSidedSysts:
                                h[s][lvl][obs][syst]["Down"] = {}
                                moments[s][lvl][obs][syst]["Down"] = {"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} 

                        h[s][lvl][obs][syst]["Up"][m] = f.Get("%s/%s%s%d_%sUp" % (obs,s,"actual" if lvl == "actual" else "",m,syst)).Clone()
                        h[s][lvl][obs][syst]["Up"][m].SetDirectory(0)
                        m1,m1err,m2,m2err,m3,m3err,m4,m4err = calcMoments(h[s][lvl][obs][syst]["Up"][m])
                        moments[s][lvl][obs][syst]["Up"]["m1"].append(m1)
                        moments[s][lvl][obs][syst]["Up"]["m1err"].append(m1err)
                        moments[s][lvl][obs][syst]["Up"]["m2"].append(m2)
                        moments[s][lvl][obs][syst]["Up"]["m2err"].append(m2err)
                        moments[s][lvl][obs][syst]["Up"]["m3"].append(m3)
                        moments[s][lvl][obs][syst]["Up"]["m3err"].append(m3err)
                        moments[s][lvl][obs][syst]["Up"]["m4"].append(m4)
                        moments[s][lvl][obs][syst]["Up"]["m4err"].append(m4err)

                        if syst not in oneSidedSysts:
                            h[s][lvl][obs][syst]["Down"][m] = f.Get("%s/%s%s%d_%sDown" % (obs,s,"actual" if lvl == "actual" else "",m,syst)).Clone()
                            h[s][lvl][obs][syst]["Down"][m].SetDirectory(0)
                            m1,m1err,m2,m2err,m3,m3err,m4,m4err = calcMoments(h[s][lvl][obs][syst]["Down"][m])
                            moments[s][lvl][obs][syst]["Down"]["m1"].append(m1)
                            moments[s][lvl][obs][syst]["Down"]["m1err"].append(m1err)
                            moments[s][lvl][obs][syst]["Down"]["m2"].append(m2)
                            moments[s][lvl][obs][syst]["Down"]["m2err"].append(m2err)
                            moments[s][lvl][obs][syst]["Down"]["m3"].append(m3)
                            moments[s][lvl][obs][syst]["Down"]["m3err"].append(m3err)
                            moments[s][lvl][obs][syst]["Down"]["m4"].append(m4)
                            moments[s][lvl][obs][syst]["Down"]["m4err"].append(m4err)

    f.Close()

    if not args.nocalibration:
        global calib
        calib = {}
        for s in signal:
            calib[s] = {}
            for lvl in ['actual','morph']:
                calib[s][lvl] = {}
                for obs in observables:
                    calib[s][lvl][obs] = {}
                    for moment in ["m1","m2","m3","m4"]:
                        calib[s][lvl][obs][moment] = Calibrator(masses[s][lvl], moments[s][lvl][obs]["nominal"][moment], moments[s][lvl][obs]["nominal"][moment+"err"],debugPrint="Calibration of %s %s %s %s" % (s,lvl,obs,moment))


    if args.debug:
        return

    if not args.noplots:
        txt=TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(43)
        txt.SetTextSize(20)
        txt.SetTextAlign(12)
        #txt.DrawLatex(0.1,0.92,'#bf{CMS} #it{Work in Progress}  %3.1f fb^{-1} (13 TeV)' % (35.9) )
        

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



    g = {}
   
    for s in signal:
        g[s]= {"actual":{},"morph":{}}
       
        if not args.noplots:
            for obs in observables:
                #for m in mtactual:
                for m in masses[s]['actual']: 
                    canvasRatio.cd()
                    canvasRatio.ResetDrawn()
                    canvasRatio.Draw()
                    canvasRatio.cd()

                    pad1.Draw()
                    pad2.Draw()

                    pad1.cd()
                    
                    l = TLegend(0.7,0.7, 0.85,0.88)
                    actualH = h[s]['actual'][obs]['nominal'][m].Clone()
                    actualH.Rebin(args.rebin)
                    actualH.SetLineColor(kBlack)
                    actualH.SetLineWidth(2)
                    morphH = h[s]['morph'][obs]['nominal'][m].Clone()
                    morphH.Scale(actualH.Integral()/morphH.Integral())
                    morphH.SetLineColor(kGreen-2)
                    morphH.SetLineWidth(2)
                    morphH.Rebin(args.rebin)

                    l.AddEntry(actualH, "actual")
                    l.AddEntry(morphH, "morphed")
                        
                    maxH = {actualH:actualH.GetMaximum(), morphH:morphH.GetMaximum()}
                    maxHsorted = sorted(maxH.items(), key=lambda kv: kv[1])
                    maxHsorted.reverse()

                    histsSorted = [hist[0] for hist in maxHsorted if hist[0] is not None]
                    for i,hist in enumerate(histsSorted):
                        if i == 0:
                            hist.Draw("hist 9")
                            hist.SetTitle("%s  %s%s  Actual vs Morphed templates  m_{t} = %.1f GeV" % ("t#bar{t}" if s == "tt" else "ST_tW","" if args.interp == "" else args.interp+"  ",obs, m/10.))
                        else:
                            hist.Draw("hist 9 same")

                    l.Draw("same")
                    
                    pad2.cd()
                    ratio = morphH.Clone()
                    ratio.Divide(actualH)   
                    ratio.SetTitle("")
                    ratio.GetYaxis().SetTitle("morphed / actual")
                    ratio.GetYaxis().SetTitleOffset(1.2)
                    #ratio.GetYaxis().SetTitleSize(0.2)

                    ratio.Draw("hist e9")
                    line = TLine(ratio.GetXaxis().GetBinLowEdge(1), 1.0, ratio.GetXaxis().GetBinUpEdge(ratio.GetNbinsX()), 1.0)
                    line.SetLineWidth(2)
                    line.Draw("same")
                    canvasRatio.SaveAs("%s/histplots/%s_comp_actual_vs_morphed_%s_mt%d.png" % (args.outDir,"TTbar" if s == "tt" else "ST_tW",obs,m))
            
            
            # Histograms
            for lvl in ['actual','morph']:
                for obs in observables:
                    nom = h[s][lvl][obs]["nominal"][1725].Clone()
#        print "Before rebinning: ", nom.GetMaximum()
                    nom.Rebin(args.rebin)
#        print " After rebinning: ", nom.GetMaximum()
                    nom.SetLineColor(kBlack)
                    nom.SetLineWidth(2)
                    for syst in systematics:
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        up = h[s][lvl][obs][syst]["Up"][1725].Clone()
                        up.Rebin(args.rebin)
                        up.SetLineColor(kRed)
                        up.SetLineWidth(2)
                        
                        try:
                            down = h[s][lvl][obs][syst]["Down"][1725].Clone()
                            down.Rebin(args.rebin)
                            down.SetLineColor(kBlue)
                            down.SetLineWidth(2)
                        except KeyError:
                            down = None

                        maxH = {up:(0 if up == None else up.GetMaximum()), down:(0 if down == None else down.GetMaximum()), nom:nom.GetMaximum()}
                        maxHsorted = sorted(maxH.items(), key=lambda kv: kv[1])
                        maxHsorted.reverse()

                        histsSorted = [hist[0] for hist in maxHsorted if hist[0] is not None]
                        
                        canvasRatio.cd()
                        canvasRatio.ResetDrawn()
                        canvasRatio.Draw()
                        canvasRatio.cd()

                        pad1.Draw()
                        pad2.Draw()

                        pad1.cd()
                        
                        for i,hist in enumerate(histsSorted):
                            if i == 0:
                                hist.Draw("histe9")
                                hist.SetTitle("%s  %s  %s  %s%s templates" % ("t#bar{t}" if s == "tt" else "ST_tW",obs,syst,"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed"))
                            else:
                                hist.Draw("histe9 same")

                        

                        txt.DrawLatex(0.69, 0.85, "Moment 1")
                        txt.DrawLatex(0.63, 0.8, "#color[2]{Up: %.4f #pm %.4f}" % (up.GetMean(), up.GetMeanError()))
                        txt.DrawLatex(0.611, 0.75,  "Nom: %.4f #pm %.4f" % (nom.GetMean(), nom.GetMeanError()))
                        if down is not None:
                            txt.DrawLatex(0.6, 0.7, "#color[4]{Down: %.4f #pm %.4f}" % (down.GetMean(), down.GetMeanError()))
                        
                        pad2.cd()
                        ratioUp = up.Clone()
                        ratioUp.Divide(nom)   
                        ratioUp.SetTitle("")
                        ratioUp.GetYaxis().SetTitle("ratio wrt nominal")
                        ratioUp.GetYaxis().SetTitleOffset(1.2)
                        ratioUp.GetYaxis().SetTitleSize(0.2)
                        ratioUp.Draw("hist e9")

                        if down is not None:
                            ratioDown = down.Clone()
                            ratioDown.Divide(nom)
                            ratioDown.Draw("hist e9 same")

                        line = TLine(ratioUp.GetXaxis().GetBinLowEdge(1), 1.0, ratioUp.GetXaxis().GetBinUpEdge(ratioUp.GetNbinsX()), 1.0)
                        line.SetLineWidth(2)
                        line.Draw("same")

                        canvasRatio.SaveAs("%s/histplots/%s_%s_hist_%s_%s.png" % (args.outDir,"TTbar" if s == "tt" else "ST_tW",lvl,obs,syst))


        for lvl in ['actual','morph']:
            # Graphs
            N = len(masses[s][lvl])
            mtmasses = [float(m)/10. for m in masses[s][lvl]]
            for obs in observables:
                g[s][lvl][obs] = { "nominal":{} }
                g[s][lvl][obs]["nominal"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m1err"]))
                g[s][lvl][obs]["nominal"]["m1"].SetName("%s_%s_%s_nominal_m1"%(s,lvl,obs))
                g[s][lvl][obs]["nominal"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m2err"]))
                g[s][lvl][obs]["nominal"]["m2"].SetName("%s_%s_%s_nominal_m2"%(s,lvl,obs))
                g[s][lvl][obs]["nominal"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m3err"]))
                g[s][lvl][obs]["nominal"]["m3"].SetName("%s_%s_%s_nominal_m3"%(s,lvl,obs))
                g[s][lvl][obs]["nominal"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m4err"]))
                g[s][lvl][obs]["nominal"]["m4"].SetName("%s_%s_%s_nominal_m4"%(s,lvl,obs))

                for syst in systematics:
                    if s == "tt" and syst in tWOnlySysts: continue
                    if s == "tW" and syst in ttOnlySysts: continue
                    g[s][lvl][obs][syst] = {"Up":{}}
                    g[s][lvl][obs][syst]["Up"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m1err"]))
                    g[s][lvl][obs][syst]["Up"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m2err"]))
                    g[s][lvl][obs][syst]["Up"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m3err"]))
                    g[s][lvl][obs][syst]["Up"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m4err"]))

                    if syst not in oneSidedSysts:
                        g[s][lvl][obs][syst]["Down"] = {}
                        g[s][lvl][obs][syst]["Down"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m1err"]))
                        g[s][lvl][obs][syst]["Down"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m2err"]))
                        g[s][lvl][obs][syst]["Down"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m3err"]))
                        g[s][lvl][obs][syst]["Down"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m4err"]))

    
    if not args.noplots:
        # Write to output root file
        f = TFile.Open("%s/plots.root"%args.outDir, "recreate")

        for s in signal:
            for lvl in ['actual','morph']:
                for obs in observables:
                    g[s][lvl][obs]["nominal"]["m1"].Write("%s_%s_%s_nominal_m1"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m2"].Write("%s_%s_%s_nominal_m2"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m3"].Write("%s_%s_%s_nominal_m3"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m4"].Write("%s_%s_%s_nominal_m4"%(s,lvl,obs))
                    
                    for syst in systematics:
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        g[s][lvl][obs][syst]["Up"]["m1"].Write("%s_%s_%s_%sUp_m1"%(s,lvl,obs,syst))
                        g[s][lvl][obs][syst]["Up"]["m2"].Write("%s_%s_%s_%sUp_m2"%(s,lvl,obs,syst))
                        g[s][lvl][obs][syst]["Up"]["m3"].Write("%s_%s_%s_%sUp_m3"%(s,lvl,obs,syst))
                        g[s][lvl][obs][syst]["Up"]["m4"].Write("%s_%s_%s_%sUp_m4"%(s,lvl,obs,syst))
                        if syst not in oneSidedSysts:
                            g[s][lvl][obs][syst]["Down"]["m1"].Write("%s_%s_%s_%sDown_m1"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Down"]["m2"].Write("%s_%s_%s_%sDown_m2"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Down"]["m3"].Write("%s_%s_%s_%sDown_m3"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Down"]["m4"].Write("%s_%s_%s_%sDown_m4"%(s,lvl,obs,syst))

    with gzip.open("moments.pklz", "wb") as mf:
        pickle.dump(moments, mf, protocol=pickle.HIGHEST_PROTOCOL)


    c.cd()

    if not args.noplots: print "About to draw plots"
    mg = {}
    fitLines = {}
    for s in signal:
        mg[s] = {"actual":{}, "morph":{}}
        fitLines[s] = {"actual":{}, "morph":{}}
        for lvl in ['actual','morph']:
            for obs in observables:
                mg[s][lvl][obs] = {}
                fitLines[s][lvl][obs] = {"nominal":{}}
                
                for moment in ["m1","m2","m3","m4"]:
                    if not args.nocalibration and not args.noplots:
                        # Draw calibration curve
                        reclvl,obsname = obs.split("_")
                        calib[s][lvl][obs][moment].calibG.SetTitle("%s %s  %s  Moment %s  Calibration curve" % (reclvl, obsTitle[obsname], lvl, moment[1]))
                        calib[s][lvl][obs][moment].calibG.Draw("ap")
                        txt.DrawLatex(0.6, 0.35, "Slope")
                        txt.DrawLatex(0.75, 0.35, "Offset")
                        txt.DrawLatex(0.57, 0.3, "%.4f #pm %.4f" % (calib[s][lvl][obs][moment].slope, calib[s][lvl][obs][moment].slopeError))
                        txt.DrawLatex(0.72, 0.3, "%.4f #pm %.4f" % (calib[s][lvl][obs][moment].offset, calib[s][lvl][obs][moment].offsetError))

                        
                        c.SaveAs("%s/calibration/calibration_%s_%s_%s_%s.png" % (args.outDir,s,lvl,obs,moment))
                        
                    # Get nominal fit line
                    g[s][lvl][obs]["nominal"][moment].Fit("pol1", "Q")
                    fitLines[s][lvl][obs]["nominal"][moment] = g[s][lvl][obs]["nominal"][moment].GetFunction("pol1") 
                    fitLines[s][lvl][obs]["nominal"][moment].SetLineColor(kBlack)    
                

                for syst in systematics:
                    if lvl == "actual" and syst in separateSystSamples: continue
                    if s == "tt" and syst in tWOnlySysts: continue
                    if s == "tW" and syst in ttOnlySysts: continue
                    mg[s][lvl][obs][syst] = {}
                    fitLines[s][lvl][obs][syst] = {"Up":{}}
                    if syst not in oneSidedSysts:
                        fitLines[s][lvl][obs][syst]["Down"] = {}
                    
                    for moment in ["m1","m2","m3","m4"]:
                        l = TLegend(0.12,0.7, 0.3,0.88)
                        l.SetBorderSize(0)

#                        l.AddEntry(g[s][lvl][obs]["nominal"][moment], "nominal", "l")
#                        g[s][lvl][obs]["nominal"][moment].Fit("pol1","Q")
#                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetLineColor(kBlack)
#                        nomFit = g[s][lvl][obs]["nominal"][moment].GetFunction("pol1")

                        g[s][lvl][obs][syst]["Up"][moment].Fit("pol1","Q")
                        g[s][lvl][obs][syst]["Up"][moment].SetLineColor(kRed)
                        g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetLineColor(kRed)
#                        upFit = g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1")
                        fitLines[s][lvl][obs][syst]["Up"][moment] = g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1")
                        
                        
                        l.AddEntry(g[s][lvl][obs][syst]["Up"][moment], "%s%s" % (syst, " Up" if syst not in oneSidedSysts else ""), "l")

                        mg[s][lvl][obs][syst][moment] = TMultiGraph()
                        mg[s][lvl][obs][syst][moment].SetName("%s_%s_%s_%s"%(s,obs,syst,moment))

                        mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs]["nominal"][moment].Clone(), "p")
                        mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs][syst]["Up"][moment].Clone(), "p")
                        if syst not in oneSidedSysts:
                            g[s][lvl][obs][syst]["Down"][moment].Fit("pol1","Q")
                            g[s][lvl][obs][syst]["Down"][moment].SetLineColor(kBlue)
                            g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetLineColor(kBlue)
#                            downFit = g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1")
                            fitLines[s][lvl][obs][syst]["Down"][moment] = g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1")
                            mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs][syst]["Down"][moment].Clone(), "p")
                            l.AddEntry(g[s][lvl][obs][syst]["Down"][moment], "%s Down" % syst, "l")

                        
                        if not args.noplots:
                            mg[s][lvl][obs][syst][moment].SetTitle("%s  %s  %s  %s%s  Moment %s" % ("t#bar{t}" if s == "tt" else "St tW",obs,syst,"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed",moment[1]))
                            mg[s][lvl][obs][syst][moment].Draw("a9")

                            mg[s][lvl][obs][syst][moment].GetXaxis().SetTitle("m_{t} [GeV]")
                            mg[s][lvl][obs][syst][moment].GetYaxis().SetTitle("%s [GeV]%s" % (moment, "" if moment == "m1" else "^{%s}"%moment[1]))
                            mg[s][lvl][obs][syst][moment].GetYaxis().SetTitleOffset(1.2)

                            l.Draw("same")
                            txt.DrawLatex(0.44, 0.35, "Moment %s Slopes"%moment[1])
                            txt.DrawLatex(0.42, 0.3, "#color[2]{Up: %.4f #pm %.4f}" % (fitLines[s][lvl][obs][syst]["Up"][moment].GetParameter(1), fitLines[s][lvl][obs][syst]["Up"][moment].GetParError(1)))
                            txt.DrawLatex(0.407, 0.25,  "Nom: %.4f #pm %.4f" % (fitLines[s][lvl][obs]["nominal"][moment].GetParameter(1), fitLines[s][lvl][obs]["nominal"][moment].GetParError(1)))
                            if syst not in oneSidedSysts:
                                txt.DrawLatex(0.4, 0.2, "#color[4]{Down: %.4f #pm %.4f}" % (fitLines[s][lvl][obs][syst]["Down"][moment].GetParameter(1), fitLines[s][lvl][obs][syst]["Down"][moment].GetParError(1)))
                            
                            txt.DrawLatex(0.76, 0.35, "Offsets")
                            txt.DrawLatex(0.72, 0.3, "#color[2]{Up: %.4f #pm %.4f}" % (fitLines[s][lvl][obs][syst]["Up"][moment].GetParameter(0), fitLines[s][lvl][obs][syst]["Up"][moment].GetParError(0)))
                            txt.DrawLatex(0.707, 0.25,  "Nom: %.4f #pm %.4f" % (fitLines[s][lvl][obs]["nominal"][moment].GetParameter(0), fitLines[s][lvl][obs]["nominal"][moment].GetParError(0)))
                            if syst not in oneSidedSysts:
                                txt.DrawLatex(0.7, 0.2, "#color[4]{Down: %.4f #pm %.4f}" % (fitLines[s][lvl][obs][syst]["Down"][moment].GetParameter(0), fitLines[s][lvl][obs][syst]["Down"][moment].GetParError(0)))

                            c.SaveAs("%s/%s_%s_%s_%s_%s.png" % (args.outDir,"TTbar" if s == "tt" else "ST_tW",lvl,obs,syst,moment))

    if not args.noplots:
        f.Close()

    # Extracted masses
    uncalibrated_mt = {}
    calibrated_mt = {}
    for s in ["tt"]:
        uncalibrated_mt[s] = {}
        calibrated_mt[s] = {}

        for lvl in ["actual","morph"]:
            uncalibrated_mt[s][lvl] = {}
            calibrated_mt[s][lvl] = {}


            for obs in observables:
                uncalibrated_mt[s][lvl][obs] = {"nominal":{}} 
                calibrated_mt[s][lvl][obs] = {"nominal":{}} 
                #uncalibrated_mt[s][lvl][obs] = { "nominal":{},"stat":{"Up":{}, "Down":{}} }
                #calibrated_mt[s][lvl][obs] = { "nominal":{},"stat":{"Up":{}, "Down":{}} } 
                

                 
                for moment in ["m1", "m2", "m3", "m4"]:
                    uncalibrated_mt[s][lvl][obs]["nominal"][moment] = {"mass":{}, "error":{}}
                    if not args.nocalibration: 
                        calibrated_mt[s][lvl][obs]["nominal"][moment] = {"mass":{}, "error":{}}
                
                    # 172.5 value is midpoint of list
                    nomEntry = len(moments[s][lvl][obs]["nominal"][moment]) // 2

                    # Nominal moment and statistical uncertainty
                    nom_mt, nom_mterr = massFromMoment(fitLines[s][lvl][obs]["nominal"][moment], moments[s][lvl][obs]["nominal"][moment][nomEntry], moments[s][lvl][obs]["nominal"][moment+"err"][nomEntry])
                    uncalibrated_mt[s][lvl][obs]["nominal"][moment]["mass"] = nom_mt
                    uncalibrated_mt[s][lvl][obs]["nominal"][moment]["error"] = nom_mterr

                    if not args.nocalibration:
                        calibrated_nom_mt = calib[s][lvl][obs][moment].calibrate(nom_mt)
                        calibrated_mt[s][lvl][obs]["nominal"][moment]["mass"] = calibrated_nom_mt
                        calibrated_mt[s][lvl][obs]["nominal"][moment]["error"] = calib[s][lvl][obs][moment].calibrateError(nom_mterr)


                    for syst in systematics:
                        if lvl == "actual" and syst in separateSystSamples: continue
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        if syst not in uncalibrated_mt[s][lvl][obs]:
                            uncalibrated_mt[s][lvl][obs][syst] = {"Up":{}}
                            if not args.nocalibration:
                                calibrated_mt[s][lvl][obs][syst] = {"Up":{}}
                            #uncalibrated_mt[s][lvl][obs][syst] = {"Up":{moment:{}}}
                            #calibrated_mt[s][lvl][obs][syst] = {"Up":{moment:{}}}
                            
                            if sys not in oneSidedSysts:
                                uncalibrated_mt[s][lvl][obs][syst]["Down"] = {}
                                if not args.nocalibration:
                                    calibrated_mt[s][lvl][obs][syst]["Down"] = {}
                                #uncalibrated_mt[s][lvl][obs][syst]["Down"] = {moment:{}}
                                #calibrated_mt[s][lvl][obs][syst]["Down"] = {moment:{}}

                        uncalibrated_mt[s][lvl][obs][syst]["Up"][moment] = {} 
                        if not args.nocalibration:
                            calibrated_mt[s][lvl][obs][syst]["Up"][moment] = {} 
                        if sys not in oneSidedSysts:
                            uncalibrated_mt[s][lvl][obs][syst]["Down"][moment] = {} 
                            if not args.nocalibration:
                                calibrated_mt[s][lvl][obs][syst]["Down"][moment] = {} 

                       
                        nomEntry = len(moments[s][lvl][obs][syst]["Up"][moment]) // 2
#                        print "Now on %s %s %s" % (lvl,syst,moment)
                        
                        try:
                            #mtUp = massFromMoment(fitLines[s][lvl][obs][syst]["Up"][moment], moments[s][lvl][obs][syst]["Up"][moment][nomEntry] )
                            #mtUp = massFromMoment(fitLines[s][lvl][obs]['nominal'][moment], moments[s][lvl][obs][syst]['Up'][moment][nomEntry] )
                            # Fit nominal moment to systematic fit line
                            mtUp = massFromMoment(fitLines[s][lvl][obs][syst]["Up"][moment], moments[s][lvl][obs]['nominal'][moment][nomEntry] )
                        except ZeroDivisionError:
                            # Skip
                            continue

                        except Exception as e:
                            pprint(e)
                            print "s =", s
                            print "lvl =", lvl
                            print "obs =", obs
                            print "syst =", syst
                            print "moment =", moment
                            print "nomEntry =", nomEntry
                            sys.exit()
                        
                        uncalibrated_mt[s][lvl][obs][syst]["Up"][moment]["mass"] = mtUp
                        uncalibrated_mt[s][lvl][obs][syst]["Up"][moment]["impact"] = mtUp - nom_mt

                        if not args.nocalibration:
                            calibrated_mtUp = calib[s][lvl][obs][moment].calibrate(mtUp)
                            calibrated_mt[s][lvl][obs][syst]["Up"][moment]["mass"] = calibrated_mtUp
                            calibrated_mt[s][lvl][obs][syst]["Up"][moment]["impact"] = calibrated_mtUp - calibrated_nom_mt
                        
                        if syst not in oneSidedSysts:
                            #mtDown = massFromMoment(fitLines[s][lvl][obs][syst]["Down"][moment], moments[s][lvl][obs][syst]["Down"][moment][nomEntry])
                            
                            mtDown = massFromMoment(fitLines[s][lvl][obs][syst]["Down"][moment], moments[s][lvl][obs]['nominal'][moment][nomEntry])
                            #mtDown = massFromMoment(fitLines[s][lvl][obs]['nominal'][moment], moments[s][lvl][obs][syst]['Down'][moment][nomEntry] )
                            uncalibrated_mt[s][lvl][obs][syst]["Down"][moment]["mass"] = mtDown
                            uncalibrated_mt[s][lvl][obs][syst]["Down"][moment]["impact"] = nom_mt - mtDown

                            if not args.nocalibration:
                                calibrated_mtDown = calib[s][lvl][obs][moment].calibrate(mtDown)
                                calibrated_mt[s][lvl][obs][syst]["Down"][moment]["mass"] = calibrated_mtDown
                                calibrated_mt[s][lvl][obs][syst]["Down"][moment]["impact"] = calibrated_nom_mt - calibrated_mtDown
                        
#        sys.exit()

    # Write moment json file
    for moment in ["m1","m2","m3","m4"]:
        if args.nocalibration:
            nominal_mt = uncalibrated_mt['tt']['morph']['rec_ptll']['nominal'][moment]['mass']
            stat_err = uncalibrated_mt['tt']['morph']['rec_ptll']['nominal'][moment]['error']
        else:
            nominal_mt = calibrated_mt['tt']['morph']['rec_ptll']['nominal'][moment]['mass']
            stat_err = calibrated_mt['tt']['morph']['rec_ptll']['nominal'][moment]['error']

        nominal_mt *= 10.
        stat_err *= 10.

        jsondata = {"POIs":[ {"fit":[nominal_mt - stat_err , nominal_mt, nominal_mt + stat_err], "name":"MT"} ] }
        jsondata["params"] = []
        for syst in systematics:
            if syst == "DS": continue
            if args.nocalibration:
                mtUp = uncalibrated_mt['tt']['morph']['rec_ptll'][syst]["Up"][moment]["mass"]
            else:
                mtUp = calibrated_mt['tt']['morph']['rec_ptll'][syst]["Up"][moment]["mass"]
            if syst not in oneSidedSysts:
                if args.nocalibration:
                    mtDown = uncalibrated_mt['tt']['morph']['rec_ptll'][syst]["Down"][moment]["mass"]
                else:
                    mtDown = calibrated_mt['tt']['morph']['rec_ptll'][syst]["Down"][moment]["mass"]
            else:
                mtDown = nominal_mt / 10.
            
            mtDown *= 10.
            mtUp *= 10.

            info = {}
            info["MT"] = [mtDown,nominal_mt,mtUp]
            info["impact_MT"] = max(abs(nominal_mt - mtDown), abs(nominal_mt - mtUp))
            info["name"] = syst
            info["groups"] = []
            info["type"] = "Gaussian"
            info["fit"] = [0.,0.,0.]
            info["prefit"] = [-1., 0., 1.]


            jsondata["params"].append(info)

        outF = "%s_%s" % (moment,args.json)
        with open(outF, "w") as f:
            json.dump(jsondata, f, indent=4, separators=(',', ': '))

        print "Impact data on moment %s written to %s" % (moment,outF)

if __name__ == "__main__":
    sys.exit(main())



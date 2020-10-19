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
try:
    from Style import *
except ImportError:
    # Ingore if in job
    pass

try: 
    import CMS_lumi
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText = "Work in Progress"
except ImportError:
    # Ignore if in job
    pass

gStyle.SetOptStat(0)
gStyle.SetLegendTextSize(0.03)
#_diffDists = ["ptll_M0_E0", "ptll_M0_E1", "ptll_M0_E2", "ptll_M1_E0", "ptll_M1_E1", "ptll_M1_E2", "ptll_M2_E0", "ptll_M2_E1", "ptll_M2_E2"]
_diffDists = ["rec_ptll_M0_E0", "rec_ptll_M0_E1", "rec_ptll_M0_E2", "rec_ptll_M1_E0", "rec_ptll_M1_E1", "rec_ptll_M1_E2", "rec_ptll_M2_E0", "rec_ptll_M2_E1", "rec_ptll_M2_E2"]
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

try:
    breakpoint
except NameError:
    from pdb import set_trace as breakpoint


#thestyle = Style()
#
#HasCMSStyle = False
#style = None
#if os.path.isfile('tdrstyle.C'):
#    ROOT.gROOT.ProcessLine('.L tdrstyle.C')
#    ROOT.setTDRStyle()
#    print "Found tdrstyle.C file, using this style."
#    HasCMSStyle = True
#    if os.path.isfile('CMSTopStyle.cc'):
#        gROOT.ProcessLine('.L CMSTopStyle.cc+')
#        style = CMSTopStyle()
#        style.setupICHEPv1()
#        print "Found CMSTopStyle.cc file, use TOP style if requested in xml file."
#if not HasCMSStyle:
#    print "Using default style defined in cuy package."
#    thestyle.SetStyle()

#ROOT.gROOT.ForceStyle()





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

massColor = {1665:kRed, 1695:kOrange, 1715:kGreen, 1725:kBlack, 1735:(kAzure+8), 1755:kBlue, 1785:kViolet}
MEscaleColor = {"MEscale1":kRed, "MEscale2":kOrange, "MEscale3":kGreen, "MEscale4":(kAzure+8), "MEscale5":kBlue, "MEscale6":kViolet} 

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
    try:
        m4err = ((sums[8] - sums[4]**2)/neff)**0.5
    except ValueError:
        m4err = 0.0

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

#        print "extractedMasses =", extractedMasses
#        print "extractedMassErrors =", extractedMassErrors

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

ttOnlySysts = ['toppt','Pdf','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','MEscale1','MEscale2','MEscale3','MEscale4','MEscale5','MEscale6']
tWOnlySysts = ['DS']


# File name
signalName = {"tt":"TTbar", "tW":"ST_tW", "tttW":"sumtttW"}

# Plot title
signalTitle = {"tt":"t#bar{t}", "tW":"ST_tW", "tttW":"t#bar{t} + tW"}
#signal = ["tt"]

def main():
    MEscaleSysts = ['MEscale1','MEscale2','MEscale3','MEscale4','MEscale5','MEscale6']
    systematics = ["pileup", "BkgNorm", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS']#,'MEscale1','MEscale2','MEscale3','MEscale4','MEscale5','MEscale6']
    
    kinematics = ["rec_ptll", "rec_Mll", "rec_ptpos", "rec_Epos", "rec_ptp_ptm", "rec_Ep_Em"]
    observables = ["rec_ptll", "rec_Mll", "rec_ptpos", "rec_Epos", "rec_ptp_ptm", "rec_Ep_Em", "rec_leadLepPt", "rec_leadJetPt"] + _diffDists
    signal = ["tt", "tW", "tttW"]
    backgrounds = ['DY','TTV','Diboson','ST_bkgd'] 
    bkgNormScaling = {\
        "DY"     : 0.300, 
        "ST_bkgd": 0.055, 
        }


    parser = ArgumentParser()
    parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input root file")
    parser.add_argument("-o", "--outDir", default="", help="output directory")
    parser.add_argument("--interp", default="", help="interpolation method used (for label only!)")
    parser.add_argument("-j", "--json", default="impacts.json", help="output json file for moment impacts")
    parser.add_argument("-f", "--format", nargs="+", default=["png"], help="output file format(s)")
    parser.add_argument("--nomoments", action="store_true", default=False, help="ignore moments")
    #parser.add_argument("--applycalibration", action="store_true", default=False, help="apply calibration to extracted masses")
    parser.add_argument("--nogensysts", action="store_true", default=False, help="don't include alternate MC generator systs")
#    parser.add_argument("--obs", default="rec_ptll", choices=(observables+["all"]), help="reco_observable")
    parser.add_argument("-a", "--asimov", action="store_true", default=False, help="use asimov set instead of data_obs for moment mass extraction")
    parser.add_argument("-s", "--sig", default="tt", choices=["tt","tW","tttW"], help="signal in data_obs distribution")
    parser.add_argument("--systs", default="", nargs="*", choices=(systematics + ["none","None"]), help="ONLY plot these systematics") 
    parser.add_argument("--scaleGen", action="store_true", default=False, help="scale gen rates to reco level")
    parser.add_argument("--obs", nargs="+", default=["rec_ptll"])#, choices=(observables+["all","kin","diff"]), help="reco_observable")
    parser.add_argument("--reco", default="rec", choices=["rec","gen"])
    parser.add_argument("--minmt", type=int, default=1665, help="min mt")
    parser.add_argument("--maxmt", type=int, default=1785, help="max mt")
    
    parser.add_argument("--noBinScaling", action="store_true", default=False, help="don't scale by 1/binWidth for variable binned templates")
    parser.add_argument("--mtscan", action="store_true", default=False, help="only do mt scan")
    parser.add_argument("--noplots", action="store_true", default=False, help="only create json file, don't make any plots")
    parser.add_argument("--displayhistmoments", action="store_true", default=False, help="display moments on histograms")
    parser.add_argument("--rebin", default=1, help="rebin value")
    parser.add_argument("--norm", action="store_true", default=False, help="normalize templates to unity")
    parser.add_argument("--debug", action="store_true", default=False, help="break out after loading templates")
    parser.add_argument("--test", action="store_true", default=False, help="only draw one plot for each type")
    args = parser.parse_args()
    
    f = TFile.Open(args.inF, "read")
 
    useDataObs = not args.asimov
    if not args.nomoments:
        if useDataObs:
            print "Moment masses extracted from data_obs"
            if args.sig == "tt":
                print "Assuming data_obs signal is just tt"
            else:
                print "WARNING! tW subtraction not implemented for data_obs yet!"
        else:
            print "Moment masses extracted from asimov set"
    scaleGen = args.scaleGen
    

    #useMomentCalibration = args.applycalibration
    useMomentCalibration = True 
    noBinScaling = args.noBinScaling 
    __mtactual = [1665,1695,1715,1725,1735,1755,1785]
    __tWactual = [1695,1725,1755]
    mtactual = []
    tWactual = []
    for m in __mtactual:
        if m >= args.minmt and m <= args.maxmt:
            mtactual.append(m)
    for m in __tWactual:
        if m >= args.minmt and m <= args.maxmt:
            tWactual.append(m)
    mtmorphed = [mtactual[0] + i for i in xrange(mtactual[-1]-mtactual[0]+1) ]

    masses = {}
    masses["tt"] = {"actual":mtactual, "morph":mtmorphed}
    masses["tW"] = {"actual":tWactual, "morph":mtmorphed}
    masses["tttW"] = {"actual":tWactual, "morph":mtmorphed}

    if args.outDir == "":
        if args.inF.find("mtTemplatesForCH.root") >= 0:
            args.outDir = args.inF.replace("mtTemplatesForCH.root", "plots")
        else:
            args.outDir = "template_plots"
    elif args.outDir[-1] == "/": args.outDir = args.outDir[:-1]
#    if args.obs != "all":
#        observables = [args.obs]
#        os.system("mkdir -p %s/histplots/mtscan" % args.outDir)
#        os.system("mkdir -p %s/rates" % args.outDir)
#        os.system("mkdir -p %s/calibration" % args.outDir)
#    else:
#        for obs in observables:
#            os.system("mkdir -p %s/%s/histplots/mtscan" % (args.outDir,obs))
#            os.system("mkdir -p %s/%s/rates" % (args.outDir,obs))
#            os.system("mkdir -p %s/%s/calibration" % (args.outDir,obs))
    if "kin" in args.obs:
        observables = kinematics
    elif "diff" in args.obs:
        observables = _diffDists
    elif "all" not in args.obs:
        observables = args.obs
 
    if args.reco == "gen":
        observables = [obs.replace("rec","gen") for obs in observables]
    if args.mtscan:
        print "mt scan only"
        signal = ["tt"]

    outFormatStr = ""
    for _fmt in args.format:
        outFormatStr += "%s, " % _fmt
    outFormatStr = outFormatStr[:-2]
    print "Saving plots in the following formats: %s" % outFormatStr

    for obs in observables:
        os.system("mkdir -p %s/%s/histplots/mtscan" % (args.outDir,obs))
        if not args.mtscan:
            os.system("mkdir -p %s/%s/rates" % (args.outDir,obs))
            if useMomentCalibration:
                os.system("mkdir -p %s/%s/calibration" % (args.outDir,obs))


    if "none" in args.systs or "None" in args.systs:
        systematics = []
    elif args.systs != "":
        # Only plot these systematics
        systematics = args.systs
        print "Plotting the following systematics:", systematics
    if args.nomoments:
        print "Ignoring moments"

    if args.test:
        # Only make these systematic plots
        print "Running in test mode"
        #systematics = ["Q2"]
        systematics = ["Q2","toppt"]
        print "Using the following systematics:", systematics

#    else:
#        gROOT.SetBatch(True)
    gROOT.SetBatch(True)

    c = TCanvas("c1","c1",1200,800)
    if args.nogensysts:
        print "Not including alternate MC generators: herwigpp, madgraph, amcanlo"
        for syst in ["herwigpp","madgraph","amcanlo"]:
            systematics.remove(syst)


    print "Loading templates from %s" % args.inF

    variableBinning = {}    # Whether a given observable is using variable binning

    global h,moments,g,uncalibrated_mt,calibrated_mt,fitLines,diff,rates,rateG,rateMG,dataobs,dataobs_bkgsub,total_background,dataobs_moments
    h = {}
    moments = {}
    rates = {}
    dataobs = {}
    dataobs_bkgsub = {} # Dataobs with background subtracted
    dataobs_moments = {}
    total_background = {}
    for b in backgrounds:
        h[b] = {"actual":{}}
        for obs in observables:
            h[b]["actual"][obs] = f.Get("%s/%s" % (obs,b))
            h[b]["actual"][obs].SetDirectory(0)

            if obs not in total_background:
                total_background[obs] = h[b]["actual"][obs].Clone("%s_backgrounds" % obs)
                total_background[obs].SetDirectory(0)
            else:
                total_background[obs].Add(h[b]["actual"][obs])

    for obs in observables:
        dataobs[obs] = f.Get("%s/data_obs" % obs)
        dataobs[obs].SetDirectory(0)

        dataobs_bkgsub[obs] = dataobs[obs].Clone("%s_dataobs_bkgsub" % obs)
        dataobs_bkgsub[obs].Add(total_background[obs], -1)

    if useDataObs and not args.nomoments:
        for obs in observables:
            dataobs_moments[obs] = {}
            m1,m1err,m2,m2err,m3,m3err,m4,m4err = calcMoments(dataobs_bkgsub[obs])
            dataobs_moments[obs]["m1"] = m1
            dataobs_moments[obs]["m1err"] = m1err
            dataobs_moments[obs]["m2"] = m2
            dataobs_moments[obs]["m2err"] = m2err
            dataobs_moments[obs]["m3"] = m3
            dataobs_moments[obs]["m3err"] = m3err
            dataobs_moments[obs]["m4"] = m4
            dataobs_moments[obs]["m4err"] = m4err

    for s in signal:
        h[s] = {"actual":{},"morph":{}}
        moments[s] = {"actual":{},"morph":{}}
        rates[s] = {"actual":{},"morph":{}}
        for obs in observables:
            h[s]['actual'][obs] = {"nominal":{}}
            h[s]['morph'][obs] = {"nominal":{}}
            rates[s]['actual'][obs] = {"nominal":{"rate":[], "rateError":[]}}
            rates[s]['morph'][obs] = {"nominal":{"rate":[], "rateError":[]}}
            #moments[obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[]} }
            moments[s]['actual'][obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }
            moments[s]['morph'][obs] = { "nominal":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }

            for lvl in ['morph','actual']:
                for _i,m in enumerate(masses[s][lvl]):
#                    print "Now loading %s %s %s %d" % (obs,s,lvl,m)
                    h[s][lvl][obs]["nominal"][m] = f.Get("%s/%s%s%d" % (obs,s,"actual" if lvl == "actual" else "",m)).Clone()
                    h[s][lvl][obs]["nominal"][m].SetDirectory(0)
                    if _i == 0:
                        # Test for variable binning
                        _binWidths = [ h[s][lvl][obs]["nominal"][m].GetBinWidth(_bin) for _bin in range(1,h[s][lvl][obs]["nominal"][m].GetNbinsX()+1) ]
                        _binW = _binWidths[0]
                        variableBinning[obs] = False 
                        for _width in _binWidths:
                            if abs(_width - _binW) > 1e-3:
#                                print "**** Variable binning detected for %s ****" % obs
                                variableBinning[obs] = True
                                break
                        
#                        if len(h[s][lvl][obs]["nominal"][m].GetXaxis().GetXbins()) > 0:
#                            variableBinning[obs] = True
#                        else:
#                            variableBinning[obs] = False

                    _error = ROOT.Double(0)
                    _rate = h[s][lvl][obs]["nominal"][m].IntegralAndError(1,h[s][lvl][obs]["nominal"][m].GetNbinsX(),_error)
                    rates[s][lvl][obs]["nominal"]["rate"].append(_rate)
                    rates[s][lvl][obs]["nominal"]["rateError"].append(_error)

                    if args.norm:
                        # Normalize
                        h[s][lvl][obs]["nominal"][m].Scale(1./h[s][lvl][obs]["nominal"][m].Integral())

                    if not noBinScaling and variableBinning[obs]: 
                        # Variable binning. Scale by 1/binWidth
                        print "Scaling %s by bin width" % obs
                        h[s][lvl][obs]["nominal"][m].Scale(1., "width")

                    if not args.nomoments:
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
                        #if lvl == "actual" and syst in separateSystSamples and m != 1725: continue
                       
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        if syst not in h[s][lvl][obs].keys():
                            h[s][lvl][obs][syst] = {"Up":{}}
                            rates[s][lvl][obs][syst] = {"Up":{"rate":[], "rateError":[]}}
                            moments[s][lvl][obs][syst] = {"Up":{"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} }
                            if syst not in oneSidedSysts:
                                h[s][lvl][obs][syst]["Down"] = {}
                                rates[s][lvl][obs][syst]["Down"] = {"rate":[], "rateError":[]}
                                moments[s][lvl][obs][syst]["Down"] = {"m1":[], "m1err":[], "m2":[], "m2err":[], "m3":[], "m3err":[], "m4":[], "m4err":[]} 
                        if syst == "BkgNorm":
                            h[s][lvl][obs][syst]["Up"][m] = h[s][lvl][obs]["nominal"][m].Clone(h[s][lvl][obs]["nominal"][m].GetName()+"_BkgNormUp")
                            h[s][lvl][obs][syst]["Up"][m].SetDirectory(0)
                            
                            
                            #if h[s][lvl][obs][syst]["Up"][m].GetNbinsX() > 10:
                            #    print "h[%s][%s][%s][%s]['Up'][%d].GetNbinsX() = %d" % (s,lvl,obs,syst,m,h[s][lvl][obs][syst]["Up"][m].GetNbinsX())
                    
                            h[s][lvl][obs][syst]["Up"][m].Add(h["DY"]["actual"][obs], -bkgNormScaling["DY"])
                            h[s][lvl][obs][syst]["Up"][m].Add(h["ST_bkgd"]["actual"][obs], -bkgNormScaling["ST_bkgd"])
                        else:
                            h[s][lvl][obs][syst]["Up"][m] = f.Get("%s/%s%s%d_%sUp" % (obs,s,"actual" if lvl == "actual" else "",m,syst)).Clone()
                            h[s][lvl][obs][syst]["Up"][m].SetDirectory(0)

                        _error = ROOT.Double(0)
                        _rate = h[s][lvl][obs][syst]["Up"][m].IntegralAndError(1,h[s][lvl][obs][syst]["Up"][m].GetNbinsX(),_error)
                        rates[s][lvl][obs][syst]["Up"]["rate"].append(_rate)
                        rates[s][lvl][obs][syst]["Up"]["rateError"].append(_error)

                        if args.norm:
                            # Normalize
                            h[s][lvl][obs][syst]["Up"][m].Scale(1./h[s][lvl][obs][syst]["Up"][m].Integral())
                        
                        if not noBinScaling and variableBinning[obs]:
                            h[s][lvl][obs][syst]["Up"][m].Scale(1., "width")

                        if not args.nomoments:
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
                            if syst == "BkgNorm":
                                h[s][lvl][obs][syst]["Down"][m] = h[s][lvl][obs]["nominal"][m].Clone(h[s][lvl][obs]["nominal"][m].GetName()+"_BkgNormDown")
                                h[s][lvl][obs][syst]["Down"][m].SetDirectory(0)
                                h[s][lvl][obs][syst]["Down"][m].Add(h["DY"]["actual"][obs], bkgNormScaling["DY"]/(1+bkgNormScaling["DY"]))
                                h[s][lvl][obs][syst]["Down"][m].Add(h["ST_bkgd"]["actual"][obs], bkgNormScaling["ST_bkgd"]/(1+bkgNormScaling["ST_bkgd"]))
                            else: 
                                h[s][lvl][obs][syst]["Down"][m] = f.Get("%s/%s%s%d_%sDown" % (obs,s,"actual" if lvl == "actual" else "",m,syst)).Clone()
                                h[s][lvl][obs][syst]["Down"][m].SetDirectory(0)
                            _error = ROOT.Double(0)
                            _rate = h[s][lvl][obs][syst]["Down"][m].IntegralAndError(1,h[s][lvl][obs][syst]["Down"][m].GetNbinsX(),_error)
                            rates[s][lvl][obs][syst]["Down"]["rate"].append(_rate)
                            rates[s][lvl][obs][syst]["Down"]["rateError"].append(_error)
                            
                            if args.norm:
                                # Normalize
                                h[s][lvl][obs][syst]["Down"][m].Scale(1./h[s][lvl][obs][syst]["Down"][m].Integral())
                            
                            if not noBinScaling and variableBinning[obs]: 
                                h[s][lvl][obs][syst]["Down"][m].Scale(1., "width")
                           
                            if not args.nomoments:
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


    if useMomentCalibration and not args.nomoments and not args.mtscan:
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

        #CMS_lumi.CMS_lumi(canvasRatio, 4, 11)

    g = {}
    rateG = {}
    for s in signal:
        g[s]= {"actual":{},"morph":{}}
        rateG[s] = {"actual":{},"morph":{}}

        if not args.noplots:
            for obs in observables:
                if args.mtscan: break   # Skip if only doing mt scan
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
                    l.SetBorderSize(0)
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
                            hist.SetTitle("%s  %s%s  Actual vs Morphed templates  m_{t} = %.1f GeV" % (signalTitle[s],"" if args.interp == "" else args.interp+"  ",obs, m/10.))
                            hist.GetYaxis().SetTitle("%sEvents / %s" % ("Normalized " if args.norm else "", "%.0f GeV" % (hist.GetBinCenter(2) - hist.GetBinCenter(1)) if not variableBinning[obs] else "GeV") )
                            hist.GetYaxis().SetTitleOffset(1.3)
                        else:
                            hist.Draw("hist 9 same")

                    l.Draw("same")
                    ####CMS_lumi.CMS_lumi(pad1, 4, 12)
                    CMS_lumi.CMS_lumi(pad1, 4, 11)
                    
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
                    for fmt in args.format:
                        canvasRatio.SaveAs("%s/%s/histplots/%s_comp_actual_vs_morphed_%s_mt%d.%s" % (args.outDir,obs,signalName[s],obs,m,fmt))
            
            print "Now on alternate mass histograms"
            # Alternate mass histograms        
            for lvl in ['actual','morph']:
                for obs in observables:
                    for syst in (["nominal"] + systematics):
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        for var in ["Up","Down"]:
                            if var == "Down" and (syst == "nominal" or syst in oneSidedSysts): continue
                            l = TLegend(0.75, 0.68, 0.88, 0.88)
                            #ratioLegend = TLegend(0.12, 0.12, 0.25, 0.35)
                            #ratioLegend = TLegend(0.3, 0.68, 0.53, 0.88)
                            
                            #ratioLegend = TLegend(0.45, 0.68, 0.65, 0.88)
                            ratioLegend = TLegend(0.45, 0.78, 0.88, 0.88)
                            l.SetBorderSize(0)
                            ratioLegend.SetBorderSize(0)
                            ratioLegend.SetNColumns(4)

                            if syst == "nominal":
                                nom = h[s][lvl][obs]["nominal"][1725].Clone()
                            else:
                                nom = h[s][lvl][obs][syst][var][1725].Clone()
                            nom.Rebin(args.rebin)
                            nom.SetMarkerStyle(20)
                            nom.SetMarkerSize(0)
                            nom.SetLineColor(kBlack)
                            nom.SetLineWidth(2)
#                            nom.Scale(1./nom.Integral())

                            nomNoErr = nom.Clone(nom.GetName()+"_noerr")
                            for b in range(1, nomNoErr.GetNbinsX()+1):
                                nomNoErr.SetBinError(b, 0.)
                            canvasRatio.cd()
                            canvasRatio.ResetDrawn()
                            canvasRatio.Draw()
                            canvasRatio.cd()

                            pad1.Draw()
                            pad2.Draw()

                            pad1.cd()

                            
                            nom.SetTitle("%s %s%s %s Mass Scan" % (obs[:3],obsTitle[obs[4:]], "" if syst == "nominal" else (" " + syst + ( (" "+var) if syst not in oneSidedSysts else "")), lvl ))
                            # Dummy draw for creating axis objects
                            nom.Draw("hist 9")
                            nom.GetYaxis().SetTitle("%sEvents / %s" % ("Normalized " if args.norm else "", "%.0f GeV" % (nom.GetBinCenter(2) - nom.GetBinCenter(1)) if not variableBinning[obs] else "GeV") )
                            nom.GetYaxis().SetTitleOffset(1.4)
                            c.cd()
                            #nom.Draw("hist 9")
                            pad1.cd()
                            mhist = {}
                            ratios = {}

                            #minRatioY = -0.3
                            #maxRatioY = 0.3
                            #minRatioY = 99.9
                            #maxRatioY = -99.9
                            
                            minRatioY = 0.8
                            maxRatioY = 1.2

                            maxY = 1.0

                            #for i,m in enumerate(masses[s][lvl]):
                            for i,m in enumerate(masses[s]['actual']):
                                if m == 1725: 
                                    mhist[m] = nom
                                else:
                                    if syst == "nominal":
                                        mhist[m] = h[s][lvl][obs]['nominal'][m].Clone()
                                        nomRate = h[s][lvl][obs]['nominal'][1725].Integral()
                                    else:
                                        mhist[m] = h[s][lvl][obs][syst][var][m].Clone()
                                        nomRate = h[s][lvl][obs][syst][var][1725].Integral()
                                    mhist[m].SetDirectory(0)
                                    mhist[m].SetLineColor(massColor[m])
                                    mhist[m].SetLineWidth(2)
                                    mhist[m].SetMarkerStyle(20)
                                    mhist[m].SetMarkerSize(0)
                                    
                                    mhist[m].Scale(nomRate/mhist[m].Integral())
                                    maxY = max(maxY, mhist[m].GetMaximum())
                                
                                ratios[m] = mhist[m].Clone(mhist[m].GetName()+"_ratio")
                                ratios[m].Divide(nomNoErr)
                                
                                minRatioY = min(minRatioY, ratios[m].GetMinimum())
                                maxRatioY = max(maxRatioY, ratios[m].GetMaximum())
                                
                                l.AddEntry(mhist[m], "m_{t} = %.1f" % (m/10.))
                                ratioLegend.AddEntry(mhist[m], "m_{t} = %.1f" % (m/10.))
                            
                            nom.GetYaxis().SetRangeUser(0, maxY*1.05)
                            c.cd()
                            nom.Draw("hist 9")
                            pad1.cd()
                            nom.Draw("hist 9")
                            for i,m in enumerate(masses[s]['actual']):
                                if m == 1725: continue
                                c.cd()
                                mhist[m].Draw("hist same 9")
                                l.Draw("same")
                                pad1.cd()
                                mhist[m].Draw("hist same 9")

                            c.cd()

                            #minRatioY = min(minRatioY, 0.7)
                            #maxRatioY = max(maxRatioY, 1.3)

                            #padding = 0.05*(maxRatioY - minRatioY)
                            #padding = 0.15*(maxRatioY - minRatioY)
                            padding = 0.1*(maxRatioY - minRatioY)

                            ###CMS_lumi.CMS_lumi(pad1, 4, 12)
                            CMS_lumi.CMS_lumi(pad1, 4, 11)
                            nom.GetXaxis().SetTitle(obsTitle[obs[4:]] + " [GeV]")
                            c.Update()
                            for fmt in args.format:
                                c.SaveAs("%s/%s/histplots/mtscan/noratio_mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "_"+ syst + (var if (syst not in oneSidedSysts and syst != "nominal") else ""),obs,fmt))
                                #c.SaveAs("%s/%s/histplots/mtscan/noratio_mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "" if syst == "nominal" else ("_"+ syst + (var if syst not in oneSidedSysts else "")),obs,fmt))
                            pad2.cd()
                            nom.GetXaxis().SetTitle("")


                            for i,m in enumerate(masses[s]['actual']):
                                if i == 0:
                                    ratios[m].SetTitle("")
                                    ratios[m].GetXaxis().SetTitle(obsTitle[obs[4:]] + " [GeV]")
                                    #ratios[m].GetYaxis().SetTitle("% ratios to 172.5")
                                    ratios[m].GetYaxis().SetTitle("ratio to 172.5")
                                    ratios[m].GetYaxis().SetRangeUser(minRatioY-padding, maxRatioY+padding)
                                    ratios[m].Draw("hist e1 9")
                                    c.cd()
                                    ratios[m].Draw("hist e1 9")
                                    ratios[m].GetXaxis().SetTitleSize(0.04)
                                    ratios[m].GetXaxis().SetTitleOffset(1.0)
                                    ratios[m].GetYaxis().SetTitleSize(0.04)
                                    ratios[m].GetYaxis().SetTitleOffset(1.0)

                                    pad2.cd()
#                        elif m == 1725:
#                            continue
                                else:
                                    ratios[m].Draw("hist same e1 9")
                                    c.cd()
                                    ratios[m].Draw("hist same e1 9")
                                    pad2.cd()

                            #line = TLine(nom.GetXaxis().GetBinLowEdge(1), 0.0, nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 0.0)
                            #line.SetLineWidth(2)
                            #line.Draw("same")

                            c.cd()
                            ratios[masses[s]['actual'][0]].SetTitle("%s %s%s %s Mass Scan%s" % (obs[:3],obsTitle[obs[4:]], "" if syst == "nominal" else (" " + syst + ( (" "+var) if syst not in oneSidedSysts else "")) ,lvl," %.0f GeV Binning" % (ratios[m].GetBinCenter(2)-ratios[m].GetBinCenter(1)) if not variableBinning[obs] else "") )
                            #line.Draw("same")
                            ratioLegend.Draw("same")
                            c.Update()
                            for fmt in args.format:
                                c.SaveAs("%s/%s/histplots/mtscan/ratio_mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "_"+ syst + (var if (syst not in oneSidedSysts and syst != "nominal") else "") ,obs,fmt))
                                #c.SaveAs("%s/%s/histplots/mtscan/ratio_mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "" if syst == "nominal" else ("_"+ syst + (var if syst not in oneSidedSysts else "")) ,obs,fmt))

                            pad1.cd()
                            ratios[masses[s]['actual'][0]].GetXaxis().SetTitleSize(0.1)
                            ratios[masses[s]['actual'][0]].GetXaxis().SetTitleOffset(1.1)
                            ratios[masses[s]['actual'][0]].GetYaxis().SetTitleSize(0.073)
                            ratios[masses[s]['actual'][0]].GetYaxis().SetTitleOffset(0.5)
                            l.Draw("same")
                            canvasRatio.cd()
                            ratios[masses[s]['actual'][0]].SetTitle("")
                            for fmt in args.format:
                                canvasRatio.SaveAs("%s/%s/histplots/mtscan/mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "_"+ syst + (var if (syst not in oneSidedSysts and syst != "nominal") else ""), obs,fmt))
                                #canvasRatio.SaveAs("%s/%s/histplots/mtscan/mtscan_%s_%s%s_%s.%s" % (args.outDir,obs,signalName[s],lvl, "" if syst == "nominal" else ("_"+ syst + (var if syst not in oneSidedSysts else "")), obs,fmt))

                    
                    nom = h[s][lvl][obs]["nominal"][1725].Clone()
                    nomNoErr = nom.Clone(nom.GetName()+"__noErr")
                    for b in range(1,nomNoErr.GetNbinsX()+1):
                        nomNoErr.SetBinError(b,0.)
                    
                    nomRatio = nom.Clone(nom.GetName()+"_nomratio")
                    nomRatio.Divide(nomNoErr)

                    # Systematic variation plots
                    for syst in systematics:
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        
                        l = TLegend(0.7,0.7, 0.88,0.88)
                        l.SetBorderSize(0)
                        up = h[s][lvl][obs][syst]["Up"][1725].Clone()
                        up.Rebin(args.rebin)
                        up.SetLineColor(kRed)
                        up.SetLineWidth(2)
                        up.SetMarkerStyle(20)
                        up.SetMarkerSize(0)
                        up.Add(nom, -1)
                        #up.Scale(1./up.Integral())

                        try:
                            down = h[s][lvl][obs][syst]["Down"][1725].Clone()
                            down.Rebin(args.rebin)
                            down.SetLineColor(kBlue)
                            down.SetLineWidth(2)
                            down.SetMarkerStyle(20)
                            down.SetMarkerSize(0)
                            down.Add(nom, -1)
                            #down.Scale(1./down.Integral())

                            l.AddEntry(up, "%s Up" % syst)
                            #l.AddEntry(nom, "nominal")
                            l.AddEntry(down, "%s Down" % syst)
                        except KeyError:
                            down = None
                            l.AddEntry(up, "%s" % syst)
                            l.AddEntry(nom, "nominal")


                        #maxH = {up:(0 if up == None else up.GetMaximum()), down:(0 if down == None else down.GetMaximum()), nom:nom.GetMaximum()}
                        maxH = {up:(0 if up == None else up.GetMaximum()), down:(0 if down == None else down.GetMaximum())}
                        maxHsorted = sorted(maxH.items(), key=lambda kv: kv[1])
                        maxHsorted.reverse()
                        
                        minY = up.GetMinimum()
                        maxY = up.GetMaximum()

                        if down is not None:
                            minY = min(minY, down.GetMinimum())
                            maxY = max(maxY, down.GetMaximum())

                        padding = 0.05*(maxY - minY)

                        histsSorted = [hist[0].Clone("_"+hist[0].GetName()) for hist in maxHsorted if hist[0] is not None]
                        
                        canvasRatio.cd()
                        canvasRatio.ResetDrawn()
                        canvasRatio.Draw()
                        canvasRatio.cd()

                        pad1.Draw()
                        pad2.Draw()

                        pad1.cd()
                        
                        for i,hist in enumerate(histsSorted):
                            if i == 0:
                                hist.Draw("hist 9")
                                hist.GetYaxis().SetRangeUser(minY-padding, maxY+padding)
                                #hist.SetTitle("%s  %s  %s  %s%s" % (signalTitle[s],obs,syst,"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed"))
                                hist.SetTitle("%s  %s  %s  %s" % (signalTitle[s],obs,syst,"" if args.interp =="" else args.interp+"  "))
                                hist.GetYaxis().SetTitle("%sEvents / %s" % ("Normalized " if args.norm else "", "%.0f GeV" % (hist.GetBinCenter(2) - hist.GetBinCenter(1)) if not variableBinning[obs] else "GeV") )
                                hist.GetYaxis().SetTitleOffset(1.3)
                            else:
                                hist.Draw("hist 9 same")
                        
                        nom.Draw("hist 9 same")
                        
                        ###CMS_lumi.CMS_lumi(pad1, 4, 12)
                        CMS_lumi.CMS_lumi(pad1, 4, 11)
                        l.Draw("same")

                        pad2.cd()
                        ratioUp = up.Clone()
                        #ratioUp.Add(nom,-1)
                        ratioUp.Divide(nomNoErr)   
                        ratioUp.SetTitle("")
                        #ratioUp.GetYaxis().SetTitle("ratio wrt nominal")
                        #ratioUp.GetYaxis().SetTitleOffset(1.2)
                        #ratioUp.GetYaxis().SetTitleSize(0.2)
                        
                        ratioUp.GetXaxis().SetTitle(obsTitle[obs[4:]] + " [GeV]")
                        ratioUp.GetXaxis().SetTitleSize(0.1)
                        ratioUp.GetXaxis().SetTitleOffset(1.1)
                        ratioUp.GetYaxis().SetTitleSize(0.073)
                        ratioUp.GetYaxis().SetTitleOffset(0.5)
                        ratioUp.GetYaxis().SetTitle("ratio to nominal")

                        #minRatioY = -0.2
                        #maxRatioY = 0.2
                        minRatioY = 0.2
                        maxRatioY = -0.2

                        minRatioY = min(minRatioY, ratioUp.GetMinimum())
                        maxRatioY = max(maxRatioY, ratioUp.GetMaximum())

                        if down is not None:
                            ratioDown = down.Clone()
                            #ratioDown.Add(nom,-1)
                            ratioDown.Divide(nomNoErr)
                            minRatioY = min(minRatioY, ratioDown.GetMinimum())
                            maxRatioY = max(maxRatioY, ratioDown.GetMaximum())

                        
                        #minRatioY = min(minRatioY, 0.0)
                        #maxRatioY = max(maxRatioY, 0.0)

                        padding = 0.05*(maxRatioY - minRatioY)
                        ratioUp.GetYaxis().SetRangeUser(minRatioY-padding, maxRatioY+padding)
                        ratioUp.Draw("hist e9")


                        if down is not None:
                            ratioDown.Draw("hist e9 same")

                        #line = TLine(ratioUp.GetXaxis().GetBinLowEdge(1), 0.0, ratioUp.GetXaxis().GetBinUpEdge(ratioUp.GetNbinsX()), 0.0)
                        #line.SetLineWidth(2)
                        #line.Draw("same")
                        nomRatio.Draw("hist e9 same")
                        for fmt in args.format:
                            canvasRatio.SaveAs("%s/%s/histplots/reldiff_%s_%s_hist_%s_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,syst,fmt))
    

            # Histograms
            for lvl in ['actual','morph']:
                if args.mtscan: break
                for obs in observables:
                    nom = h[s][lvl][obs]["nominal"][1725].Clone()
                    nom.Rebin(args.rebin)
                    nom.SetLineColor(kBlack)
                    nom.SetLineWidth(2)

                    nomNoErr = nom.Clone(nom.GetName()+"__noErr")
                    for b in range(1,nomNoErr.GetNbinsX()+1):
                        nomNoErr.SetBinError(b,0.)
                    
                    nomRatio = nom.Clone(nom.GetName()+"_nomratio")
                    nomRatio.Divide(nomNoErr)
                    # Gen/rec plots
                    #reclvl = obs[:4]
                    #otherlvl = "gen_" if revlvl == "rec_" else "rec_"
                    #otherobs = otherlvl + obs[4:]
                    #if otherobs in observables:
                    if obs[:4] == "rec_" and obs.replace("rec_","gen_") in observables:
                        gen = h[s][lvl][obs.replace("rec_","gen_")]["nominal"][1725]
                        gen.Rebin(args.rebin)
                        gen.SetLineColor(kRed)
                        gen.SetLineWidth(2)

                        rec = nom.Clone("rec_%s" % nom.GetName())

                        if scaleGen:
                            # Scale gen rate to rec rate
                            gen.Scale(rec.Integral() / gen.Integral())

                        l = TLegend(0.75,0.7, 0.88,0.88)
                        l.SetBorderSize(0)
                        l.AddEntry(nom, "rec")
                        l.AddEntry(gen, "gen")

                        ratioGen = gen.Clone("%s_gen_ratio" % obs[4:])
                        ratioGen.Divide(nom)

                        canvasRatio.cd()
                        canvasRatio.ResetDrawn()
                        canvasRatio.Draw()
                        canvasRatio.cd()

                        pad1.Draw()
                        pad2.Draw()

                        pad1.cd()
                        rec.SetTitle("%s  %s  m_{t} = %.1f GeV" % (obsTitle[obs[4:]], signalTitle[s], 172.5))
                        rec.GetYaxis().SetTitle("%sEvents / %s" % ("Normalized " if args.norm else "", "%.0f GeV" % (rec.GetBinCenter(2) - rec.GetBinCenter(1)) if not variableBinning[obs] else "GeV" ) )
                        rec.Draw("hist")
                        gen.Draw("hist same")

                        CMS_lumi.CMS_lumi(pad1, 4, 11)
                        l.Draw("same")

                        pad2.cd()
                        ratioGen.SetTitle("")
                        ratioGen.GetYaxis().SetTitle("gen / rec")

                        ratioGen.GetXaxis().SetTitle(obsTitle[obs[4:]] + " [GeV]")
                        ratioGen.GetXaxis().SetTitleSize(0.1)
                        ratioGen.GetXaxis().SetTitleOffset(1.1)
                        ratioGen.GetYaxis().SetTitleSize(0.073)
                        ratioGen.GetYaxis().SetTitleOffset(0.5)

                        ratioGen.Draw("hist")

                        if scaleGen:
                            line = TLine(ratioGen.GetXaxis().GetBinLowEdge(1), 1.0, ratioGen.GetXaxis().GetBinUpEdge(ratioGen.GetNbinsX()), 1.0)
                            line.SetLineWidth(2)
                            line.Draw("same")

                        for fmt in args.format:
                            canvasRatio.SaveAs("%s/%s/histplots/recgen_%s_%s_hist_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,fmt))

                    for syst in systematics:
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        
                        l = TLegend(0.75,0.7, 0.88,0.88)
                        l.SetBorderSize(0)

                        up = h[s][lvl][obs][syst]["Up"][1725].Clone()
                        up.Rebin(args.rebin)
                        up.SetLineColor(kRed)
                        up.SetLineWidth(2)
                        try:
                            down = h[s][lvl][obs][syst]["Down"][1725].Clone()
                            down.Rebin(args.rebin)
                            down.SetLineColor(kBlue)
                            down.SetLineWidth(2)
                            l.AddEntry(up, "%s Up" % syst)
                            l.AddEntry(nom, "nominal")
                            l.AddEntry(down, "%s Down" % syst)
                        except KeyError:
                            down = None
                            l.AddEntry(up, "%s" % syst)
                            l.AddEntry(nom, "nominal")

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
                                hist.Draw("hist 9")
                                #hist.SetTitle("%s  %s  %s  %s%s" % (signalTitle[s],obs,syst,"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed"))
                                hist.SetTitle("%s  %s  %s  %s" % (signalTitle[s],obsTitle[obs[4:]],syst,"" if args.interp =="" else args.interp+"  "))
                                hist.GetYaxis().SetTitle("%sEvents / %s" % ("Normalized " if args.norm else "", "%.0f GeV" % (hist.GetBinCenter(2) - hist.GetBinCenter(1)) if not variableBinning[obs] else "GeV" ) )
                                hist.GetYaxis().SetTitleOffset(1.4)
                            else:
                                hist.Draw("hist 9 same")

                        nom.Draw("hist 9 same")

                        if args.displayhistmoments:
                            txt.DrawLatex(0.69, 0.85, "Moment 1")
                            txt.DrawLatex(0.63, 0.8, "#color[2]{Up: %.4f #pm %.4f}" % (up.GetMean(), up.GetMeanError()))
                            txt.DrawLatex(0.611, 0.75,  "Nom: %.4f #pm %.4f" % (nom.GetMean(), nom.GetMeanError()))
                            if down is not None:
                                txt.DrawLatex(0.6, 0.7, "#color[4]{Down: %.4f #pm %.4f}" % (down.GetMean(), down.GetMeanError()))
                        
                        #CMS_lumi.CMS_lumi(pad1, 4, 12)
                        CMS_lumi.CMS_lumi(pad1, 4, 11)
                        l.Draw("same")

                        pad2.cd()
                        ratioUp = up.Clone()
                        ratioUp.Divide(nomNoErr)   
                        ratioUp.SetTitle("")
                        ratioUp.GetYaxis().SetTitle("ratio to nominal")

                        ratioUp.GetXaxis().SetTitle(obsTitle[obs[4:]] + " [GeV]")
                        ratioUp.GetXaxis().SetTitleSize(0.1)
                        ratioUp.GetXaxis().SetTitleOffset(1.1)
                        ratioUp.GetYaxis().SetTitleSize(0.073)
                        ratioUp.GetYaxis().SetTitleOffset(0.5)

                        if down is not None:
                            ratioDown = down.Clone()
                            ratioDown.Divide(nomNoErr)
                            minRatioY = min(ratioUp.GetMinimum(), ratioDown.GetMinimum())
                            maxRatioY = max(ratioUp.GetMaximum(), ratioDown.GetMaximum())
                            padding = 0.05*(maxRatioY - minRatioY)
                            ratioUp.GetYaxis().SetRangeUser(minRatioY - padding, maxRatioY + padding)
                        
                        
                        ratioUp.Draw("hist e9")

                        if down is not None:
                            ratioDown.Draw("hist e9 same")
                        
                        
                        #line = TLine(ratioUp.GetXaxis().GetBinLowEdge(1), 1.0, ratioUp.GetXaxis().GetBinUpEdge(ratioUp.GetNbinsX()), 1.0)
                        #line.SetLineWidth(2)
                        #line.Draw("same")
                        nomRatio.Draw("hist e9 same")
                        for fmt in args.format:
                            canvasRatio.SaveAs("%s/%s/histplots/%s_%s_hist_%s_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,syst,fmt))
                        
                    if s == "tt":                    
                        # Draw ME scale systematics on one plot:
                        canvasRatio.cd()
                        canvasRatio.ResetDrawn()
                        canvasRatio.Draw()
                        canvasRatio.cd()

                        pad1.Draw()
                        pad2.Draw()

                        pad1.cd()
                            
                        l = TLegend(0.72,0.65, 0.88,0.88)
                        l.SetBorderSize(0)
                        l.AddEntry(nom, "nominal")
                        minY = nom.GetMinimum()
                        maxY = nom.GetMaximum()

                        MEscaleHists = {}
                        for syst in MEscaleSysts:
                            if syst not in systematics: continue  # Skip this one
                            MEscaleHists[syst] = h[s][lvl][obs][syst]["Up"][1725].Clone()
                            MEscaleHists[syst].Rebin(args.rebin)
                            MEscaleHists[syst].SetLineColor(MEscaleColor[syst])
                            MEscaleHists[syst].SetLineWidth(2)
                            MEscaleHists[syst].SetMarkerStyle(20)
                            MEscaleHists[syst].SetMarkerSize(0)

                            minY = min(minY, MEscaleHists[syst].GetMinimum())
                            maxY = max(maxY, MEscaleHists[syst].GetMaximum())
                            l.AddEntry(MEscaleHists[syst], syst)
                       
                        paddingY = 0.1 * (maxY - minY)
                        nom.SetTitle("%s  %s  ME scales" % (s,obsTitle[obs[4:]]) )
                        nom.Draw("hist")
                        nom.GetYaxis().SetRangeUser(max(minY - paddingY,0), maxY + paddingY)
                        for syst in MEscaleSysts:
                            if syst not in systematics: continue  # Skip this one
                            MEscaleHists[syst].Draw("hist same")

                        ###CMS_lumi.CMS_lumi(pad1, 4, 12)
                        CMS_lumi.CMS_lumi(pad1, 4, 11)
                        l.Draw("same")

                        pad2.cd()

                        line = TLine(nom.GetXaxis().GetBinLowEdge(1), 1.0, nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()), 1.0)
                        line.SetLineWidth(2)
                       
                        minRatioY = 0.95
                        maxRatioY = 1.05
                        
                        ratioMEscale = []
                        for syst in MEscaleSysts:
                            if syst not in systematics: continue  # Skip this one
                            ratioMEscale.append(MEscaleHists[syst].Clone("ratio_%s" % syst))
                            ratioMEscale[-1].Divide(nom)
                            minRatioY = min(minRatioY, ratioMEscale[-1].GetMinimum())
                            maxRatioY = max(maxRatioY, ratioMEscale[-1].GetMaximum())

                        paddingRatioY = 0.1 * (maxRatioY - minRatioY)
                        for i,_h in enumerate(ratioMEscale):
                            if i == 0:
                                _h.SetTitle("")
                                _h.Draw("hist")
                                _h.GetYaxis().SetRangeUser(minRatioY - paddingRatioY, maxRatioY + paddingRatioY)
                            else:
                                _h.Draw("hist same")

                        line.Draw("same")
                        for fmt in args.format:
                            canvasRatio.SaveAs("%s/%s/histplots/MEscales_%s_%s_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,fmt))

        for lvl in ['actual','morph']:
            # Graphs
            N = len(masses[s][lvl])
            mtmasses = [float(m)/10. for m in masses[s][lvl]]
            for obs in observables:
                if not args.nomoments:
                    # Moment graphs
                    g[s][lvl][obs] = { "nominal":{} }
                    g[s][lvl][obs]["nominal"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m1err"]))
                    g[s][lvl][obs]["nominal"]["m1"].SetName("%s_%s_%s_nominal_m1"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m2err"]))
                    g[s][lvl][obs]["nominal"]["m2"].SetName("%s_%s_%s_nominal_m2"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m3err"]))
                    g[s][lvl][obs]["nominal"]["m3"].SetName("%s_%s_%s_nominal_m3"%(s,lvl,obs))
                    g[s][lvl][obs]["nominal"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs]["nominal"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs]["nominal"]["m4err"]))
                    g[s][lvl][obs]["nominal"]["m4"].SetName("%s_%s_%s_nominal_m4"%(s,lvl,obs))

                # Rate graphs
                rateG[s][lvl][obs] = { "nominal":{} }
                rateG[s][lvl][obs]["nominal"] = TGraphErrors(N, array('d',mtmasses), array('d',rates[s][lvl][obs]["nominal"]["rate"]), array('d', [0.]*N), array('d', rates[s][lvl][obs]["nominal"]["rateError"]))
                rateG[s][lvl][obs]["nominal"].SetName("rate_%s_%s_%s_nominal" % (s,lvl,obs))
                rateG[s][lvl][obs]["nominal"].SetFillColor(0)
                
                if lvl == "actual":
                    rateG[s]["actual"][obs]["nominal"].SetMarkerColor(kBlue)
                    rateG[s]["actual"][obs]["nominal"].SetLineColor(kBlue)


                for syst in systematics:
                    if s == "tt" and syst in tWOnlySysts: continue
                    if s == "tW" and syst in ttOnlySysts: continue
                    if not args.nomoments:
                        # Syst moment graphs
                        g[s][lvl][obs][syst] = {"Up":{}}
                        g[s][lvl][obs][syst]["Up"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m1err"]))
                        g[s][lvl][obs][syst]["Up"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m2err"]))
                        g[s][lvl][obs][syst]["Up"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m3err"]))
                        g[s][lvl][obs][syst]["Up"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Up"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Up"]["m4err"]))
                    
                    # Syst rate graphs
                    rateG[s][lvl][obs][syst] = {}
                    rateG[s][lvl][obs][syst]["Up"] = TGraphErrors(N, array('d',mtmasses), array('d',rates[s][lvl][obs][syst]["Up"]["rate"]), array('d', [0.]*N), array('d', rates[s][lvl][obs][syst]["Up"]["rateError"]))
                    rateG[s][lvl][obs][syst]["Up"].SetName("rate_%s_%s_%s_%s%s" % (s,lvl,obs,syst,"Up" if syst not in oneSidedSysts else ""))
                    rateG[s][lvl][obs][syst]["Up"].SetFillColor(0)
                    if lvl == "actual":
                        rateG[s]["actual"][obs][syst]["Up"].SetMarkerColor(kBlue)
                        rateG[s]["actual"][obs][syst]["Up"].SetLineColor(kBlue)


                    if syst not in oneSidedSysts:
                        if not args.nomoments:
                            g[s][lvl][obs][syst]["Down"] = {}
                            g[s][lvl][obs][syst]["Down"]["m1"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m1"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m1err"]))
                            g[s][lvl][obs][syst]["Down"]["m2"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m2"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m2err"]))
                            g[s][lvl][obs][syst]["Down"]["m3"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m3"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m3err"]))
                            g[s][lvl][obs][syst]["Down"]["m4"] = TGraphErrors(N, array('d',mtmasses), array('d',moments[s][lvl][obs][syst]["Down"]["m4"]), array('d', [0.]*N), array('d', moments[s][lvl][obs][syst]["Down"]["m4err"]))

                        rateG[s][lvl][obs][syst]["Down"] = TGraphErrors(N, array('d',mtmasses), array('d',rates[s][lvl][obs][syst]["Down"]["rate"]), array('d', [0.]*N), array('d', rates[s][lvl][obs][syst]["Down"]["rateError"]))
                        rateG[s][lvl][obs][syst]["Down"].SetName("rate_%s_%s_%s_%sDown" % (s,lvl,obs,syst))
                        rateG[s][lvl][obs][syst]["Down"].SetFillColor(0)
                        if lvl == "actual":
                            rateG[s]["actual"][obs][syst]["Down"].SetMarkerColor(kBlue)
                            rateG[s]["actual"][obs][syst]["Down"].SetLineColor(kBlue)
    
    if not args.noplots:
        # Write to output root file
        f = TFile.Open("%s/plots.root"%args.outDir, "recreate")

        for s in signal:
            for lvl in ['actual','morph']:
                for obs in observables:
                    if not args.nomoments:
                        g[s][lvl][obs]["nominal"]["m1"].Write("%s_%s_%s_nominal_m1"%(s,lvl,obs))
                        g[s][lvl][obs]["nominal"]["m2"].Write("%s_%s_%s_nominal_m2"%(s,lvl,obs))
                        g[s][lvl][obs]["nominal"]["m3"].Write("%s_%s_%s_nominal_m3"%(s,lvl,obs))
                        g[s][lvl][obs]["nominal"]["m4"].Write("%s_%s_%s_nominal_m4"%(s,lvl,obs))
                   
                    rateG[s][lvl][obs]["nominal"].Write()                    

                    for syst in systematics:
                        if s == "tt" and syst in tWOnlySysts: continue
                        if s == "tW" and syst in ttOnlySysts: continue
                        if not args.nomoments:
                            g[s][lvl][obs][syst]["Up"]["m1"].Write("%s_%s_%s_%sUp_m1"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Up"]["m2"].Write("%s_%s_%s_%sUp_m2"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Up"]["m3"].Write("%s_%s_%s_%sUp_m3"%(s,lvl,obs,syst))
                            g[s][lvl][obs][syst]["Up"]["m4"].Write("%s_%s_%s_%sUp_m4"%(s,lvl,obs,syst))
                        rateG[s][lvl][obs][syst]["Up"].Write()
                        if syst not in oneSidedSysts:
                            if not args.nomoments:
                                g[s][lvl][obs][syst]["Down"]["m1"].Write("%s_%s_%s_%sDown_m1"%(s,lvl,obs,syst))
                                g[s][lvl][obs][syst]["Down"]["m2"].Write("%s_%s_%s_%sDown_m2"%(s,lvl,obs,syst))
                                g[s][lvl][obs][syst]["Down"]["m3"].Write("%s_%s_%s_%sDown_m3"%(s,lvl,obs,syst))
                                g[s][lvl][obs][syst]["Down"]["m4"].Write("%s_%s_%s_%sDown_m4"%(s,lvl,obs,syst))
                            rateG[s][lvl][obs][syst]["Down"].Write()

    if not args.nomoments:
        with gzip.open("%s/moments.pklz"%args.outDir, "wb") as mf:
            pickle.dump(moments, mf, protocol=pickle.HIGHEST_PROTOCOL)


    c.cd()

    if not args.nomoments:# and not args.mtscan:
        if not args.noplots: print "About to draw moment plots"
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
                        if useMomentCalibration and not args.noplots:
                            # Draw calibration curve
#                        reclvl,obsname = obs.split("_")
                            reclvl = obs[:3]
                            obsname = obs[4:]
                            calib[s][lvl][obs][moment].calibG.SetTitle("%s %s  %s  Moment %s  Calibration curve" % (reclvl, obsTitle[obsname], lvl, moment[1]))
                            calib[s][lvl][obs][moment].calibG.Draw("ap")
                            txt.DrawLatex(0.6, 0.35, "Slope")
                            txt.DrawLatex(0.75, 0.35, "Offset")
                            txt.DrawLatex(0.57, 0.3, "%.4f #pm %.4f" % (calib[s][lvl][obs][moment].slope, calib[s][lvl][obs][moment].slopeError))
                            txt.DrawLatex(0.72, 0.3, "%.4f #pm %.4f" % (calib[s][lvl][obs][moment].offset, calib[s][lvl][obs][moment].offsetError))

                            for fmt in args.format: 
                                c.SaveAs("%s/%s/calibration/calibration_%s_%s_%s_%s.%s" % (args.outDir,obs,s,lvl,obs,moment,fmt))
                            
                        # Get nominal fit line
                        g[s][lvl][obs]["nominal"][moment].Fit("pol1", "Q")
                        fitLines[s][lvl][obs]["nominal"][moment] = g[s][lvl][obs]["nominal"][moment].GetFunction("pol1") 
                        fitLines[s][lvl][obs]["nominal"][moment].SetLineColor(kBlack)    
                        g[s][lvl][obs]["nominal"][moment].SetLineColor(kBlack)
                        g[s][lvl][obs]["nominal"][moment].SetLineWidth(2)
                        g[s][lvl][obs]["nominal"][moment].SetMarkerStyle(20)
                        g[s][lvl][obs]["nominal"][moment].SetMarkerSize(1)
                        g[s][lvl][obs]["nominal"][moment].SetMarkerColor(kBlack)
                        # Fit line 
                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetLineColor(kBlack)
                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetLineWidth(2)
                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetMarkerStyle(20)
                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetMarkerSize(1)
                        g[s][lvl][obs]["nominal"][moment].GetFunction("pol1").SetMarkerColor(kBlack)
                    

                    for syst in systematics:
                        #if lvl == "actual" and syst in separateSystSamples: continue
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
                            g[s][lvl][obs][syst]["Up"][moment].SetLineWidth(2)
                            g[s][lvl][obs][syst]["Up"][moment].SetMarkerStyle(22)
                            g[s][lvl][obs][syst]["Up"][moment].SetMarkerSize(2)
                            g[s][lvl][obs][syst]["Up"][moment].SetMarkerColor(kRed)
                            # Fit line
                            g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetLineColor(kRed)
                            g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetLineWidth(2)
                            g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetMarkerStyle(22)
                            g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetMarkerSize(2)
                            g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1").SetMarkerColor(kRed)

#                        upFit = g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1")
                            fitLines[s][lvl][obs][syst]["Up"][moment] = g[s][lvl][obs][syst]["Up"][moment].GetFunction("pol1")
                            
                            
                            #l.AddEntry(g[s][lvl][obs][syst]["Up"][moment], "%s%s" % (syst, " Up" if syst not in oneSidedSysts else ""), "l")
                            l.AddEntry(fitLines[s][lvl][obs][syst]["Up"][moment], "%s%s" % (syst, " Up" if syst not in oneSidedSysts else ""))
                            l.AddEntry(g[s][lvl][obs]["nominal"][moment], "nominal")
                            mg[s][lvl][obs][syst][moment] = TMultiGraph()
                            mg[s][lvl][obs][syst][moment].SetName("%s_%s_%s_%s_%s"%(s,lvl,obs,syst,moment))

                            mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs]["nominal"][moment].Clone(), "p")
                            mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs][syst]["Up"][moment].Clone(), "p")
                            if syst not in oneSidedSysts:
                                g[s][lvl][obs][syst]["Down"][moment].Fit("pol1","Q")
                                g[s][lvl][obs][syst]["Down"][moment].SetLineColor(kBlue)
                                g[s][lvl][obs][syst]["Down"][moment].SetLineWidth(2)
                                g[s][lvl][obs][syst]["Down"][moment].SetMarkerStyle(23)
                                g[s][lvl][obs][syst]["Down"][moment].SetMarkerSize(2)
                                g[s][lvl][obs][syst]["Down"][moment].SetMarkerColor(kBlue)
                                # For legend
                                g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetLineColor(kBlue)
                                g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetLineWidth(2)
                                g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetMarkerStyle(23)
                                g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetMarkerSize(2)
                                g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1").SetMarkerColor(kBlue)
                                fitLines[s][lvl][obs][syst]["Down"][moment] = g[s][lvl][obs][syst]["Down"][moment].GetFunction("pol1")
                                mg[s][lvl][obs][syst][moment].Add(g[s][lvl][obs][syst]["Down"][moment].Clone(), "p")
                                #l.AddEntry(g[s][lvl][obs][syst]["Down"][moment], "%s Down" % syst, "l")
                                l.AddEntry(fitLines[s][lvl][obs][syst]["Down"][moment], "%s Down" % syst)

                            
                            if not args.noplots:
                                mg[s][lvl][obs][syst][moment].SetTitle("%s  %s %s  %s  %s%s  Moment %s" % (signalTitle[s],obs[:3],obsTitle[obs[4:]],syst,"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed",moment[1]))
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
                                
                                for fmt in args.format:
                                    c.SaveAs("%s/%s/%s_%s_%s_%s_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,syst,moment,fmt))

                for obs in observables:
                    for moment in ["m1","m2","m3","m4"]:

                        # Gen/Rec plots
                        if obs[:4] == "rec_" and obs.replace("rec_","gen_") in observables:

                            mgobs = "recgen_%s" % obs[4:]
                            mg[s][lvl][mgobs] = {}
                            genobs = obs.replace("rec_","gen_")
                            l = TLegend(0.12,0.7, 0.3,0.88)
                            l.SetBorderSize(0)
                            fitLines[s][lvl][obs]["nominal"][moment].SetLineColor(kBlack)
                            fitLines[s][lvl][genobs]["nominal"][moment].SetLineColor(kRed) 
                            
                            g[s][lvl][obs]["nominal"][moment].SetLineColor(kBlack)
                            g[s][lvl][genobs]["nominal"][moment].SetLineColor(kRed)
                            
                            l.AddEntry(g[s][lvl][obs]["nominal"][moment], "rec")
                            l.AddEntry(g[s][lvl][genobs]["nominal"][moment], "gen")

                            mg[s][lvl][mgobs][moment] = TMultiGraph()
                            mg[s][lvl][mgobs][moment].SetName("%s_%s_%s_%s"%(s,lvl,mgobs,moment))

                            mg[s][lvl][mgobs][moment].Add(g[s][lvl][obs]["nominal"][moment].Clone(), "p")
                            mg[s][lvl][mgobs][moment].Add(g[s][lvl][genobs]["nominal"][moment].Clone(), "p")
                            if not args.noplots:
                                mg[s][lvl][mgobs][moment].SetTitle("%s  %s  %s%s  Moment %s" % (signalTitle[s],obsTitle[obs[4:]],"" if args.interp =="" else args.interp+"  ","actual" if lvl == "actual" else "morphed",moment[1]))
                                mg[s][lvl][mgobs][moment].Draw("a9")

                                mg[s][lvl][mgobs][moment].GetXaxis().SetTitle("m_{t} [GeV]")
                                mg[s][lvl][mgobs][moment].GetYaxis().SetTitle("%s [GeV]%s" % (moment, "" if moment == "m1" else "^{%s}"%moment[1]))
                                mg[s][lvl][mgobs][moment].GetYaxis().SetTitleOffset(1.2)

                                l.Draw("same")
                                txt.DrawLatex(0.44, 0.35, "Moment %s Slopes"%moment[1])
                                txt.DrawLatex(0.407, 0.25,  "Rec: %.4f #pm %.4f" % (fitLines[s][lvl][obs]["nominal"][moment].GetParameter(1), fitLines[s][lvl][obs]["nominal"][moment].GetParError(1)))
                                txt.DrawLatex(0.4, 0.2, "#color[2]{Gen: %.4f #pm %.4f}" % (fitLines[s][lvl][genobs]["nominal"][moment].GetParameter(1), fitLines[s][lvl][genobs]["nominal"][moment].GetParError(1)))
                                
                                txt.DrawLatex(0.76, 0.35, "Offsets")
                                txt.DrawLatex(0.707, 0.25,  "Rec: %.4f #pm %.4f" % (fitLines[s][lvl][obs]["nominal"][moment].GetParameter(0), fitLines[s][lvl][obs]["nominal"][moment].GetParError(0)))
                                txt.DrawLatex(0.7, 0.2, "#color[2]{Gen: %.4f #pm %.4f}" % (fitLines[s][lvl][genobs]["nominal"][moment].GetParameter(0), fitLines[s][lvl][genobs]["nominal"][moment].GetParError(0)))
                                
                                for fmt in args.format:
                                    c.SaveAs("%s/%s/recgen_%s_%s_%s_%s.%s" % (args.outDir,obs,signalName[s],lvl,obs,moment,fmt))

    
    # Compare morphed/actual moments and rates
    momentMG = {}
    rateMG = {}
    for s in signal:
        if args.mtscan: break
        momentMG[s] = {}
        rateMG[s] = {}
        for obs in observables:
            momentMG[s][obs] = {"nominal":TMultiGraph()}
            momentMG[s][obs]["nominal"].SetName("comp_%s_%s_moments")
            #momentMG[s][obs]["nomainl"].Add(
            
            #l = TLegend(0.75, 0.3, 0.9, 0.45)

            rateMG[s][obs] = {"nominal":TMultiGraph()}
            rateMG[s][obs]["nominal"].SetName("comp_%s_%s_rates" % (s,obs))
            rateMG[s][obs]["nominal"].Add(rateG[s]["morph"][obs]["nominal"])
            rateMG[s][obs]["nominal"].Add(rateG[s]["actual"][obs]["nominal"])

            l = TLegend(0.7, 0.3, 0.85, 0.45)
            l.SetBorderSize(0)
            l.AddEntry(rateG[s]["actual"][obs]["nominal"], "actual")
            l.AddEntry(rateG[s]["morph"][obs]["nominal"], "morph")

            if not args.noplots:
                rateMG[s][obs]["nominal"].Draw("alp9")
                rateMG[s][obs]["nominal"].SetTitle("Nominal Rates")
                rateMG[s][obs]["nominal"].GetXaxis().SetTitle("m_{t} [GeV]")
                rateMG[s][obs]["nominal"].GetYaxis().SetTitle("Rate")

                l.Draw("same")
                for fmt in args.format:
                    c.SaveAs("%s/%s/rates/rate_%s_%s_nominal.%s" % (args.outDir, obs,signalName[s],obs,fmt))

            for syst in systematics:
                if s == "tt" and syst in tWOnlySysts: continue
                if s == "tW" and syst in ttOnlySysts: continue
                rateMG[s][obs][syst] = {"Up":TMultiGraph()}
                rateMG[s][obs][syst]["Up"].SetName("comp_%s_%s_%s%s_rates" % (s,obs,syst,"Up" if syst not in oneSidedSysts else ""))
                rateMG[s][obs][syst]["Up"].Add(rateG[s]["morph"][obs][syst]["Up"])
                rateMG[s][obs][syst]["Up"].Add(rateG[s]["actual"][obs][syst]["Up"])

                l = TLegend(0.7, 0.3, 0.85, 0.45)
                l.SetBorderSize(0)
                l.AddEntry(rateG[s]["actual"][obs][syst]["Up"], "actual")
                l.AddEntry(rateG[s]["morph"][obs][syst]["Up"], "morph")

                if not args.noplots:
                    rateMG[s][obs][syst]["Up"].Draw("alp9")
                    rateMG[s][obs][syst]["Up"].SetTitle("%s%s Rates" % (syst, " Up" if syst not in oneSidedSysts else ""))
                    rateMG[s][obs][syst]["Up"].GetXaxis().SetTitle("m_{t} [GeV]")
                    rateMG[s][obs][syst]["Up"].GetYaxis().SetTitle("Rate")

                    l.Draw("same")
                    for fmt in args.format:
                        c.SaveAs("%s/%s/rates/rate_%s_%s_%s%s.%s" % (args.outDir, obs, signalName[s],obs,syst," Up" if (syst not in oneSidedSysts and syst != "nominal") else "",fmt))
                
                if syst not in oneSidedSysts and syst != "nominal":
                    rateMG[s][obs][syst]["Down"] = TMultiGraph()
                    rateMG[s][obs][syst]["Down"].SetName("comp_%s_%s_%sDown_rates" % (s,obs,syst))
                    rateMG[s][obs][syst]["Down"].Add(rateG[s]["morph"][obs][syst]["Down"])
                    rateMG[s][obs][syst]["Down"].Add(rateG[s]["actual"][obs][syst]["Down"])

                    l = TLegend(0.7, 0.3, 0.85, 0.45)
                    l.SetBorderSize(0)
                    l.AddEntry(rateG[s]["actual"][obs][syst]["Down"], "actual")
                    l.AddEntry(rateG[s]["morph"][obs][syst]["Down"], "morph")

                    if not args.noplots:
                        rateMG[s][obs][syst]["Down"].Draw("alp9")
                        rateMG[s][obs][syst]["Down"].SetTitle("%s Down Rates" % (syst))
                        rateMG[s][obs][syst]["Down"].GetXaxis().SetTitle("m_{t} [GeV]")
                        rateMG[s][obs][syst]["Down"].GetYaxis().SetTitle("Rate")

                        l.Draw("same")
                        for fmt in args.format:
                            c.SaveAs("%s/%s/rates/rate_%s_%s_%sDown.%s" % (args.outDir, obs, signalName[s],obs,syst,fmt))



    if not args.noplots:
        f.Close()


    if not args.nomoments:# and not args.mtscan:
        # Extracted masses
        uncalibrated_mt = {}
        calibrated_mt = {}
        for s in signal: 
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
                        if useMomentCalibration: 
                            calibrated_mt[s][lvl][obs]["nominal"][moment] = {"mass":{}, "error":{}}
                    
                        # 172.5 value is midpoint of list
                        nomEntry = len(moments[s][lvl][obs]["nominal"][moment]) // 2

                        # Nominal moment and statistical uncertainty
                        if useDataObs:
                            nom_mt, nom_mterr = massFromMoment(fitLines[s][lvl][obs]["nominal"][moment], dataobs_moments[obs][moment], dataobs_moments[obs][moment+"err"])
                        else:
                            nom_mt, nom_mterr = massFromMoment(fitLines[s][lvl][obs]["nominal"][moment], moments[s][lvl][obs]["nominal"][moment][nomEntry], moments[s][lvl][obs]["nominal"][moment+"err"][nomEntry])
                        uncalibrated_mt[s][lvl][obs]["nominal"][moment]["mass"] = nom_mt
                        uncalibrated_mt[s][lvl][obs]["nominal"][moment]["error"] = nom_mterr

                        if useMomentCalibration: 
                            calibrated_nom_mt = calib[s][lvl][obs][moment].calibrate(nom_mt)
                            calibrated_mt[s][lvl][obs]["nominal"][moment]["mass"] = calibrated_nom_mt
                            calibrated_mt[s][lvl][obs]["nominal"][moment]["error"] = calib[s][lvl][obs][moment].calibrateError(nom_mterr)


                        for syst in systematics:
                            #if lvl == "actual" and syst in separateSystSamples: continue
                            if s == "tt" and syst in tWOnlySysts: continue
                            if s == "tW" and syst in ttOnlySysts: continue
                            if syst not in uncalibrated_mt[s][lvl][obs]:
                                uncalibrated_mt[s][lvl][obs][syst] = {"Up":{}}
                                if useMomentCalibration: 
                                    calibrated_mt[s][lvl][obs][syst] = {"Up":{}}
                                #uncalibrated_mt[s][lvl][obs][syst] = {"Up":{moment:{}}}
                                #calibrated_mt[s][lvl][obs][syst] = {"Up":{moment:{}}}
                                
                                if sys not in oneSidedSysts:
                                    uncalibrated_mt[s][lvl][obs][syst]["Down"] = {}
                                    if useMomentCalibration: 
                                        calibrated_mt[s][lvl][obs][syst]["Down"] = {}
                                    #uncalibrated_mt[s][lvl][obs][syst]["Down"] = {moment:{}}
                                    #calibrated_mt[s][lvl][obs][syst]["Down"] = {moment:{}}

                            uncalibrated_mt[s][lvl][obs][syst]["Up"][moment] = {} 
                            if useMomentCalibration: 
                                calibrated_mt[s][lvl][obs][syst]["Up"][moment] = {} 
                            if sys not in oneSidedSysts:
                                uncalibrated_mt[s][lvl][obs][syst]["Down"][moment] = {} 
                                if useMomentCalibration: 
                                    calibrated_mt[s][lvl][obs][syst]["Down"][moment] = {} 

                           
                            nomEntry = len(moments[s][lvl][obs][syst]["Up"][moment]) // 2
#                        print "Now on %s %s %s" % (lvl,syst,moment)
                            
                            try:
                                #mtUp = massFromMoment(fitLines[s][lvl][obs][syst]["Up"][moment], moments[s][lvl][obs][syst]["Up"][moment][nomEntry] )
                                #mtUp = massFromMoment(fitLines[s][lvl][obs]['nominal'][moment], moments[s][lvl][obs][syst]['Up'][moment][nomEntry] )
                                # Fit nominal moment to systematic fit line
                                if useDataObs:
                                    mtUp = massFromMoment(fitLines[s][lvl][obs][syst]["Up"][moment], dataobs_moments[obs][moment] )
                                else:
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

                            if useMomentCalibration: 
                                calibrated_mtUp = calib[s][lvl][obs][moment].calibrate(mtUp)
                                calibrated_mt[s][lvl][obs][syst]["Up"][moment]["mass"] = calibrated_mtUp
                                calibrated_mt[s][lvl][obs][syst]["Up"][moment]["impact"] = calibrated_mtUp - calibrated_nom_mt
                            
                            if syst not in oneSidedSysts:
                                #mtDown = massFromMoment(fitLines[s][lvl][obs][syst]["Down"][moment], moments[s][lvl][obs][syst]["Down"][moment][nomEntry])
                                
                                if useDataObs:
                                    mtDown = massFromMoment(fitLines[s][lvl][obs][syst]["Down"][moment], dataobs_moments[obs][moment])
                                else:
                                    mtDown = massFromMoment(fitLines[s][lvl][obs][syst]["Down"][moment], moments[s][lvl][obs]['nominal'][moment][nomEntry])
                                #mtDown = massFromMoment(fitLines[s][lvl][obs]['nominal'][moment], moments[s][lvl][obs][syst]['Down'][moment][nomEntry] )
                                uncalibrated_mt[s][lvl][obs][syst]["Down"][moment]["mass"] = mtDown
                                uncalibrated_mt[s][lvl][obs][syst]["Down"][moment]["impact"] = nom_mt - mtDown

                                if useMomentCalibration:
                                    calibrated_mtDown = calib[s][lvl][obs][moment].calibrate(mtDown)
                                    calibrated_mt[s][lvl][obs][syst]["Down"][moment]["mass"] = calibrated_mtDown
                                    calibrated_mt[s][lvl][obs][syst]["Down"][moment]["impact"] = calibrated_nom_mt - calibrated_mtDown
                            

        # Write moment json file
        for s in signal: 
            for lvl in ["actual", "morph"]:
                for moment in ["m1","m2","m3","m4"]:
                    for obs in observables:
                        for calib in ["uncalibrated", "calibrated"]:
                            #if not useMomentCalibration: 
                            if calib == "uncalibrated": 
                                nominal_mt = uncalibrated_mt[s][lvl][obs]['nominal'][moment]['mass']
                                stat_err = uncalibrated_mt[s][lvl][obs]['nominal'][moment]['error']
                            else:
                                nominal_mt = calibrated_mt[s][lvl][obs]['nominal'][moment]['mass']
                                stat_err = calibrated_mt[s][lvl][obs]['nominal'][moment]['error']

                            nominal_mt *= 10.
                            stat_err *= 10.

#                        jsondata = {"POIs":[ {"fit":[nominal_mt - stat_err , nominal_mt, nominal_mt + stat_err], "name":"MT", "stat":[nominal_mt - stat_err , nominal_mt, nominal_mt + stat_err]} ] }
#                        jsondata["params"] = []
                            jsondata = {"POIs":[ {"fit":[nominal_mt - stat_err , nominal_mt, nominal_mt + stat_err], "name":"MT", }] }
                            statVals = \
                            {
                                u'name': u'stat',
                                u'MT': [nominal_mt - stat_err, nominal_mt, nominal_mt + stat_err],
                                u'impact_MT': stat_err,
                                u'impact_r': 0.0,
                                u'prefit': [-1.0, 0.0, 1.0],
                                u'fit': [ -1.0, 0.0, 1.0],
                                u'groups': [],
                                u'r': [1.0, 1.0, 1.0],
                                u'type': "Gaussian",
                            } 
                            jsondata["params"] = [statVals]
                            for syst in systematics:
                                if s == "tt" and syst in tWOnlySysts: continue
                                if s == "tW" and syst in ttOnlySysts: continue
                                if calib == "uncalibrated": 
                                    mtUp = uncalibrated_mt[s]['morph'][obs][syst]["Up"][moment]["mass"]
                                else:
                                    mtUp = calibrated_mt[s]['morph'][obs][syst]["Up"][moment]["mass"]
                                if syst not in oneSidedSysts:
                                    if not useMomentCalibration: 
                                        mtDown = uncalibrated_mt[s]['morph'][obs][syst]["Down"][moment]["mass"]
                                    else:
                                        mtDown = calibrated_mt[s]['morph'][obs][syst]["Down"][moment]["mass"]
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

                            outF = "%s/%s/%s_%s_%s_%s_%s_%s" % (args.outDir, obs, s, obs, lvl, moment,calib,args.json)
                            with open(outF, "w") as f:
                                json.dump(jsondata, f, indent=4, separators=(',', ': '))

                            print "Impact data on moment %s written to %s" % (moment,outF)

if __name__ == "__main__":
    sys.exit(main())



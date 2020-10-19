#!/usr/bin/env python
from ROOT import *
import os
import sys 
from argparse import ArgumentParser
from array import array
import numpy as np
import math
from copy import deepcopy

gStyle.SetOptStat(0)
gROOT.SetBatch(True)

obsUnfolding= {\
    "ptneg"   : "pt_neg",
    "ptpos"   : "pt_pos",
    "ptll"    : "pt_ll",
    "Mll"     : "m_ll",
    "Ep_Em"   : "Ep_Em",
    "Epos"    : "E_pos",
    "Eneg"    : "E_neg",
    "ptp_ptm" : "ptp_ptm",
}

obsTitle = {\
    "ptll"   : "p_{T}(ll)",
    "Mll"    : "M(ll)",
    "ptpos"  : "p_{T}(l^{+})",
    "ptneg"  : "p_{T}(l^{-})",
    "Epos"   : "E(l^{+})",
    "Eneg"   : "E(l^{-})",
    "ptp_ptm" : "p_{T}(l^{+}) + p_{T}(l^{-})",
    "Ep_Em"   : "E(l^{+}) + E(l^{-})"
}


def calcMoments(h):
    sums = {0:0., 1:0., 2:0., 4:0.}
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
    
    return m1,m1err,m2,m2err

def calcMomentFunc(h):
#    sums = {0:0., 1:0., 2:0., 3:0., 4:0., 6:0., 8:0.}
#    kvals = sorted(sums.keys())
    m1 = h.GetMean()
    m1err = h.GetMeanError()

    sigma = h.GetRMS()
    m2 = sigma**2 + m1**2
    m2err = 2 * sigma * h.GetRMSError() + 2 * m1 * m1err


    return m1,m1err,m2,m2err


def massFromMoment(fit, moment, momentError = None):
    offset = fit.GetParameter(0)
    slope = fit.GetParameter(1)

    if momentError is not None:
        return (moment - offset) / slope, momentError / slope
    else:
        return (moment - offset) / slope


parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mt1725_ttactualdataobs_ptll_mtTemplatesForCH.root")
parser.add_argument("-o", "--outDir", default="Manual_likelihoods")
parser.add_argument("-s", "--sig", default="tttW", choices=["tt","tW","tttw","tttW"], help="signal process(es)")
parser.add_argument("-nm", "--numMT", type=int, default=7, choices=[3,5,7], help="number of tt mass points to use")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reco level")
parser.add_argument("--obs", default="ptll", help="observable")
#parser.add_argument("-d", "--dataobs", "--data_obs", dest="useDataObs", action="store_true", default=False, help="Use dataobs for observed distribution instead of asimiov set")
parser.add_argument("--rebin", type=int, default=1, help="rebin quantity")
parser.add_argument("--scale", action="store_true", default=False, help="apply bin uncertainty scaling from scaling file")
#parser.add_argument("--scaleF", default="/uscms_data/d3/msaunder/combine/CMSSW_10_2_13/src/UserCode/unfolding/chi2_unfolded_tttW_allobs/scaling.root", help="bin uncertainty scaling file")
parser.add_argument("--scaleF", default="/uscms_data/d3/msaunder/combine/CMSSW_10_2_13/src/UserCode/unfolding/chi2_unfolded_tt_allobs/scaling.root", help="bin uncertainty scaling file")
parser.add_argument("-a", "--asimov", action="store_true", default=False, help="use asimov set instead of data_obs")
parser.add_argument("--toy", "--toyOpt", dest="toyOpt", default="", choices=["mc","MC","Mc","MCStat","pd","PD","Pd","Pseudodata"], help="create toy templates by fluctuating bin contents according to MC uncertainty")
parser.add_argument("--gaussErrs", action="store_true", default=False, help="use gaussian errors instead of poisson")
parser.add_argument("-g", "--gauss", "--gaussML", dest="useGaussML", action="store_true", default=False, help="use gaussian ML estimate (for when errors are not poisson")
parser.add_argument("-n", "-N", "--toyN", dest="toyN", type=int, default=10, help="number of toys to create")
parser.add_argument("--saveToys", action="store_true", default=False, help="save toys in root file")
parser.add_argument("--toySeed", type=int, default=0, help="seed for TRandom3 to calculate toy fluctuations")
parser.add_argument("--morph", action="store_true", default=False, help="use morphed templates instead of actual templates")
parser.add_argument("--norm", action="store_true", default=False, help="normalize to dataobs rate")
parser.add_argument("--floatToyRate", action="store_true", default=False, help="float toy normalization rate")

args = parser.parse_args()

saveAllToys = args.saveToys
normToyRate = not args.floatToyRate
# Number of tt mass points to use
numMassPoints = args.numMT

inF = args.inF
outDir = args.outDir
reco = args.reco
obs = args.obs
recoObs = "%s_%s" % (reco,obs)
rebin = args.rebin

if rebin < 1:
    print "Invalid rebin: %d", rebin
    sys.exit()

# http://cp3.irmp.ucl.ac.be/~delaere/BND2010/AnalysisMethods_Delaere_BND2010.pdf#page=105
useGaussML = args.useGaussML
if useGaussML:
    print "Using gaussian MLE for nll"
else:
    print "Using standard poisson MLE for nll"

signal = args.sig
if signal == "tttw": signal = "tttW"
print "Using signal: %s" % ("tt+tW" if signal == "tttW" else signal)
useUncScaling = args.scale
uncScalingFile = args.scaleF

usePoissonErrors = not args.gaussErrs

if useUncScaling:
    print "Scaling bin uncertainties %susing file: %s" % ("(assumed to be Poisson) " if usePoissonErrors else "", uncScalingFile)
    scalingF = TFile.Open(uncScalingFile)
    if usePoissonErrors:
        scaling = scalingF.Get("%s_poisson_scaling" % obsUnfolding[obs])
    else:
        scaling = scalingF.Get("%s_scaling" % obsUnfolding[obs])
    
    scaling.SetDirectory(0)
    scalingF.Close()


unfF = TFile.Open(uncScalingFile.replace("scaling.root", "unfolded.root"))
unfH = unfF.Get("%s/%s_mt1725_unfolded" % (obsUnfolding[obs], obsUnfolding[obs]))
unfH.SetDirectory(0)

poisson_unfH = unfF.Get("%s/%s_mt1725_unfolded_poisson" % (obsUnfolding[obs], obsUnfolding[obs]))
poisson_unfH.SetDirectory(0)

unfF.Close()


useDataObs = not args.asimov
useMorphed = args.morph
normalizeExpected = args.norm
if usePoissonErrors:
    print "Using poisson errors for smearing"

toyOpt = args.toyOpt
if toyOpt in ["MC","mc","Mc"]: toyOpt = "MCStat"
if toyOpt in ["PD","pd","Pd"]: toyOpt = "Pseudodata"
    
ntoys = args.toyN if toyOpt != "" else 0 
toySeed = args.toySeed
rnd = TRandom3(toySeed)


if toyOpt != "":
    if ntoys < 1:
        print "Invalid number of toys: %d" % ntoys
        sys.exit()

    print "Creating %d %s toys with %s normalization and seed: %d" % (ntoys, toyOpt, "fixed" if normToyRate else "floating", toySeed)

backgrounds = ["DY", "ST_bkgd", "TTV", "Diboson"]


bkgd = None

if signal == "tt":
    if numMassPoints == 7:
        actualMasses = [166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5]
    elif numMassPoints == 5:
        actualMasses = [169.5, 171.5, 172.5, 173.5, 175.5]
    elif numMassPoints == 3:
        actualMasses = [169.5, 172.5, 175.5]
else:
    actualMasses = [169.5, 172.5, 175.5]



if rebin > 1:
    print "Rebin: %d -> 1" % rebin
sig = {}
actualSig = {}
if useMorphed:
    masses = [166.5 + 0.1*i for i in range(121)]
    print "Using morphed %s templates:"%signal, masses
else:
    masses = actualMasses
    print "Using actual %s templates for masses:"%signal, masses

f = TFile.Open(inF)

for m in masses:
    print "Loading %s/%s%s%d" % (recoObs, signal, "" if useMorphed else "actual", int(m*10))
    sig[m] = f.Get("%s/%s%s%d" % (recoObs, signal, "" if useMorphed else "actual", int(m*10)))
    sig[m].SetDirectory(0)
    if rebin > 1:
        sig[m].Rebin(rebin)
    if m in actualMasses:
        actualSig[m] = f.Get("%s/%sactual%d" % (recoObs, signal, int(m*10)))
        actualSig[m].SetDirectory(0)
        if rebin > 1:
            actualSig[m].Rebin(rebin)

for b in backgrounds:
    if bkgd is None:
        bkgd = f.Get("%s/%s" % (recoObs,b)).Clone("%s_backgrounds" % recoObs)
        bkgd.SetDirectory(0)
    else:
        bkgd.Add(f.Get("%s/%s" % (recoObs,b)))

if rebin > 1:
    bkgd.Rebin(rebin)
    unfH.Rebin(rebin)
    poisson_unfH.Rebin(rebin)

if useDataObs:
    dataobs = f.Get("%s/data_obs" % recoObs)
    dataobs.SetDirectory(0)
    if rebin > 1:
        dataobs.Rebin(rebin)
f.Close()

exp_moments = {0:{}}
obs_moments = {0:{}}

expected = {0:{}}
observed = {0:{}}
for m in masses:
    expected[0][m] = sig[m].Clone("expected_mt%d" % (int(m*10)))
    expected[0][m].Add(bkgd)
    if normalizeExpected:
        expected[0][m].Scale((actualSig[172.5].Integral() + bkgd.Integral())/expected[0][m].Integral())

    m1,m1err,m2,m2err = calcMoments(expected[0][m])
    exp_moments[0][m] = {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}

    if toyOpt == "MCStat" and ntoys > 0:
        print "Creating MC stat toys for mt = %.1f" % m

        # Create MC stat toys here
        for t in range(1, ntoys+1):
            if t not in expected:
                expected[t] = {}
                exp_moments[t] = {}
            expected[t][m] = expected[0][m].Clone(expected[0][m].GetName()+"_toy_%d " % t)
            for b in range(1, expected[0][m].GetNbinsX()+1):
                content = expected[t][m].GetBinContent(b)
                error = expected[t][m].GetBinError(b)
                newContent = rnd.Gaus(content,error) if not usePoissonErrors else rnd.Gaus(content, content**0.5)
                newError = (newContent/content)*error
                
                expected[t][m].SetBinContent(b, newContent)
                expected[t][m].SetBinError(b, newError)

            if normToyRate:
                expected[t][m].Scale(expected[0][m].Integral() / expected[t][m].Integral())

            m1,m1err,m2,m2err = calcMoments(expected[t][m])
            exp_moments[t][m] = {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}

if useDataObs:
    observed[0][-1] = dataobs.Clone("%s_dataobs" % obs)
    if usePoissonErrors:
        old = poisson_unfH.Clone("old")
        for b in range(1, observed[0][-1].GetNbinsX()+1):
            observed[0][-1].SetBinError(b, observed[0][-1].GetBinContent(b)**0.5)
    else:
        old = unfH.Clone("old")

    old.Add(bkgd)
            

    #old = observed[0][-1].Clone("old")
    # Apply uncertainty scaling here (if applicable)
    if useUncScaling:
        # Scaling is for signal only
        # Subtract off background
        observed[0][-1].Add(bkgd, -1)

        # Scale bin errors
        for b in range(1, observed[0][-1].GetNbinsX()+1):
            observed[0][-1].SetBinError(b, scaling.GetBinContent(b) * (observed[0][-1].GetBinError(b)**2 - 2*bkgd.GetBinError(b)**2)**0.5)
            #observed[0][-1].SetBinError(b, (observed[0][-1].GetBinError(b)**2 - 2*bkgd.GetBinError(b)**2)**0.5)
            #observed[0][-1].SetBinError(b, observed[0][-1].GetBinError(b) * scaling.GetBinContent(b))
    
        # Add background back
        observed[0][-1].Add(bkgd)

    diff = old.Clone("diff")
    diffErr = old.Clone("diff")

    diff.Add(observed[0][-1], -1)
    for b in range(1, old.GetNbinsX()+1):
        diffErr.SetBinContent(b, old.GetBinError(b) - observed[0][-1].GetBinError(b))

    print "diff:", diff.Integral()
    print "diffErr:", diffErr.Integral()
    
    
    m1,m1err,m2,m2err = calcMoments(observed[0][-1])
    obs_moments[0][-1] = {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}
    
    if toyOpt == "Pseudodata" and ntoys > 0:
        print "Creating pseudodata toys"
        # Create Pseudodata toys here
        for t in range(1, ntoys+1):
            observed[t] = {-1: observed[0][-1].Clone(observed[0][-1].GetName()+"_toy_%d " % t)}
            for b in range(1, observed[0][-1].GetNbinsX()+1):
                content = observed[t][-1].GetBinContent(b)
                error = observed[t][-1].GetBinError(b)
                newContent = rnd.Gaus(content,error)
                #newContent = rnd.Gaus(content,error) if not usePoissonErrors else rnd.Gaus(content, content**0.5)
                newError = (newContent/content)*error
                
                observed[t][-1].SetBinContent(b, newContent)
                observed[t][-1].SetBinError(b, newError)
        
            m1,m1err,m2,m2err = calcMoments(observed[t][-1])
            obs_moments[t] = {-1: {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}}

else:
    for m in actualMasses:
        observed[0][m] = actualSig[m].Clone("observed_mt%d" % (int(m*10)))
        if usePoissonErrors:
            for b in range(1, observed[0][m].GetNbinsX()+1):
                observed[0][m].SetBinError(b, observed[0][m].GetBinContent(b)**0.5)

        # Apply uncertainty scaling here (if applicable)
        if useUncScaling:
            for b in range(1, observed[0][m].GetNbinsX()+1):
                observed[0][m].SetBinError(b, observed[0][m].GetBinError(b) * scaling.GetBinContent(b))


        observed[0][m].Add(bkgd)
        #if normalizeExpected:
        #    observed[0][m].Scale(actualSig[172.5].Integral()/observed[0][m].Integral())
    
        m1,m1err,m2,m2err = calcMoments(observed[0][m])
        obs_moments[0][m] = {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}
        
        if toyOpt == "Pseudodata" and ntoys > 0:
            print "Creating pseudodata toys"
            # Create Pseudodata toys here
            for t in range(1, ntoys+1):
                if t not in observed:
                    observed[t] = {}
                    obs_moments[t] = {}
                observed[t][m] = observed[0][m].Clone(observed[0][m].GetName()+"_toy_%d " % t)
                for b in range(1, observed[0][m].GetNbinsX()+1):
                    content = observed[t][m].GetBinContent(b)
                    error = observed[t][m].GetBinError(b)
                    newContent = rnd.Gaus(content,error)
                    #newContent = rnd.Gaus(content,error) if not usePoissonErrors else rnd.Gaus(content, content**0.5)
                    newError = (newContent/content)*error
                    
                    observed[t][m].SetBinContent(b, newContent)
                    observed[t][m].SetBinError(b, newError)
           
                if normToyRate:
                    observed[t][m].Scale(observed[0][m].Integral() / observed[t][m].Integral())

                m1,m1err,m2,m2err = calcMoments(observed[t][m])
                obs_moments[t][m] = {"m1":m1,"m1err":m1err,"m2":m2,"m2err":m2err}

# https://root.cern.ch/doc/master/RooNLLVar_8cxx_source.html#l00271
def calcNLL(_expected, _observed):
    _nll = 0.

    for b in range(1,_expected.GetNbinsX()+1):
        N = _observed.GetBinContent(b)# * _observed.GetBinWidth(b)
        mu = _expected.GetBinContent(b)# * _expected.GetBinWidth(b)

        #-1*(-mu + N*log(mu) - TMath::LnGamma(N+1)) ;
        tmp = mu - N * TMath.Log(mu) - TMath.LnGamma(N+1)

        _nll += tmp

    return _nll

# http://cp3.irmp.ucl.ac.be/~delaere/BND2010/AnalysisMethods_Delaere_BND2010.pdf#page=105
def calcGaussNLL(_expected, _observed):
    _nll = 0.

    for b in range(1,_expected.GetNbinsX()+1):
        N = _observed.GetBinContent(b)
        sigma2 = _observed.GetBinError(b)**2
        mu = _expected.GetBinContent(b)

        #tmp = TMath.Log(2.*TMath.Pi() * sigma2) + (N - mu)**2 / sigma2
        tmp = TMath.Log(2.*TMath.Pi() * sigma2) + (N - mu)**2 / (2.*sigma2)
        #tmp = (N - mu)**2 / (2.*sigma2)

        _nll += tmp

    return _nll


nll = {}
minNLL = {}
deltaNLL2 = {}
nllG = {}

moment1G = {}
moment2G = {}
fit1 = {}
fit2 = {}

global m1_mt, m1_mterr, m2_mt, m2_mterr
m1_mt = {}
m1_mterr = {}
m2_mt = {}
m2_mterr = {}
# t=0 is default fit
# t>1: toys (if any)
for t in range(ntoys+1):
    nll[t] = {}
    minNLL[t] = {}
    deltaNLL2[t] = {}
    nllG[t] = {}
    m1_mt[t] = {}
    m1_mterr[t] = {}
    m2_mt[t] = {}
    m2_mterr[t] = {}
    if t == 0 or (t > 0 and toyOpt == "MCStat"):
        moment1G[t] = TGraphErrors(len(actualMasses), array('d', actualMasses), array('d', [exp_moments[t][m]["m1"] for m in actualMasses]), array('d', [0.]*len(actualMasses)), array('d', [exp_moments[t][m]["m1err"] for m in actualMasses]))
        moment2G[t] = TGraphErrors(len(actualMasses), array('d', actualMasses), array('d', [exp_moments[t][m]["m2"] for m in actualMasses]), array('d', [0.]*len(actualMasses)), array('d', [exp_moments[t][m]["m2err"] for m in actualMasses]))

        moment1G[t].SetName("moment1_%d" % t)
        moment2G[t].SetName("moment2_%d" % t)
       
        moment1G[t].Fit("pol1", "Q")
        moment2G[t].Fit("pol1", "Q")

        fit1[t] = moment1G[t].GetFunction("pol1")
        fit2[t] = moment2G[t].GetFunction("pol1")

    
    for mt_dataobs in ([-1] if useDataObs else actualMasses):
        # Moments
        if toyOpt == "MCStat":
            m1_mt[t][mt_dataobs], m1_mterr[t][mt_dataobs] = massFromMoment(fit1[t], obs_moments[0][mt_dataobs]["m1"], obs_moments[0][mt_dataobs]["m1err"])
            m2_mt[t][mt_dataobs], m2_mterr[t][mt_dataobs] = massFromMoment(fit2[t], obs_moments[0][mt_dataobs]["m2"], obs_moments[0][mt_dataobs]["m2err"])
        
        elif toyOpt == "Pseudodata":
            m1_mt[t][mt_dataobs], m1_mterr[t][mt_dataobs] = massFromMoment(fit1[0], obs_moments[t][mt_dataobs]["m1"], obs_moments[t][mt_dataobs]["m1err"])
            m2_mt[t][mt_dataobs], m2_mterr[t][mt_dataobs] = massFromMoment(fit2[0], obs_moments[t][mt_dataobs]["m2"], obs_moments[t][mt_dataobs]["m2err"])
        
        else:
            m1_mt[t][mt_dataobs], m1_mterr[t][mt_dataobs] = massFromMoment(fit1[t], obs_moments[t][mt_dataobs]["m1"], obs_moments[t][mt_dataobs]["m1err"])
            m2_mt[t][mt_dataobs], m2_mterr[t][mt_dataobs] = massFromMoment(fit2[t], obs_moments[t][mt_dataobs]["m2"], obs_moments[t][mt_dataobs]["m2err"])

        # Likelihoods
        nll[t][mt_dataobs] = []
        if useGaussML:
            for m in masses:
                if toyOpt == "MCStat":
                    nll[t][mt_dataobs].append(calcGaussNLL(expected[t][m], observed[0][mt_dataobs]))
                elif toyOpt == "Pseudodata":
                    nll[t][mt_dataobs].append(calcGaussNLL(expected[0][m], observed[t][mt_dataobs]))
                else:
                    # Just default fit
                    nll[t][mt_dataobs].append(calcGaussNLL(expected[0][m], observed[0][mt_dataobs]))
        
        else:
            for m in masses:
                if toyOpt == "MCStat":
                    nll[t][mt_dataobs].append(calcNLL(expected[t][m], observed[0][mt_dataobs]))
                elif toyOpt == "Pseudodata":
                    nll[t][mt_dataobs].append(calcNLL(expected[0][m], observed[t][mt_dataobs]))
                else:
                    # Just default fit
                    nll[t][mt_dataobs].append(calcNLL(expected[0][m], observed[0][mt_dataobs]))

        minNLL[t][mt_dataobs] = min(nll[t][mt_dataobs])
        deltaNLL2[t][mt_dataobs] = [2*(n - minNLL[t][mt_dataobs]) for n in nll[t][mt_dataobs]]

        nllG[t][mt_dataobs] = TGraph(len(masses), array('d', masses), array('d', deltaNLL2[t][mt_dataobs]))
        nllG[t][mt_dataobs].SetMarkerStyle(20)
        nllG[t][mt_dataobs].SetMarkerSize(2)
        nllG[t][mt_dataobs].SetName("nll_dataobs%s%s" % ("" if useDataObs else "_"+str(int(mt_dataobs*10)), "" if t == 0 else "_toy_%d" % t) )
        nllG[t][mt_dataobs].SetTitle("Likelihood points  %s%s;m_{t} [GeV];2#Delta NLL" % ("data_obs" if useDataObs else "%.1f" % mt_dataobs, "" if t == 0 else "  toy %d" % t))
#g.Draw("alp")


tmp = deepcopy(m1_mt)

minMT = masses[0]
maxMT = masses[-1]
deltaMT = 0.01

evalX = [minMT + deltaMT * i for i in range(int(float(maxMT - minMT)/deltaMT)+1)]
evalY = {}
evalG = {}

for t in range(ntoys+1):
    evalY[t] = {}
    evalG[t] = {}
    for m in ([-1] if useDataObs else actualMasses): 
        evalY[t][m] = []
        for x in evalX:
            evalY[t][m].append(nllG[t][m].Eval(x, 0, "S"))

        evalG[t][m] = TGraph(len(evalX), array('d', evalX), array('d', evalY[t][m]))
        evalG[t][m].SetLineColor(kRed)
        evalG[t][m].SetLineWidth(2)
        evalG[t][m].SetMarkerColor(kRed)
        evalG[t][m].SetName("nll_eval%s%s" % ("_dataobs" if useDataObs else "_"+str(int(m*10)), "" if t == 0 else "_toy_%d" % t))
        evalG[t][m].SetTitle("Likelihood scan  %s%s;m_{t} [GeV];2#Delta NLL" % ("data_obs" if useDataObs else "%.1f" % m, "" if t == 0 else "  toy %d" % t) )


def find_index_nearest(arrayIn,value):
    array = np.array(arrayIn)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return (idx-1)
    else:
        return idx

#def old_find_index_nearest(arrayIn, value):
#    idx = -1
#
#    for i,x in enumerate(arrayIn):
#        if x > value:
#            break
#
#    idx = i if (abs(arrayIn[i] - value) < abs(arrayIn[i-1] - value)) else (i-1)
#    return idx

def find_crossing(mtarrayIn, nllarrayIn, value=1):
    mtarray = np.array(mtarrayIn)
    nllarray = np.array(nllarrayIn)
    minNLL = min(nllarray)
    minNLLindex = np.argmin(nllarray)

    leftNLL = nllarray[:minNLLindex]
    leftNLL = leftNLL[::-1]     # Must be in increasing order
    rightNLL = nllarray[minNLLindex:]

    leftMT = mtarray[:minNLLindex]
    leftMT = leftMT[::-1]
    rightMT = mtarray[minNLLindex:]

    rightIndex = find_index_nearest(rightNLL, minNLL+1.)
    leftIndex = find_index_nearest(leftNLL, minNLL+1.)

    #mtDown = leftMT[leftIndex][0]
    #mtUp = rightMT[rightIndex][0]
    #mtCentral = mtarray[minNLLindex][0]

    _mtDown = leftMT[leftIndex]
    _mtUp = rightMT[rightIndex]
    _mtCentral = mtarray[minNLLindex]

    return (_mtDown,_mtCentral,_mtUp)


mtDown = {0:{}}
mtCentral = {0:{}}
mtUp = {0:{}}
mtErrUp = {0:{}}
mtErrDown = {0:{}}


mtToyH = {}
moment1_mtToyH = {}
moment2_mtToyH = {}

pullH = {}
moment1_pullH = {}
moment2_pullH = {}

mg = {}

txt=TLatex()
txt.SetNDC(True)
txt.SetTextFont(43)
txt.SetTextSize(40)
txt.SetTextAlign(12)

m1_mt = deepcopy(tmp)

os.system("mkdir -p %s" % outDir)
c = TCanvas("c","c",1600,1200)

nFailedToys = 0
failedToys = []


def drawFitHist(h, outFile):
    h.Draw("hist")
    h.Fit("gaus", "Q")
    func = h.GetFunction("gaus")
    mean = func.GetParameter(1)
    meanErr = func.GetParError(1)
    sigma = func.GetParameter(2)
    sigmaErr = func.GetParError(2)
    
    func.Draw("same")
    txt.DrawLatex(0.65, 0.75, "Mean: %.3f #pm %.3f" % (mean, meanErr))
    txt.DrawLatex(0.65, 0.7, " RMS: %.3f #pm %.3f" % (sigma, sigmaErr))

    c.SaveAs(outFile)

    return mean,meanErr,sigma,sigmaErr


if useDataObs:
    mtToyH[0] = TH1D("%s_likelihood_mtToys" % obs, "%s  toy masses" % obsTitle[obs], 200, 169.5, 175.5)
    mtToyH[0].GetXaxis().SetTitle("m_{t} [GeV]")
    mtToyH[0].GetYaxis().SetTitle("Toys")

    moment1_mtToyH[0] = TH1D("%s_moment1_mtToys" % obs, "%s  Moment 1 toy masses" % obsTitle[obs], 200, 169.5, 175.5)
    moment1_mtToyH[0].GetXaxis().SetTitle("m_{t} [GeV]")
    moment1_mtToyH[0].GetYaxis().SetTitle("Toys")
    
    moment2_mtToyH[0] = TH1D("%s_moment2_mtToys" % obs, "%s  Moment 2 toy masses" % obsTitle[obs], 200, 169.5, 175.5)
    moment2_mtToyH[0].GetXaxis().SetTitle("m_{t} [GeV]")
    moment2_mtToyH[0].GetYaxis().SetTitle("Toys")
    
    pullH[0] = TH1D("%s_likelihood_pulls" % obs, "%s  Likelihood Toy Pulls" % obsTitle[obs], 200, -5, 5)
    pullH[0].GetXaxis().SetTitle("#frac{m_{t}^{toy} - m_{t}^{fit}}{#sigma_{m_{t}^{fit}}}")
    pullH[0].GetYaxis().SetTitle("Toys")

    moment1_pullH[0] = TH1D("%s_moment1_pulls" % obs, "%s  Moment 1 Toy Pulls" % obsTitle[obs], 200, -5, 5)
    moment1_pullH[0].GetXaxis().SetTitle("#frac{m_{t}^{toy} - m_{t}^{fit}}{#sigma_{m_{t}^{fit}}}")
    moment1_pullH[0].GetYaxis().SetTitle("Toys")
    
    moment2_pullH[0] = TH1D("%s_moment2_pulls" % obs, "%s  Moment 2 Toy Pulls" % obsTitle[obs], 200, -5, 5)
    moment2_pullH[0].GetXaxis().SetTitle("#frac{m_{t}^{toy} - m_{t}^{fit}}{#sigma_{m_{t}^{fit}}}")
    moment2_pullH[0].GetYaxis().SetTitle("Toys")
   

    for t in range(ntoys+1):
        mtDown[t] = {}
        mtCentral[t] = {}
        mtUp[t] = {}
        mtErrDown[t] = {}
        mtErrUp[t] = {}
        try:
            mtDown[t][-1], mtCentral[t][-1], mtUp[t][-1] = find_crossing(evalX, evalY[t][-1], 1.)
        except IndexError:
            # Cant find crossing
            print "Can't evaluate toy %d" % t
            nFailedToys += 1
            continue
        mtErrDown[t][-1] = abs(mtDown[t][-1] - mtCentral[t][-1])
        mtErrUp[t][-1] = abs(mtUp[t][-1] - mtCentral[t][-1])
        if t > 0:
            if mtErrDown[t][-1] > 0.0 and mtErrUp[t][-1] > 0.0:
                mtToyH[0].Fill(mtCentral[t][-1])
                moment1_mtToyH[0].Fill(m1_mt[t][-1])
                moment2_mtToyH[0].Fill(m2_mt[t][-1])
               
                res = mtCentral[t][-1] - mtCentral[0][-1]
                pull = res/mtErrUp[0][-1] if res > 0 else res/mtErrDown[0][-1]
                pullH[0].Fill(pull)

                moment1_pullH[0].Fill((m1_mt[t][-1] - m1_mt[0][-1]) / m1_mterr[0][-1])
                moment2_pullH[0].Fill((m2_mt[t][-1] - m2_mt[0][-1]) / m2_mterr[0][-1])
            else:
                print "Toy %d has failed %s error" % (t, "DOWN" if mtErrDown[t][-1] < 0.0001 else "UP")
                nFailedToys += 1
                failedToys.append(t)

    if ntoys > 0:
        if nFailedToys > 0:
            print "Toy fits failed: %d" % nFailedToys
        else:
            print "All toy fits succeeded!"


    print "Likelihood:\t%.3f  -%.3f / +%.3f" % (mtCentral[0][-1], mtErrDown[0][-1], mtErrUp[0][-1])
    print "Moment 1:\t%.3f  -%.3f / +%.3f" % (m1_mt[0][-1], m1_mterr[0][-1], m1_mterr[0][-1])
    print "Moment 2:\t%.3f  -%.3f / +%.3f" % (m2_mt[0][-1], m2_mterr[0][-1], m2_mterr[0][-1])
    print ""
    if ntoys > 0:
        #print "Likelihood toys:\t%.3f +- %.3f" % (mtToyH[0].GetMean(), mtToyH[0].GetRMS())
#        mtToyH[0].Draw("hist")
#        mtToyH[0].Fit("gaus")
#        likelihoodFunc = mtToyH[0].GetFunction("gaus")
#        likelihoodMean = likelihoodFunc.GetParameter(1)
#        likelihoodMeanErr = likelihoodFunc.GetParError(1)
#        likelihoodSigma = likelihoodFunc.GetParameter(2)
#        likelihoodSigmaErr = likelihoodFunc.GetParError(2)
#    
#        likelihoodFunc.Draw("same")
#        txt.DrawLatex(0.55, 0.55, "Mean: %.2f #pm %.2f" % (likelihoodMean, likelihoodMeanErr))
#        txt.DrawLatex(0.55, 0.45, " RMS: %.2f #pm %.2f" % (likelihoodSigma, likelihoodSigmaErr))
#
#        c.SaveAs("%s/likelihood_mtToys.png" % outDir)
        lhMean,lhMeanErr,lhSigma,lhSigmaErr = drawFitHist(mtToyH[0], "%s/likelihood_mtToys.png" % outDir)
        print "Likelihood toys:\t%.3f +- %.3f" % (lhMean, lhSigma)
        
        lhPullMean,lhPullMeanErr,lhPullSigma,lhPullSigmaErr = drawFitHist(pullH[0], "%s/likelihood_pulls.png" % outDir)
        print "Likelihood pulls:\t%.3f +- %.3f" % (lhPullMean, lhPullSigma)

#        print "moment1 toys:\t%.3f +- %.3f" % (moment1_mtToyH[0].GetMean(), moment1_mtToyH[0].GetRMS())
#        moment1_mtToyH[0].Draw("hist")
#        c.SaveAs("%s/moment1_mtToys.png" % outDir)
        
        m1Mean,m1MeanErr,m1Sigma,m1SigmaErr = drawFitHist(moment1_mtToyH[0], "%s/moment1_mtToys.png" % outDir)
        print "moment1 toys:\t%.3f +- %.3f" % (m1Mean, m1Sigma)
        m1PullMean,m1PullMeanErr,m1PullSigma,m1PullSigmaErr = drawFitHist(moment1_pullH[0], "%s/moment1_pulls.png" % outDir)
        print "moment1 pulls:\t%.3f +- %.3f" % (m1PullMean, m1PullSigma)
        
#        print "moment2 toys:\t%.3f +- %.3f" % (moment2_mtToyH[0].GetMean(), moment2_mtToyH[0].GetRMS())
#        moment2_mtToyH[0].Draw("hist")
#        c.SaveAs("%s/moment2_mtToys.png" % outDir)
        
        m2Mean,m2MeanErr,m2Sigma,m2SigmaErr = drawFitHist(moment2_mtToyH[0], "%s/moment2_mtToys.png" % outDir)
        print "moment2 toys:\t%.3f +- %.3f" % (m2Mean, m2Sigma)
        m2PullMean,m2PullMeanErr,m2PullSigma,m2PullSigmaErr = drawFitHist(moment2_pullH[0], "%s/moment2_pulls.png" % outDir)
        print "moment2 pulls:\t%.3f +- %.3f" % (m2PullMean, m2PullSigma)

    
    mg[0] = TMultiGraph()
    mg[0].SetName("mg_dataobs")
    mg[0].SetTitle("Likelihood scan  data_obs;m_{t} [GeV]; 2#Delta NLL")
    mg[0].Add(nllG[0][-1], "p")
    mg[0].Add(evalG[0][-1], "lp")

    if ntoys > 0:
        for t in range(1,ntoys+1):
            mg[t] = TMultiGraph()
            mg[t].SetName("mg_dataobs_toy_%d" % t)
            mg[t].SetTitle("Likelihood scan  data_obs  Toy %d;m_{t} [GeV]; 2#Delta NLL" % t)
            mg[t].Add(nllG[t][-1], "p")
            mg[t].Add(evalG[t][-1], "lp")




    mg[0].Draw("a")

    txt.DrawLatex(0.35, 0.75, "m_{t} = %.3f  -%.3f / +%.3f GeV" % (mtCentral[0][-1], mtErrDown[0][-1], mtErrUp[0][-1]))
    

    c.SaveAs("%s/dataobs_manualLikelihood.png" % outDir)

    
else: 
    for m in actualMasses[1:-1]:
        #print "About to find crossing for mt = %.1f" % m
        mtDown[m], mtCentral[m], mtUp[m] = find_crossing(evalX, evalY[0][m], 1.)
        mtErrDown[m] = abs(mtDown[m] - mtCentral[m])
        mtErrUp[m] = abs(mtUp[m] - mtCentral[m])
        print "At %.3f:\t%.3f  -%.3f / +%.3f" % (m, mtCentral[m], mtErrDown[m], mtErrUp[m])

        mt10 = int(m*10)
        mg[m] = TMultiGraph()
        mg[m].SetName("mg_%d" % mt10)
        mg[m].SetTitle("Likelihood scan  %.1f;m_{t} [GeV]; 2#Delta NLL" % m)
        mg[m].Add(nllG[0][m], "p")
        mg[m].Add(evalG[0][m], "lp")

        mg[m].Draw("a")

        c.SaveAs("%s/mt%d_manualLikelihood.png" % (outDir,mt10))


# Save masses
with open("%s/masses.txt" % outDir, "w") as f:
    if useDataObs:
        f.write("data_obs\n")
        f.write("Likelihood:\t%.3f  -%.3f / +%.3f\n" % (mtCentral[0][-1], mtErrDown[0][-1], mtErrUp[0][-1]))
        f.write("Moment 1:\t%.3f  -%.3f / +%.3f\n" % (m1_mt[0][-1], m1_mterr[0][-1], m1_mterr[0][-1]))
        f.write("Moment 2:\t%.3f  -%.3f / +%.3f\n" % (m2_mt[0][-1], m2_mterr[0][-1], m2_mterr[0][-1]))

        if ntoys > 0:
            f.write("\n%d %s Toys\n" % (ntoys, toyOpt))
            lhMean,lhMeanErr,lhSigma,lhSigmaErr
            f.write("Likelihood:\t%.3f +- %.3f\n" % (lhMean, lhSigma))
            f.write("Moment 1:\t%.3f +- %.3f\n" % (m1Mean, m1Sigma))
            f.write("Moment 2:\t%.3f +- %.3f\n" % (m2Mean, m2Sigma))

#            f.write("Likelihood:\t%.3f +- %.3f\n" % (mtToyH[0].GetMean(), mtToyH[0].GetRMS()))
#            f.write("Moment 1:\t%.3f +- %.3f\n" % (moment1_mtToyH[0].GetMean(), moment1_mtToyH[0].GetRMS()))
#            f.write("Moment 2:\t%.3f +- %.3f\n" % (moment2_mtToyH[0].GetMean(), moment2_mtToyH[0].GetRMS()))
    else:
        for m in actualMasses[1:-1]:
            f.write("data_obs = %.3f" % m)
            f.write("Likelihood:\t%.3f  -%.3f / +%.3f\n" % (mtCentral[m], mtErrDown[m], mtErrUp[m]))
            f.write("Moment 1:\t%.3f  -%.3f / +%.3f\n" % (m1_mt[0][-1], m1_mterr[0][-1], m1_mterr[0][-1]))
            f.write("Moment 2:\t%.3f  -%.3f / +%.3f\n\n" % (m2_mt[0][-1], m2_mterr[0][-1], m2_mterr[0][-1]))

# Save toys
print "Saving toys...",
if ntoys > 0:
    if saveAllToys:
        toyF = TFile.Open("%s/alltoys.root" % outDir, "recreate")
        for t in range(1,ntoys+1):
            observed[t][-1].Write()
        toyF.Close()
    

    toyF = TFile.Open("%s/toyH.root" % outDir, "recreate")
    mtToyH[0].Write()
    moment1_mtToyH[0].Write()
    moment2_mtToyH[0].Write()
    pullH[0].Write()
    moment1_pullH[0].Write()
    moment2_pullH[0].Write()
    toyF.Close()

    if nFailedToys > 0:
        with open("%s/failedToys.txt" % outDir, "w") as f:
            for t in failedToys:
                f.write("%d\n" % t)

print "done"
print "\ndeltaNLL2[0][1] =", deltaNLL2[0][-1]
sys.stdout.flush()
sys.exit()

calib_masses = actualMasses[1:-1] 
#calib_masses = masses[1:-1]


# Calibration curve
calibG = TGraphAsymmErrors(len(calib_masses), array('d', calib_masses), array('d', [mtCentral[mt] for mt in calib_masses]), array('d', [0.]*len(calib_masses)), array('d', [0.]*len(calib_masses)), array('d', [mtErrDown[mt] for mt in calib_masses]), array('d', [mtErrUp[mt] for mt in calib_masses]))

biasG = TGraphAsymmErrors(len(calib_masses), array('d', calib_masses), array('d', [(mtCentral[mt] - mt) for mt in calib_masses]), array('d', [0.]*len(calib_masses)), array('d', [0.]*len(calib_masses)), array('d', [mtErrDown[mt] for mt in calib_masses]), array('d', [mtErrUp[mt] for mt in calib_masses]))

calibG.SetName("calibration")
calibG.SetTitle("Calibration;m_{t}^{MC} [GeV];m_{t}^{fit} [GeV]")
calibG.SetMarkerStyle(20)
calibG.SetMarkerSize(2)

biasG.SetName("bias")
biasG.SetTitle("Bias;m_{t}^{MC} [GeV];m_{t}^{fit} - m_{t}^{MC} [GeV]")
biasG.SetMarkerStyle(20)
biasG.SetMarkerSize(2)


os.system("mkdir -p %s" % outDir)

f = TFile.Open("%s/graphs.root" % outDir, "recreate")
calibG.Write()
biasG.Write()
for m in ([-1] if useDataObs else actualMasses):
    nllG[0][m].Write()
    evalG[0][m].Write()

f.Close()


txt=TLatex()
txt.SetNDC(True)
txt.SetTextFont(43)
txt.SetTextSize(35)
txt.SetTextAlign(12)


calibG.Draw("ap")
calibG.Fit("pol1", "Q")
calibFunc = calibG.GetFunction("pol1")

calibInt = calibFunc.GetParameter(0)
calibIntErr = calibFunc.GetParError(0)
calibSlope = calibFunc.GetParameter(1)
calibSlopeErr = calibFunc.GetParError(1)

txt.DrawLatex(0.45, 0.25, "m_{t}^{fit} = (%.3f #pm %.3f) m_{t}^{MC} + (%.3f #pm %.3f)" % (calibSlope, calibSlopeErr, calibInt, calibIntErr))

c.SaveAs("%s/calibration.png" % outDir)

biasG.Draw("ap")
biasG.Fit("pol1", "Q")
biasFunc = biasG.GetFunction("pol1")

biasInt = biasFunc.GetParameter(0)
biasIntErr = biasFunc.GetParError(0)
biasSlope = biasFunc.GetParameter(1)
biasSlopeErr = biasFunc.GetParError(1)

#txt.DrawLatex(0.35, 0.85, "Bias = (%.3f #pm %.3f) m_{t}^{MC} + (%.3f #pm %.3f)" % (biasSlope, biasSlopeErr, biasInt, biasIntErr))
txt.DrawLatex(0.45, 0.75, "slope = %.3f #pm %.3f    offset = %.3f #pm %.3f" % (biasSlope, biasSlopeErr, biasInt, biasIntErr))

c.SaveAs("%s/bias.png" % outDir)



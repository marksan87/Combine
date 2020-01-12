#!/usr/bin/env python
from __future__ import division
import ROOT
from ROOT import *
from argparse import ArgumentParser
from array import array
import os
import sys
from time import sleep
from pprint import pprint
import pickle
import gzip
from copy import deepcopy
from binning10GeV import templateBinning

gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.PyConfig.DisableRootLogon = True

allSystematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS',"MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"]

oneSidedSysts = ["CRerdON", "CRGluon", "CRQCD", "amcanlo", "madgraph", "herwigpp" ]

#systematicsToScale = {"TTbar":["Q2","MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"],
#                      "ST_tW":["Q2"]}

systematicsToScale = {\
                    "TTbar": ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", "hdamp","UE","CRerdON","CRGluon","CRQCD","amcanlo","madgraph","herwigpp"],
                    "ST_tW": ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "Q2", "isr", "fsr", "hdamp", "DS"],
        }

_diffDists = ["ptll_M0_E0", "ptll_M0_E1", "ptll_M0_E2", "ptll_M1_E0", "ptll_M1_E1", "ptll_M1_E2", "ptll_M2_E0", "ptll_M2_E1", "ptll_M2_E2"]
_observables = ['ptll', 'Mll', 'ptpos', 'ptneg', 'Epos', 'Eneg', 'ptp_ptm', 'Ep_Em', "leadLepPt", "leadJetPt"] + _diffDists
obsTitle = {"ptll":"p_{T}(ll)", 
            "ptpos":"p_{T}(l^{+})", 
            "ptneg":"p_{T}(l^{-})", 
            "Epos":"E(l^{+})", 
            "Eneg":"E(l^{-})", 
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

signal = ["TTbar", "ST_tW"]
background = ["DY", "TTV", "WJets", "Diboson", "ST_bkgd"]
systematics = {"nominal":"hists", \
               "topptUp":"histstoppt", \
               "topptDown":"hists", \
               "pileupUp":"histsPU_up", \
               "pileupDown":"histsPU_down", \
               "Q2Up":"histsQ2_up", \
               "Q2Down":"histsQ2_down", \
               "PdfUp":"histsPdf_up", \
               "PdfDown":"histsPdf_down", \
               "LumiUp":"histsLumi_up", \
               "LumiDown":"histsLumi_down", \
               "EleIDEffUp":"histsEleIDEff_up", \
               "EleIDEffDown":"histsEleIDEff_down", \
               "EleRecoEffUp":"histsEleRecoEff_up", \
               "EleRecoEffDown":"histsEleRecoEff_down", \
               "EleScaleUp":"histsEleScale_up", \
               "EleScaleDown":"histsEleScale_down", \
               "EleSmearUp":"histsEleSmear_up", \
               "EleSmearDown":"histsEleSmear_down", \
               "MuIDEffUp":"histsMuIDEff_up", \
               "MuIDEffDown":"histsMuIDEff_down", \
               "MuIsoEffUp":"histsMuIsoEff_up", \
               "MuIsoEffDown":"histsMuIsoEff_down", \
               "MuTrackEffUp":"histsMuTrackEff_up", \
               "MuTrackEffDown":"histsMuTrackEff_down", \
               "MuScaleUp":"histsMuScale_up", \
               "MuScaleDown":"histsMuScale_down", \
               "TrigEffUp":"histsTrigEff_up", \
               "TrigEffDown":"histsTrigEff_down", \
               "BTagSFUp":"histsBTagSF_up", \
               "BTagSFDown":"histsBTagSF_down", \
               "JECUp":"histsJEC_up", \
               "JECDown":"histsJEC_down", \
               "JERUp":"histsJER_up", \
               "JERDown":"histsJER_down", \
               "isrUp":"histsisr_up", \
               "isrDown":"histsisr_down", \
               "fsrUp":"histsfsr_up", \
               "fsrDown":"histsfsr_down", \
               "DSUp":"histsDS", \
               "DSDown":"hists", \
               "hdampUp":"histshdamp_up", \
               "hdampDown":"histshdamp_down", \
               "UEUp":"histsUE_up", \
               "UEDown":"histsUE_down", \
               "CRerdONUp":"histsCRerdON", \
               "CRerdONDown":"hists", \
               "CRGluonUp":"histsCRGluon", \
               "CRGluonDown":"hists", \
               "CRQCDUp":"histsCRQCD", \
               "CRQCDDown":"hists", \
               "amcanloUp":"histsamcanlo", \
               "amcanloDown":"hists", \
               "madgraphUp":"histsmadgraph", \
               "madgraphDown":"hists", \
               "herwigppUp":"histsherwigpp", \
               "herwigppDown":"hists", \
#               "MEscale1Up":"histsMEscale1", \
#               "MEscale1Down":"hists", \
#               "MEscale2Up":"histsMEscale2", \
#               "MEscale2Down":"hists", \
#               "MEscale3Up":"histsMEscale3", \
#               "MEscale3Down":"hists", \
#               "MEscale4Up":"histsMEscale4", \
#               "MEscale4Down":"hists", \
#               "MEscale5Up":"histsMEscale5", \
#               "MEscale5Down":"hists", \
#               "MEscale6Up":"histsMEscale6", \
#               "MEscale6Down":"hists", \
               }

separateSystSamples = {
    "TTbar":['isr','fsr','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp'],
    "ST_tW":['isr','fsr','hdamp','DS','hdamp', 'Q2']
}
ttOnlySysts = ['toppt','Pdf','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp',"MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"]
tWOnlySysts = ['DS']

PU_mt1725only = False 

def ttxsec(mt, xsecref = 803, a1 = -0.745047, a2 = 0.127417):
    mtref = 172.5
    return xsecref * (mtref/mt)**4 * (1. + a1 * (mt - mtref)/mtref + a2 * ( (mt - mtref)/mtref )**2)


def clamp(val, minV, maxV):
    # Clamps a value between min and max 
    return max(minV, min(val,maxV))


def systMorph(histMT, nominalMeanBin, diffHist, name, useVariableBinning = False):
    mtMeanBin = histMT.FindBin(histMT.GetMean())
    binShift = 0
    if mtMeanBin != nominalMeanBin and not useVariableBinning:
        binShift = mtMeanBin - nominalMeanBin
#        print "Mean shifted by %s" % binShift

    Nbins = diffHist.GetNbinsX()
    morphed = histMT.Clone(name)
    for b in xrange(1, Nbins+1):
        morphed.SetBinContent(b, morphed.GetBinContent(b) * diffHist.GetBinContent( clamp(b-binShift, 1, Nbins) ))
        morphed.SetBinError(b, morphed.GetBinError(b) * diffHist.GetBinContent( clamp(b-binShift, 1, Nbins) ))
#        morphed.SetBinContent(b, morphed.GetBinContent(b) * diffHist.GetBinContent( b ))
#        morphed.SetBinError(b, morphed.GetBinError(b) * diffHist.GetBinContent( b ))

    return morphed


def makeHist(fName, name, title, config, masses = [], systematic = "", weight="weight"):
    hist = {}
    if len(masses) > 0:
        for m in masses:
            f = TFile.Open("%s%d.root" % (fName,m), "read")
            tree = f.Get("goodEvents")
            for obs,vals in config.items():
                if obs not in hist: hist[obs] = {}
                hist[obs][m] = TH1F("%s%s%d%s" % (obs,name,m,systematic), "%s  %s  m_{t} = %.1f GeV" % (vals["title"], title, m/decimalScaling), vals["nbins"], vals["min"], vals["max"])
                tree.Draw("%s>>%s" % (obs,hist[obs][m].GetName()), "%s * (%s > %f && %s < %f)" % (weight, obs, vals["min"], obs, vals["max"]))
                hist[obs][m].SetDirectory(0)
            
            f.Close() 
    else:
        f = TFile.Open("%s.root" % fName, "read")
        tree = f.Get("goodEvents")
        for obs,vals in config.items():
            hist[obs] = TH1F("%s%s%s" % (obs,name,systematic), "%s  %s" % (vals["title"], title), vals["nbins"], vals["min"], vals["max"])
            tree.Draw("%s>>%s" % (obs,hist[obs].GetName()), "%s * (%s > %f && %s < %f)" % (weight, obs, vals["min"], obs, vals["max"]))
            hist[obs].SetDirectory(0)
        
        f.Close() 
    return hist


def cdfMorphTemplates(templates, morph_masses, name, title, Neff = 0, precision=1, systematic="", interp = "pol3", morphRates = False, verbosity=1):
    morph = {}
    fit = {}
    fitFunc = {}
    binG = {}
    binMorphG = {}
    rates = {}
    normalized = {}
    cdf = {}
    cdfBinG = {}
    cdfMorph = {}
    signalType = name[-2:]
    
    actual_masses = sorted(templates.keys())
    if Neff == 0:
        # Use effective entries from 172.5 template
        Neff = templates[int(172.5 * decimalScaling)].GetEffectiveEntries()

    if verbosity > 0: print "CDF Morphing  %s  %s with Neff = %d" % (name,systematic[1:],Neff)
    obs = name[4:-3]
    for m in actual_masses:
        if verbosity > 1: print "Now on %s  %s  mt = %.1f" % (name,systematic,m)
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone("__normalized_%s_%s_mt%d" % (name,systematic,m))
        normalized[m].SetDirectory(0)
        if not morphRates:
            normalized[m].Scale(1.0 / rates[m])

        # Create cdf
        cdf[m] = normalized[m].GetCumulative()

    # Morph bins of the cdf
    for b in xrange(1, cdf[actual_masses[0]].GetNbinsX()+1):
        binG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [normalized[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [normalized[m].GetBinError(b) for m in actual_masses]))
        binG[b].SetName("%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        binG[b].SetTitle("%s  Bin %d" % (obs,b))
        binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binG[b].GetYaxis().SetTitle("Entries")
        binG[b].GetYaxis().SetTitleOffset(1.3)
        
        cdfBinG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [cdf[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [cdf[m].GetBinError(b) for m in actual_masses]))
        cdfBinG[b].SetName("cdf_%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        cdfBinG[b].SetTitle("CDF  %s  Bin %d" % (obs,b))
        cdfBinG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        cdfBinG[b].GetYaxis().SetTitle("Entries")
        cdfBinG[b].GetYaxis().SetTitleOffset(1.3)
        cdfBinG[b].Draw()
        fit[b] = cdfBinG[b].Fit(interp, "S"+ ("" if verbosity > 1 else "Q"))
        fitFunc[b] = cdfBinG[b].GetFunction(interp)
        errors = array('d', [0.] * len(morph_masses))
        fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/decimalScaling for m in morph_masses]), errors, 2./3., False)
        for i,m in enumerate(morph_masses):
            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
            if m not in cdfMorph:
                if len(cdf[actual_masses[0]].GetXaxis().GetXbins()) == 0:
                    cdfMorph[m] = TH1F("__cdf_%s%d%s" % (name,m,systematic), "CDF %s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
                else:
                    cdfMorph[m] = TH1F("__cdf_%s%d%s" % (name,m,systematic), "CDF %s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
            cdfMorph[m].SetBinContent(b, fitFunc[b].Eval(m/decimalScaling))
   
    # 'Differentiate' back to pdf
    for i,m in enumerate(morph_masses):
        if m not in morph:
            if len(cdf[actual_masses[0]].GetXaxis().GetXbins()) == 0:
                morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
            else:
                morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
        
        # Set bin contents
        morph[m].SetBinContent(1, cdfMorph[m].GetBinContent(1))
        for b in xrange(2, cdf[actual_masses[0]].GetNbinsX()+1):
            morph[m].SetBinContent(b, (cdfMorph[m].GetBinContent(b) - cdfMorph[m].GetBinContent(b-1)) )
        
        # Set bin errors
        N = morph[m].Integral()
        for b in xrange(1, cdfMorph[actual_masses[0]].GetNbinsX()+1):
            morph[m].SetBinError(b, (morph[m].GetBinContent(b)*N/Neff)**0.5)

    for b in xrange(1, morph[morph_masses[0]].GetNbinsX()+1):
        binMorphG[b] = TGraphErrors(len(morph_masses), array('d', [m/decimalScaling for m in morph_masses]), array('d', [morph[m].GetBinContent(b) for m in morph_masses]), array('d', [0.] * len(morph_masses)), array('d', [morph[m].GetBinError(b) for m in morph_masses]))
        binMorphG[b].SetName("%s%smorphed_bin_%d" % (name,"" if systematic == "" else systematic+"_",b))
        binMorphG[b].SetTitle("%s Morphed Bin %d" % ("t#bar{t}" if signalType == "tt" else signalType, b) )
        binMorphG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binMorphG[b].GetYaxis().SetTitle("Entries")
        binMorphG[b].GetYaxis().SetTitleOffset(1.4)

    
    #print "Done with %s  %s" % (name,systematic[1:])    
    if not morphRates:
        for m in morph_masses:
            morph[m].Scale(rates[int(decimalScaling*172.5)])


    return morph,binG,binMorphG,cdfMorph,cdfBinG


def avgCdfMorphTemplates(templates, morph_masses, name, title, Neff = 0, precision=1, systematic="", interp = "pol3", morphRates = False, verbosity=1):
    # Morph using average of forward & backward cdf 
    global morph, forward_cdfMorph, backward_cdfMorph, forward_cdf, backward_cdf
    morph = {}
    forward_fit = {}
    forward_fitFunc = {}
    backward_fit = {}
    backward_fitFunc = {}
    binG = {}
    binMorphG = {}
    rates = {}
    normalized = {}
    forward_cdf = {}
    forward_cdfBinG = {}
    forward_cdfMorph = {}
    backward_cdf = {}
    backward_cdfBinG = {}
    backward_cdfMorph = {}
    
    signalType = name[-2:]
    nbins = templates[int(172.5 * decimalScaling)].GetNbinsX()

    actual_masses = sorted(templates.keys())
    if Neff == 0:
        # Use effective entries from 172.5 template
        Neff = templates[int(172.5 * decimalScaling)].GetEffectiveEntries()

    if verbosity > 0: print "CDF Morphing  %s  %s with Neff = %d" % (name,systematic[1:],Neff)
    obs = name[4:-3]
    for m in actual_masses:
        if verbosity > 1: print "Now on %s  %s  mt = %.1f" % (name,systematic,m)
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone("__normalized_%s_%s_mt%d" % (name,systematic,m))
        normalized[m].SetDirectory(0)
        if not morphRates:
            normalized[m].Scale(1.0 / rates[m])

        # Create cdf
        forward_cdf[m]  = normalized[m].GetCumulative(True, "_fw")
        backward_cdf[m] = normalized[m].GetCumulative(False, "_bw")

    # Morph bins of the cdf
    for b in xrange(1, nbins+1):
        binG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [normalized[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [normalized[m].GetBinError(b) for m in actual_masses]))
        binG[b].SetName("%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        binG[b].SetTitle("%s  Bin %d" % (obs,b))
        binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binG[b].GetYaxis().SetTitle("Entries")
        binG[b].GetYaxis().SetTitleOffset(1.3)
        
        forward_cdfBinG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [forward_cdf[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [forward_cdf[m].GetBinError(b) for m in actual_masses]))
        forward_cdfBinG[b].SetName("fw_cdf_%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        forward_cdfBinG[b].SetTitle("Forward CDF  %s  Bin %d" % (obs,b))
        forward_cdfBinG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        forward_cdfBinG[b].GetYaxis().SetTitle("Entries")
        forward_cdfBinG[b].GetYaxis().SetTitleOffset(1.3)
        forward_cdfBinG[b].Draw()
        forward_fit[b] = forward_cdfBinG[b].Fit(interp, "S"+ ("" if verbosity > 1 else "Q"))
        forward_fitFunc[b] = forward_cdfBinG[b].GetFunction(interp)
        forward_errors = array('d', [0.] * len(morph_masses))
        forward_fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/decimalScaling for m in morph_masses]), forward_errors, 2./3., False)
        
        backward_cdfBinG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [backward_cdf[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [backward_cdf[m].GetBinError(b) for m in actual_masses]))
        backward_cdfBinG[b].SetName("bw_cdf_%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        backward_cdfBinG[b].SetTitle("Backward CDF  %s  Bin %d" % (obs,b))
        backward_cdfBinG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        backward_cdfBinG[b].GetYaxis().SetTitle("Entries")
        backward_cdfBinG[b].GetYaxis().SetTitleOffset(1.3)
        backward_cdfBinG[b].Draw()
        backward_fit[b] = backward_cdfBinG[b].Fit(interp, "S"+ ("" if verbosity > 1 else "Q"))
        backward_fitFunc[b] = backward_cdfBinG[b].GetFunction(interp)
        backward_errors = array('d', [0.] * len(morph_masses))
        backward_fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/decimalScaling for m in morph_masses]), backward_errors, 2./3., False)
        for i,m in enumerate(morph_masses):
            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
            if m not in forward_cdfMorph:
                if len(forward_cdf[actual_masses[0]].GetXaxis().GetXbins()) == 0:
                    forward_cdfMorph[m] = TH1F("__fw_cdf_%s%d%s" % (name,m,systematic), "Forward CDF %s  m_{t} = %s GeV" % (title, massStr), nbins, normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
                    backward_cdfMorph[m] = TH1F("__bw_cdf_%s%d%s" % (name,m,systematic), "Backward CDF %s  m_{t} = %s GeV" % (title, massStr), nbins, normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
                else:
                    forward_cdfMorph[m] = TH1F("__fw_cdf_%s%d%s" % (name,m,systematic), "Forward CDF %s  m_{t} = %s GeV" % (title, massStr), nbins, array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
                    backward_cdfMorph[m] = TH1F("__bw_cdf_%s%d%s" % (name,m,systematic), "Backward CDF %s  m_{t} = %s GeV" % (title, massStr), nbins, array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
            forward_cdfMorph[m].SetBinContent(b, forward_fitFunc[b].Eval(m/decimalScaling))
            backward_cdfMorph[m].SetBinContent(b, backward_fitFunc[b].Eval(m/decimalScaling))
  
    
    # 'Differentiate' back to pdf
    for i,m in enumerate(morph_masses):
        if m not in morph:
            if len(forward_cdf[actual_masses[0]].GetXaxis().GetXbins()) == 0:
                morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), nbins, normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
            else:
                morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), nbins, array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
        
        # Set bin contents as average of forward and backward
        morph[m].SetBinContent(1, 0.5 * (forward_cdfMorph[m].GetBinContent(1) + backward_cdfMorph[m].GetBinContent(1) - backward_cdfMorph[m].GetBinContent(2)) )
        morph[m].SetBinContent(nbins, 0.5 * (forward_cdfMorph[m].GetBinContent(nbins) - forward_cdfMorph[m].GetBinContent(nbins-1) + backward_cdfMorph[m].GetBinContent(nbins)) )
        for b in xrange(2, nbins):
            morph[m].SetBinContent(b, 0.5 * (forward_cdfMorph[m].GetBinContent(b) - forward_cdfMorph[m].GetBinContent(b-1) + backward_cdfMorph[m].GetBinContent(b) - backward_cdfMorph[m].GetBinContent(b+1)) )
        
        # Set bin errors
        N = morph[m].Integral()
        for b in xrange(1, nbins+1):
            morph[m].SetBinError(b, (morph[m].GetBinContent(b)*N/Neff)**0.5)

    for b in xrange(1, nbins+1):
        binMorphG[b] = TGraphErrors(len(morph_masses), array('d', [m/decimalScaling for m in morph_masses]), array('d', [morph[m].GetBinContent(b) for m in morph_masses]), array('d', [0.] * len(morph_masses)), array('d', [morph[m].GetBinError(b) for m in morph_masses]))
        binMorphG[b].SetName("%s%smorphed_bin_%d" % (name,"" if systematic == "" else systematic+"_",b))
        binMorphG[b].SetTitle("%s Morphed Bin %d" % ("t#bar{t}" if signalType == "tt" else signalType, b) )
        binMorphG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binMorphG[b].GetYaxis().SetTitle("Entries")
        binMorphG[b].GetYaxis().SetTitleOffset(1.4)

    
    #print "Done with %s  %s" % (name,systematic[1:])    
    if not morphRates:
        for m in morph_masses:
            morph[m].Scale(rates[int(decimalScaling*172.5)])


    return morph,binG,binMorphG,forward_cdfMorph,forward_cdfBinG,backward_cdfMorph,backward_cdfBinG



def morphTemplates(templates, morph_masses, name, title, precision=1, systematic="", interp = "pol3", morphRates = False, verbosity=1):
    morph = {}
    fit = {}
    fitFunc = {}
    binG = {}
    binMorphG = {}
    rates = {}
    normalized = {}
    
    signalType = name[-2:]
    
    actual_masses = sorted(templates.keys())
    if verbosity > 0: print "Morphing  %s  %s" % (name,systematic[1:])
    obs = name[4:-3]
    #print "obs =", obs
    for m in actual_masses:
        if verbosity > 1: print "Now on %s  %s  mt = %.1f" % (name,systematic,m)
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone()
        normalized[m].SetDirectory(0)
        if not morphRates:
            normalized[m].Scale(1.0 / rates[m])


    for b in xrange(1, normalized[actual_masses[0]].GetNbinsX()+1):
        binG[b] = TGraphErrors(len(actual_masses), array('d', [m/decimalScaling for m in actual_masses]), array('d', [normalized[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [normalized[m].GetBinError(b) for m in actual_masses]))
        binG[b].SetName("%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        binG[b].SetTitle("%s  Bin %d" % (obs,b))
        binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binG[b].GetYaxis().SetTitle("Entries")
        binG[b].GetYaxis().SetTitleOffset(1.3)
        binG[b].Draw()
        fit[b] = binG[b].Fit(interp, "S"+ ("" if verbosity > 1 else "Q"))
        fitFunc[b] = binG[b].GetFunction(interp)
        errors = array('d', [0.] * len(morph_masses))
        fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/decimalScaling for m in morph_masses]), errors, 2./3., False)
        for i,m in enumerate(morph_masses):
            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
            if m not in morph:
                if len(normalized[actual_masses[0]].GetXaxis().GetXbins()) == 0:
                    morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
                else:
                    morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %s GeV" % (title, massStr), normalized[actual_masses[0]].GetNbinsX(), array('d',normalized[actual_masses[0]].GetXaxis().GetXbins())) 
            morph[m].SetBinContent(b, fitFunc[b].Eval(m/decimalScaling))
            #morph[m].SetBinError(b, errors[i])
            if useNewErrors:
                #morph[m].SetBinError(b, (morph[m].GetBinContent(b)/Neff)**0.5)
                morph[m].SetBinError(b, 1000.*(morph[m].GetBinContent(b)/Neff)**0.5)
            else:
                morph[m].SetBinError(b, errors[i])
    
    for b in xrange(1, morph[morph_masses[0]].GetNbinsX()+1):
        binMorphG[b] = TGraphErrors(len(morph_masses), array('d', [m/decimalScaling for m in morph_masses]), array('d', [morph[m].GetBinContent(b) for m in morph_masses]), array('d', [0.] * len(morph_masses)), array('d', [morph[m].GetBinError(b) for m in morph_masses]))
        binMorphG[b].SetName("%s%smorphed_bin_%d" % (name,"" if systematic == "" else systematic+"_",b))
        binMorphG[b].SetTitle("%s Morphed Bin %d" % ("t#bar{t}" if signalType == "tt" else signalType, b) )
        binMorphG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binMorphG[b].GetYaxis().SetTitle("Entries")
        binMorphG[b].GetYaxis().SetTitleOffset(1.3)

    
    #print "Done with %s  %s" % (name,systematic[1:])    
    if not morphRates:
        for m in morph_masses:
            morph[m].Scale(rates[int(decimalScaling*172.5)])


    return morph,binG,binMorphG



def morphTemplates2D(templates, morph_masses, name, title, systematic="", interp = "pol3", morphRates = False, verbosity=1):
    morph = {}
    fit = {}
    fitFunc = {}
    binG = {}

    rates = {}
    normalized = {}
    actual_masses = sorted(templates.keys())
    for m in actual_masses:
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone()
        normalized[m].SetDirectory(0)
        if not morphRates:
            normalized[m].Scale(1.0 / rates[m])

    #print "morphTemplates2D"
    obs = name[4:-3]
    nbins_obs = templates[actual_masses[0]].GetNbinsX()
    nbins_mass = len(morph_masses)

    graph2D = TGraph2DErrors(nbins_obs*nbins_mass)
    for m in actual_masses:
        for b in xrange(1, nbins_obs):
            N = graph2D.GetN()+1
            graph2D.SetPoint(N, m / decimalScaling, normalized[m].GetBinCenter(b), normalized[m].GetBinContent(b))
            graph2D.SetPointError(N, 0.01, 0.01, normalized[m].GetBinError(b))
  
    graph2D.SetName("%s_Graph2D%s" % (name, systematic))
    graph2D.SetTitle("%s  %s" % (name, "" if len(systematic) == 0 else systematic[1:]))
    graph2D.GetXaxis().SetTitle("m_{t} [GeV]")
    graph2D.GetXaxis().SetTitleOffset(1.2)
    graph2D.GetYaxis().SetTitle("%s [GeV]" % obsTitle[obs])
    graph2D.GetYaxis().SetTitleOffset(1.2)
    graph2D.GetXaxis().SetRangeUser(166, 179)
    graph2D.Draw("ERR")
    #ROOT.c1.Update()
    #graph2D = TGraph2DErrors(len(morph_masses) * actual_masses[0].GetNbinsX(), array('d', [m for m in morph_masses]*nbins_obs]), array('d', [
    #print "About to fit 2D graph %s" % graph2D.GetName()
    fit = graph2D.Fit(interp, "S" + ("" if verbosity > 1 else "Q"))
    fitFunc = graph2D.FindObject(interp)
    #fitFunc = graph2D.GetFunction(interp)
    errors = array('d', [0.] * len(morph_masses))
    fit.GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/decimalScaling for m in morph_masses]), errors, 2./3., False)
    return graph2D, fit, errors


def applySmoothing(templates, masses, method, methodArgs = None, verbosity=1):
    # Apply smoothing to templates bin-by-bin
    # methodArgs will be a dictionary with method args. If None, then use defaults
    nbins = templates[masses[0]].GetNbinsX()
    name = templates[masses[0]].GetName().replace(str(masses[0]), "")

    smoothG = {}
    smoother = {}
    bins = {}
    binErr = {}
    binG = {}

    # Fill lists of bin contents & errors
    for m in masses:
        bins[m] = []
        binErr[m] = []
        for b in xrange(1,nbins+1):
            bins[m].append(templates[m].GetBinContent(b))
            binErr[m].append(templates[m].GetBinError(b))

    for _n,b in enumerate(xrange(1,nbins+1)):
        # Bin graphs for original templates
        binG[b] = TGraphErrors(len(masses), array('d', [_m/decimalScaling for _m in masses]), array('d', [bins[m][b-1] for m in masses]), array('d', [0]*len(masses)), array('d', [binErr[m][b-1] for m in masses]) )
        binG[b].SetName("%s_actual_bin_%d" % (name, b))
        binG[b].SetTitle("Actual bin %d" % b)
        binG[b].SetLineColor(kBlue)
        binG[b].SetLineWidth(2)
        binG[b].SetMarkerStyle(22)
        binG[b].SetMarkerColor(kBlue)

        smoother[b] = TGraphSmooth()

        if method.lower() == "lowess":
            if methodArgs is None:
                # Default args
                smoothG[b] = smoother[b].SmoothLowess(binG[b],'')   # Second option is not used
            else:
                if _n == 0:
                    if "span" in methodArgs:
                        span = methodArgs["span"]
                    else:
                        span = 0.67
                        if verbosity > 1:
                            print "span not found in methodArgs. Using default value %.2f" % span
                    
                    if "iter" in methodArgs:
                        it = methodArgs["iter"]
                    else:
                        it = 3
                        if verbosity > 1:
                            print "iter not found in methodArgs. Using default value %d" % it
                    
                    if "delta" in methodArgs:
                        delta = methodArgs["delta"]
                    else:
                        delta = 0.
                        if verbosity > 1:
                            print "delta not found in methodArgs. Using default value %.2f" % delta
                smoothG[b] = smoother[b].SmoothLowess(binG[b],'', span, it, delta) 
                smoothG[b].SetName("%s_smoothed_bin_%d" % (name, b))
    

    # Apply smoothing corrections to templates
    for i,m in enumerate(masses):
        for b in xrange(1,nbins+1):
            scaling = smoothG[b].GetY()[i] / templates[m].GetBinContent(b)
            templates[m].SetBinContent(b, scaling * templates[m].GetBinContent(b))
            templates[m].SetBinError(b, scaling * templates[m].GetBinError(b))

    return binG, smoothG, smoother


def create_templates(inDir, includedSysts, rebin, isToy, toyFunc, toySeed, cutMin, cutMax, massMin, massMax, deltaMT, rateScaling, observables, recoLvls, interp, outF, bins, makePlots, plotDir, debug, debugOut, includeGraphs, morphRates, useNewtoppt, useNewtW, useNewMorphing, useMorphFile, extMorphFile, useSmoothing, addBinStats, binStatsFixNorm, scaleToNominal, normalize, useAsimov, binFile, verbosity=1): 
    binning = None
    if binFile != "":
        binning = {}
        with open(binFile, "r") as _binF:
            exec("allBinning = %s" % _binF.read())
        for obs in observables:
            binning[obs] = allBinning[obs]["rec"]
        print "Using rec binning:"
        print binning

    elif bins != "":
        _bins = eval(bins)
        binning = {}
        for obs in observables:
            binning[obs] = _bins
        print "Using binning:", binning
    else:
        binning = {}
        allBinning = templateBinning
        print "\nUsing pre-defined binning:"
        for obs in observables:
            if obs in allBinning:
                binning[obs] = allBinning[obs]["rec"]
                print "%s\t" % obs, binning[obs]
        print ""

    if useNewtoppt:
        systematics["topptUp"]   = "histstoppt_up"
        systematics["topptDown"] = "histstoppt_down"
    else:
        systematics["topptUp"] = "histstoppt"
        systematics["topptDown"] = "hists"

    # If includedSysts is not an empty string, remove all other systs from systematics dictionary
    if "none" in includedSysts or "None" in includedSysts:
        # No systematics
        #systematics = {"nominal":"hists"}
        systsToRemove = deepcopy(allSystematics)
        for removeSys in systsToRemove:
            try:
                systematics.pop(removeSys+"Up")
                systematics.pop(removeSys+"Down")
            except KeyError:
                continue

        print "No systematics included!"
    elif includedSysts != "":
        systsToRemove = deepcopy(allSystematics)
        for incSys in includedSysts:
            systsToRemove.remove(incSys)

        for removeSys in systsToRemove:
            try:
                systematics.pop(removeSys+"Up")
                systematics.pop(removeSys+"Down")
            except KeyError:
                pass    

        if verbosity > 0:
            print "Only including the following systematics"
            pprint(systematics)

    if morphRates:
        print "Morphing templates to true rates"
    
    if scaleToNominal:
        print "Scaling the following systematics to the nominal rate: "
        for _s,_systs in systematicsToScale.items():
            print "%s\t%s" % (_s,_systs)

    if useMorphFile:
        print "Using morphed tt templates from %s" % args.extMorphFile
        extMorphF = TFile.Open(extMorphFile, "read")
        if morphRates:
            ttnomRate = 113799
            print "tt nominal rate hardcoded to %d" % ttnomRate



    masses = {}
    #masses["TTbar"] = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
    #masses["ST_tW"] = [1695, 1725, 1755]
    masses["TTbar"] = [int(166.5*decimalScaling), int(169.5*decimalScaling), int(171.5*decimalScaling), int(172.5*decimalScaling), int(173.5*decimalScaling), int(175.5*decimalScaling), int(178.5*decimalScaling)]
    masses["ST_tW"] = [int(169.5*decimalScaling), int(172.5*decimalScaling), int(175.5*decimalScaling)]
    #deltaM = 1
    #morph_masses = range(masses["TTbar"][0], masses["TTbar"][-1] + deltaM, deltaM)
    
    morph_masses = range(int(massMin*decimalScaling), int((massMax+deltaMT)*decimalScaling), int(deltaMT*decimalScaling)) 
    print "\nProcessing templates from %s" % inDir 

    for _m in masses["TTbar"]:
        if _m < morph_masses[0] or _m > morph_masses[-1]:
            # Remove actual mass point outside of morphing range
            masses["TTbar"].remove(_m)
            if verbosity > 0: print "Removing %s from actual TT mass range" % _m

    for _m in masses["ST_tW"]:
        if _m < morph_masses[0] or _m > morph_masses[-1]:
            # Remove actual mass point outside of morphing range
            masses["ST_tW"].remove(_m)
            if verbosity > 0: print "Removing %s from actual ST_tW mass range" % _m

    exec('minstr = "%.' + str(precision) + 'f" % (massMin)')
    exec('maxstr = "%.' + str(precision) + 'f" % (massMax)')
    
    if verbosity >= 0: print "Will produce %d morphed templates between %s and %s at %s GeV increments\n" % (len(morph_masses), minstr, maxstr, str(round(deltaMT,precision)))


    global templates
    templates = {}
    tt = {}
    tW = {}
    DY = {}
    WW = {}
    WZ = {}
    ZZ = {}
    TTW = {}
    TTZ = {}
    WJets = {}
    ST_bkgd = {}
    data_obs = {}

#    if len(obsList) == 0 or "all" in obsList:
#        # Use all observables
#        observables = _observables
#    else:
#        # Use only these observables
#        observables = obsList

    print "Observables:",
    for _obs in observables: print (" "+_obs),
    print "\n" 
    
    diffSyst = {}   # Ratio of syst/nominal @ mt=172.5 for separate sample systematics

    if binning is not None:
        useVariableBinning = True 
        variableBins = True
    else:
        useVariableBinning = None 
        variableBins = None 
    
    # Backgrounds
    for b in background:
        templates[b] = {}
        f = TFile.Open("%s/%s/%s.root" % (inDir, systematics["nominal"], b), "read")
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)  # rec_ptll, etc..
                if binning is not None:
                    tmp = f.Get("%s_%s" % (recoObs, b)).Clone("__"+b)
                    if obs in binning:
                        templates[b][recoObs] = tmp.Rebin(len(binning[obs])-1, b, array('d',binning[obs]))
                    else:
                        templates[b][recoObs] = tmp.Clone(b)
                    templates[b][recoObs].SetDirectory(0)
                    #templates[b][recoObs].Scale(1., "width")
                elif cutMin == 0 and cutMax == 0:
                    templates[b][recoObs] = f.Get("%s_%s" % (recoObs, b)).Clone(b)
                    templates[b][recoObs].SetDirectory(0)
                    templates[b][recoObs].SetTitle("%s %s %s" % (reco, obsTitle[obs], b))
                else:
                    tmp = f.Get("%s_%s" % (recoObs, b)).Clone("__"+b)
                    binW = tmp.GetBinCenter(2) - tmp.GetBinCenter(1)
                    templates[b][recoObs] = TH1F(b, tmp.GetTitle(), tmp.GetNbinsX()-cutMin-cutMax, tmp.GetBinCenter(1)+cutMin-binW/2., tmp.GetBinCenter(tmp.GetNbinsX())-cutMax+binW/2.)
                    templates[b][recoObs].SetDirectory(0)
                    for _bin in xrange(templates[b][recoObs].GetNbinsX()+1):
                        templates[b][recoObs].SetBinContent(_bin, tmp.GetBinContent(_bin+cutMin))
                        templates[b][recoObs].SetBinError(_bin, tmp.GetBinError(_bin+cutMin))


                # Rate scaling
                templates[b][recoObs].Scale(rateScaling)


                if useVariableBinning is None and binning is None:
                    # Do once on the first template to check for variable binning
                    if templates[b][recoObs].GetNbinsX() <= 1:
                        useVariableBinning = False
                    else:
                        useVariableBinning = False
                        binW = templates[b][recoObs].GetBinCenter(2) - templates[b][recoObs].GetBinCenter(1)
                        
                        # Test for variable binning
                        for _bin in xrange(2,templates[b][recoObs].GetNbinsX()):
                            nextBinW = templates[b][recoObs].GetBinCenter(_bin+1) - templates[b][recoObs].GetBinCenter(_bin)
                            if abs(binW - nextBinW) > 1e-4:                                
                                # Variable binning detected
                                useVariableBinning = True
                                # Don't rebin if variable binning:
                                #rebin=1
                                break

                    if useVariableBinning:
                        _bins = templates[b][recoObs].GetXaxis().GetXbins()
                        variableBins = [_bins[_b] for _b in xrange(len(_bins))]
                        if verbosity > 0:
                            print "Variable binning detected:"
                            print variableBins
                        
                # Rebin
   #             templates[b][recoObs].Rebin(rebin)
        f.Close()

    for s in signal:
        sample = "ttactual" if s == "TTbar" else "tWactual"
        templates[s] = {}
        for syst,systDir in systematics.iteritems():
            if verbosity > 0: print "Now loading %s: %s" % (syst, systDir) 
            if syst.find("Up") >= 0:
                systType = syst[:syst.find("Up")]
            elif syst.find("Down") >= 0:
                systType = syst[:syst.find("Down")]
            else:
                systType = syst
           
            if s == "TTbar" and systType in tWOnlySysts: continue
            if s == "ST_tW" and systType in ttOnlySysts: continue
            templates[s][syst] = {}

            if systType in separateSystSamples[s]:
                # Samples only exist for mt = 172.5
#                if s == "TTbar":
#                    print "%s is a separate sample" % syst
                fName = s 
                f = TFile.Open("%s/%s/%s.root" % (inDir, systDir, fName), "read")
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        if not recoObs in templates[s][syst]:
                            templates[s][syst][recoObs] = {}
#                        newTemplateName = ""
#                        if s == "TTbar":
#                            newTemplateName = "tt"
                        _histname = "%s%d%s" % (sample,int(decimalScaling*172.5),"" if syst=="nominal" else "_%s" % syst)
                        if binning is not None: 
                            tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__"+_histname)
                            if obs in binning:
                                templates[s][syst][recoObs][int(decimalScaling*172.5)] = tmp.Rebin(len(binning[obs])-1, _histname, array('d', binning[obs]))
                            else:
                                templates[s][syst][recoObs][int(decimalScaling*172.5)] = tmp.Clone(_histname)
                            
                            templates[s][syst][recoObs][int(decimalScaling*172.5)].SetDirectory(0)
                        elif cutMin == 0 and cutMax == 0:
                            templates[s][syst][recoObs][int(decimalScaling*172.5)] = f.Get("%s_%s" % (recoObs,fName)).Clone(_histname)
                            templates[s][syst][recoObs][int(decimalScaling*172.5)].SetDirectory(0)
                        else:
                            tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__"+_histname) 
                            binW = tmp.GetBinCenter(2) - tmp.GetBinCenter(1)
                            templates[s][syst][recoObs][int(decimalScaling*172.5)] = TH1F(_histname, tmp.GetTitle(), tmp.GetNbinsX()-cutMin-cutMax, tmp.GetBinCenter(1)+cutMin-binW/2., tmp.GetBinCenter(tmp.GetNbinsX())-cutMax+binW/2.)
                            templates[s][syst][recoObs][int(decimalScaling*172.5)].SetDirectory(0)
                            for _bin in xrange(templates[s][syst][recoObs][int(decimalScaling*172.5)].GetNbinsX()+1):
                                templates[s][syst][recoObs][int(decimalScaling*172.5)].SetBinContent(_bin, tmp.GetBinContent(_bin+cutMin))
                                templates[s][syst][recoObs][int(decimalScaling*172.5)].SetBinError(_bin, tmp.GetBinError(_bin+cutMin))
#                        templates[s][syst][recoObs][int(decimalScaling*172.5)] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,int(decimalScaling*172.5),"" if syst=="nominal" else "_%s" % syst))
#                        templates[s][syst][recoObs][int(decimalScaling*172.5)].SetDirectory(0)

                        exec('massStr = "%.' + str(precision) + 'f" % (172.5)')
                        templates[s][syst][recoObs][int(decimalScaling*172.5)].SetTitle("%s %s %s  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t}" if s == "TTbar" else "tW", massStr)) 


                        # Rate scaling
                        templates[s][syst][recoObs][int(decimalScaling*172.5)].Scale(rateScaling)

            else:
                for m in masses[s]:
                    fName = s if m == int(decimalScaling*172.5) else "%s_mt%d" % (s,int(float(m)/decimalScaling*10) )
                    f = TFile.Open("%s/%s/%s.root" % (inDir, systDir, fName), "read")
                    if verbosity > 0: print "Now loading file %s" % ("%s/%s/%s.root" % (inDir, systDir, fName))
                    for obs in observables:
                        for reco in recoLvls:
                            recoObs = "%s_%s" % (reco,obs)
                            if not recoObs in templates[s][syst]:
                                templates[s][syst][recoObs] = {}
                            
                            _histname = "%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst)
                            if binning is not None:
                                tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__"+_histname)
                                if obs in binning:
                                    templates[s][syst][recoObs][m] = tmp.Rebin(len(binning[obs])-1, _histname, array('d', binning[obs]))
                                else:
                                    templates[s][syst][recoObs][m] = tmp.Clone(_histname)
                                templates[s][syst][recoObs][m].SetDirectory(0)
                            elif cutMin == 0 and cutMax == 0:
                                templates[s][syst][recoObs][m] = f.Get("%s_%s" % (recoObs,fName)).Clone(_histname)
                                templates[s][syst][recoObs][m].SetDirectory(0)
                            else:
                                tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__" + _histname)
                                binW = tmp.GetBinCenter(2) - tmp.GetBinCenter(1)
                                templates[s][syst][recoObs][m] = TH1F(_histname, tmp.GetTitle(), tmp.GetNbinsX()-cutMin-cutMax, tmp.GetBinCenter(1)+cutMin-binW/2., tmp.GetBinCenter(tmp.GetNbinsX())-cutMax+binW/2.)
                                templates[s][syst][recoObs][m].SetDirectory(0)
                                for _bin in xrange(templates[s][syst][recoObs][m].GetNbinsX()+1):
                                    templates[s][syst][recoObs][m].SetBinContent(_bin, tmp.GetBinContent(_bin+cutMin))
                                    templates[s][syst][recoObs][m].SetBinError(_bin, tmp.GetBinError(_bin+cutMin))                            
                            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)') 
                            templates[s][syst][recoObs][m].SetTitle("%s %s %s  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t}" if s == "TTbar" else "tW", massStr)) 
                            
#                            templates[s][syst][recoObs][m] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst))
#                            templates[s][syst][recoObs][m].SetDirectory(0)
                            
                            
                            # Rate scaling
                            templates[s][syst][recoObs][m].Scale(rateScaling)
        
                    f.Close()

    # If new tW method selected:
    # nominal = 0.5 * (DR + DS)
    # up   = DS
    # down = DR
    if useNewtW and "ST_tW" in signal and "DSUp" in systematics:
        m = int(decimalScaling*172.5)
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)
                templates["ST_tW"]["DSUp"][recoObs][m].Scale(templates["ST_tW"]["nominal"][recoObs][m].Integral()/templates["ST_tW"]["DSUp"][recoObs][m].Integral())          
                DR = templates["ST_tW"]["nominal"][recoObs][m].Clone(templates["ST_tW"]["nominal"][recoObs][m].GetName()+"DR")
                DS = templates["ST_tW"]["DSUp"][recoObs][m].Clone(templates["ST_tW"]["DSUp"][recoObs][m].GetName()+"DS")
                templates["ST_tW"]["nominal"][recoObs][m].Add(DS)
                templates["ST_tW"]["nominal"][recoObs][m].Scale(0.5)



    # Create systematic variations for alternate mass points
    for s in signal:
        diffSyst[s] = {}
        sample = "ttactual" if s == "TTbar" else "tWactual"
        for syst,systDir in systematics.iteritems():
            if syst.find("Up") >= 0:
                systType = syst[:syst.find("Up")]
            elif syst.find("Down") >= 0:
                systType = syst[:syst.find("Down")]
            else:
                systType = syst

            
            if s == "TTbar" and systType in tWOnlySysts: continue
            if s == "ST_tW" and systType in ttOnlySysts: continue

            if systType in separateSystSamples[s]:
                diffSyst[s][syst] = {}
                for m in masses[s]:
                    if m == int(decimalScaling*172.5): continue

                    for obs in observables:
                        for reco in recoLvls:
                            recoObs = "%s_%s" % (reco,obs)
                        
                            diffSyst[s][syst][recoObs] = templates[s][syst][recoObs][int(decimalScaling*172.5)].Clone()
                            diffSyst[s][syst][recoObs].Divide(templates[s]['nominal'][recoObs][int(decimalScaling*172.5)])

                            templates[s][syst][recoObs][m] = systMorph(templates[s]["nominal"][recoObs][m], templates[s]["nominal"][recoObs][int(decimalScaling*172.5)].FindBin(templates[s]["nominal"][recoObs][int(decimalScaling*172.5)].GetMean()), diffSyst[s][syst][recoObs], name="%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst), useVariableBinning=useVariableBinning)


    # Rebin prior to morphing
    if rebin > 1:
        for b in background:
            for obs in observables:
                for reco in recoLvls:
                    recoObs = "%s_%s" % (reco,obs)
                    templates[b][recoObs].Rebin(rebin)
        
        for s in signal:
            for syst,systDir in systematics.iteritems():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst
                
                if s == "TTbar" and systType in tWOnlySysts: continue
                if s == "ST_tW" and systType in ttOnlySysts: continue
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        for m in masses[s]:
                            templates[s][syst][recoObs][m].Rebin(rebin)

    
    # Ref: TOP-17-001, https://indico.cern.ch/event/761804/contributions/3160985/attachments/1733339/2802398/Defranchis_template_constraints.pdf
    # If isToy == True, create toy templates by smearing templates prior to morphing
    # 
    # Each template is fluctuated according to the MC stat uncertainty
    # Systematics that come from per-event weights are fluctuated coherently with the nominal templates
    # Separate sample systematics are fluctuated independently

    global toyFluctuations,toyBinSF
    toyFluctuations = {}
    toyBinSF = {}
    _rnd = TRandom3(toySeed)

    if isToy:
        # Create data_obs here before the templates get smeared!
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)
                data_obs[recoObs] = templates["TTbar"]["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
                data_obs[recoObs].Add(templates["ST_tW"]["nominal"][recoObs][int(decimalScaling*172.5)])
                data_obs[recoObs].SetTitle("data_obs")
                for b in background:
                    # Omit WJets
                    if b != "WJets":
                        data_obs[recoObs].Add(templates[b][recoObs])

        for s in signal:
            toyFluctuations[s] = {}
            toyBinSF[s] = {}
            # Get variations used to create toys
            for syst,systDir in systematics.iteritems():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst
                
                if s == "TTbar" and systType in tWOnlySysts: continue
                if s == "ST_tW" and systType in ttOnlySysts: continue
                
                # Only calculate variations for nominal and separate sample templates
                if systType not in separateSystSamples[s] and syst != "nominal": 
                    if verbosity > 0: print "Skipping %s syst: %s" % (s,syst)
                    continue
                
                # Ignore down variations for one-sided systematics (these just use the nominal templates anyways)
                if systType in oneSidedSysts and syst.find("Down") >= 0: 
                    if verbosity > 0: print "Skipping %s down variation for syst: %s" % (s,syst)
                    continue

                toyFluctuations[s][syst] = {}
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        toyFluctuations[s][syst][recoObs] = {}
                        for m in masses[s]:
                            toyFluctuations[s][syst][recoObs][m] = {}
                            
                            # Loop over each bin
                            for b in xrange(1, templates[s][syst][recoObs][m].GetNbinsX()+1):
                                _binContent = templates[s][syst][recoObs][m].GetBinContent(b)
                                _binError   = templates[s][syst][recoObs][m].GetBinError(b)
                                if toyFunc == "gaussian":
                                    # Gaussian with mean = binContent and sigma = binError
                                    newBinContent = _rnd.Gaus(_binContent, _binError) 
                                else:
                                    # Poisson with mean = binContent
                                    newBinContent = _rnd.PoissonD(_binContent)

                                toyFluctuations[s][syst][recoObs][m][b] = newBinContent - _binContent
                                

            # Apply fluctuations to all templates
            for syst,systDir in systematics.iteritems():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst
                
                if s == "TTbar" and systType in tWOnlySysts: continue
                if s == "ST_tW" and systType in ttOnlySysts: continue
                
                toyBinSF[s][syst] = {}
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        toyBinSF[s][syst][recoObs] = {}
                        for m in masses[s]:
                            toyBinSF[s][syst][recoObs][m] = {}
                            for b in xrange(1, templates[s][syst][recoObs][m].GetNbinsX()+1):
                                _oldBinContent = templates[s][syst][recoObs][m].GetBinContent(b)
                                if _oldBinContent <= 1E-4: continue     # Skip if bin content is already zero

                                _fluc = toyFluctuations[s]["nominal"][recoObs][m][b] if (syst == "nominal" or systType not in separateSystSamples[s] or (syst.find("Down") >= 0 and systType in oneSidedSysts) ) else toyFluctuations[s][syst][recoObs][m][b] 
                                
                                # Scale factor = (oldBinContent + fluctuation)/oldBinContent = 1 + fluctuation/oldBinContent
                                _toySF = (1 + _fluc/_oldBinContent)
                                toyBinSF[s][syst][recoObs][m][b] = _toySF 

                                # Scale bin contents and errors
                                templates[s][syst][recoObs][m].SetBinContent(b, _toySF * _oldBinContent)
                                templates[s][syst][recoObs][m].SetBinError(b, _toySF * templates[s][syst][recoObs][m].GetBinError(b))




    # If scaleToNominal is selected, scale those systematics to the nominal rate
    if scaleToNominal:
        for s in signal:
            for syst,systDir in systematics.iteritems():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst
               
                if s not in systematicsToScale or systType not in systematicsToScale[s]: continue
                if syst not in templates[s]: continue
                if s == "TTbar" and systType in tWOnlySysts: continue
                if s == "ST_tW" and systType in ttOnlySysts: continue
                
                if verbosity > 0: print "Scaling %s %s" % (s, syst)
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        for m in masses[s]:
                            templates[s][syst][recoObs][m].Scale(templates[s]["nominal"][recoObs][m].Integral() / templates[s][syst][recoObs][m].Integral())



    if useSmoothing:
        # Create smoothed templates prior to morphing
        print "Applying smoothing prior to morphing"
        global originalTemplates, tt_originalBinG, tW_originalBinG, tt_smoothBinG, tW_smoothBinG, tt_smoother, tW_smoother
        # Save a separate copy of the templates pre-smoothing to originalTemplates
        # templates dict will now contain smoothed templates
        originalTemplates = {}
        tt_originalBinG = {}
        tW_originalBinG = {}
        tt_smoothBinG = {}
        tW_smoothBinG = {}
        tt_smoother = {}
        tW_smoother = {}
        for s in signal:
            originalTemplates[s] = {}
            sample = "ttsmooth" if s == "TTbar" else "tWsmooth"
            for syst, systDir in systematics.iteritems():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst

                if s == "TTbar" and systType in tWOnlySysts: continue
                if s == "ST_tW" and systType in ttOnlySysts: continue

                originalTemplates[s][syst] = {}
                if s == "TTbar":
                    tt_originalBinG[syst] = {}
                    tt_smoothBinG[syst] = {}
                    tt_smoother[syst] = {}
                else:
                    tW_originalBinG[syst] = {}
                    tW_smoothBinG[syst] = {}
                    tW_smoother[syst] = {}
                
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        originalTemplates[s][syst][recoObs] = {}
                        for m in masses[s]:
                            originalName = templates[s][syst][recoObs][m].GetName()
                            templates[s][syst][recoObs][m].SetName(originalName.replace("actual","smooth"))
                            originalTemplates[s][syst][recoObs][m] = templates[s][syst][recoObs][m].Clone(originalName)

                        # Smooth templates in place
                        if s == "TTbar":
                            tt_originalBinG[syst][recoObs], tt_smoothBinG[syst][recoObs], tt_smoother[syst][recoObs] = applySmoothing(templates[s][syst][recoObs], masses=masses[s], method="lowess", methodArgs={"span":0.9}, verbosity=verbosity)
                        else:
                            tW_originalBinG[syst][recoObs], tW_smoothBinG[syst][recoObs], tW_smoother[syst][recoObs] = applySmoothing(templates[s][syst][recoObs], masses=masses[s], method="lowess", verbosity=verbosity)

#        gROOT.SetBatch(False)
#        c2 = TCanvas("c2","c2", 1200, 1200)
#        c2.cd()
#
#        testBin = 2
#        testSyst = "nominal"
#        tW_smoothBinG[testSyst]['rec_ptll'][testBin].SetLineColor(kRed)
#        tW_smoothBinG[testSyst]['rec_ptll'][testBin].SetLineWidth(2)
#        tW_smoothBinG[testSyst]['rec_ptll'][testBin].SetMarkerStyle(23)
#        tW_smoothBinG[testSyst]['rec_ptll'][testBin].SetMarkerColor(kRed)
#        tW_originalBinG[testSyst]['rec_ptll'][testBin].SetLineWidth(2)
#        tW_originalBinG[testSyst]['rec_ptll'][testBin].SetMarkerStyle(22)
#        tW_originalBinG[testSyst]['rec_ptll'][testBin].SetMarkerColor(kBlue)
#        tW_originalBinG[testSyst]['rec_ptll'][testBin].Draw("alp")
#        tW_smoothBinG[testSyst]['rec_ptll'][testBin].Draw("lp same")
#        sys.exit()


    global tt_BinG, tt_BinMorphG,tttW_morphed, tt_morphed,tW_morphed, tt_fw_cdfMorph, tt_bw_cdfMorph, tW_fw_cdfMorph, tW_bw_cdfMorph
    # Create morphed templates for tt, tW
    tt_morphed = {}     # Morphed templates
    tt_BinG = {}        # Per-bin graphs for actual templates
    tt_BinMorphG = {}   # Per-bin graphs for morphed templates
    tt_G2D = {}         # 2D plot of morphed templates
    tt_GFit = {}
    tt_GErrors = {}
    tW_morphed = {}
    tW_BinG = {}
    tW_BinMorphG = {}
    tW_G2D = {}
    tW_GFit = {}
    tW_GErrors = {}
    tt_fw_cdfMorph = {}
    tW_fw_cdfMorph = {}
    tt_fw_cdfBinG = {}
    tW_fw_cdfBinG = {}
    tt_bw_cdfMorph = {}
    tW_bw_cdfMorph = {}
    tt_bw_cdfBinG = {}
    tW_bw_cdfBinG = {}

    tttW_morphed = {}
    tttW_G2D = {}

    print "Morphing templates...",
    sys.stdout.flush()
    for syst in systematics.keys():
        tt_morphed[syst] = {}
        tW_morphed[syst] = {}
        tt_BinG[syst] = {}
        tW_BinG[syst] = {}
        tt_BinMorphG[syst] = {}
        tW_BinMorphG[syst] = {}
        tt_G2D[syst] = {}
        tW_G2D[syst] = {}
        tt_GFit[syst] = {}
        tW_GFit[syst] = {}
        tt_GErrors[syst] = {}
        tW_GErrors[syst] = {}
        tt_fw_cdfMorph[syst] = {}
        tW_fw_cdfMorph[syst] = {}
        tt_fw_cdfBinG[syst] = {}
        tW_fw_cdfBinG[syst] = {}
        tt_bw_cdfMorph[syst] = {}
        tW_bw_cdfMorph[syst] = {}
        tt_bw_cdfBinG[syst] = {}
        tW_bw_cdfBinG[syst] = {}

        if syst.find("Up") >= 0:
            systType = syst[:syst.find("Up")]
        elif syst.find("Down") >= 0:
            systType = syst[:syst.find("Down")]
        else:
            systType = syst

        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco, obs)
                if useMorphFile:
                    tt_morphed[syst][recoObs] = {}
                    for _m in morph_masses:
                        tt_morphed[syst][recoObs][_m] = extMorphF.Get("tt%d%s" % (_m, "" if syst == "nominal" else "_"+syst))
                        tt_morphed[syst][recoObs][_m].SetDirectory(0)
                        if morphRates:
                            tt_morphed[syst][recoObs][_m].Scale(ttnomRate * ttxsec(mt=_m/10.) / 803. /tt_morphed[syst][recoObs][_m].Integral())

                        else:
                            tt_morphed[syst][recoObs][_m].Scale(1./tt_morphed[syst][recoObs][_m].Integral())
                else:
                    if useNewMorphing:
                        try:
                            #tt_morphed[syst][recoObs],tt_BinG[syst][recoObs],tt_BinMorphG[syst][recoObs], tt_fw_cdfMorph[syst][recoObs], tt_fw_cdfBinG[syst][recoObs] = cdfMorphTemplates(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp=interp, morphRates=morphRates, verbosity=verbosity)
                            tt_morphed[syst][recoObs],tt_BinG[syst][recoObs],tt_BinMorphG[syst][recoObs],tt_fw_cdfMorph[syst][recoObs], tt_fw_cdfBinG[syst][recoObs],tt_bw_cdfMorph[syst][recoObs], tt_bw_cdfBinG[syst][recoObs] = avgCdfMorphTemplates(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp=interp, morphRates=morphRates, verbosity=verbosity)
                        except KeyError:
                            pass
                    else:
                        try:
                            tt_morphed[syst][recoObs],tt_BinG[syst][recoObs],tt_BinMorphG[syst][recoObs] = morphTemplates(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp=interp, morphRates=morphRates, verbosity=verbosity)
                        except KeyError:
                            pass

                # tW morphing
                if useNewMorphing:
                    if systType not in ttOnlySysts:
                        #tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs],tW_fw_cdfMorph[syst][recoObs], tW_fw_cdfBinG[syst][recoObs] = cdfMorphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", morphRates=morphRates, verbosity=verbosity)
                        tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs],tW_fw_cdfMorph[syst][recoObs], tW_fw_cdfBinG[syst][recoObs],tW_bw_cdfMorph[syst][recoObs], tW_bw_cdfBinG[syst][recoObs] = avgCdfMorphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", morphRates=morphRates, verbosity=verbosity)
                        
#                    try:
#                        tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs],tW_cdfMorph[syst][recoObs] = cdfMorphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", morphRates=morphRates, verbosity=verbosity)
#                    except KeyError:        
#                        pass 
                else:
                    try:
                        tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs] = morphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", morphRates=morphRates, verbosity=verbosity)
                    except KeyError:
                        pass

#                if not useNewMorphing:
#                    try:
#                        tt_G2D[syst][recoObs],tt_GFit[syst][recoObs],tt_GErrors[syst][recoObs] = morphTemplates2D(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp=interp, morphRates=morphRates, verbosity=verbosity)
#                    except KeyError:
#                        pass
#                    
#                    try:
#                        tW_G2D[syst][recoObs],tW_GFit[syst][recoObs],tW_GErrors[syst][recoObs] = morphTemplates2D(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", morphRates=morphRates, verbosity=verbosity)
#                    except KeyError:
#                        pass

    print "done"

    # Create sum tt+tW templates
    if "TTbar" in signal and "ST_tW" in signal:
        print "Summing tt and tW templates...",
        sys.stdout.flush()
        templates["tt+tW"] = {}
        if useSmoothing:
            originalTemplates["tt+tW"] = {}
        for syst,systDir in systematics.iteritems():
            
            if syst.find("Up") >= 0:
                systType = syst[:syst.find("Up")]
            elif syst.find("Down") >= 0:
                systType = syst[:syst.find("Down")]
            else:
                systType = syst
           
            templates["tt+tW"][syst] = {}
            if useSmoothing:
                originalTemplates["tt+tW"][syst] = {}
            tttW_morphed[syst] = {}

            for obs in observables:
                for reco in recoLvls:
                    recoObs = "%s_%s" % (reco,obs)
                    if not recoObs in templates["tt+tW"][syst]:
                        templates["tt+tW"][syst][recoObs] = {}
                        if useSmoothing:
                            originalTemplates["tt+tW"][syst][recoObs] = {}
                    if not recoObs in tttW_morphed[syst]:
                        tttW_morphed[syst][recoObs] = {}

                    # Actual tt+tW templates
                    if useSmoothing:
                        for m in masses["ST_tW"]:
                            if systType not in tWOnlySysts:
                                templates["tt+tW"][syst][recoObs][m] = templates["TTbar"][syst][recoObs][m].Clone(templates["TTbar"][syst][recoObs][m].GetName().replace("ttsmooth","tttWsmooth") )
                                if useSmoothing:
                                    originalTemplates["tt+tW"][syst][recoObs][m] = originalTemplates["TTbar"][syst][recoObs][m].Clone(originalTemplates["TTbar"][syst][recoObs][m].GetName().replace("ttactual","tttWactual") )
                                if systType not in ttOnlySysts:
                                    templates["tt+tW"][syst][recoObs][m].Add(templates["ST_tW"][syst][recoObs][m])
                                    if useSmoothing:
                                        originalTemplates["tt+tW"][syst][recoObs][m].Add(originalTemplates["ST_tW"][syst][recoObs][m])
                                else:
                                    # Add nominal tW
                                    templates["tt+tW"][syst][recoObs][m].Add(templates["ST_tW"]["nominal"][recoObs][m])
                                    if useSmoothing:
                                        originalTemplates["tt+tW"][syst][recoObs][m].Add(originalTemplates["ST_tW"]["nominal"][recoObs][m])

                            else:
                                templates["tt+tW"][syst][recoObs][m] = templates["ST_tW"][syst][recoObs][m].Clone(templates["ST_tW"][syst][recoObs][m].GetName().replace("tWsmooth","tttWsmooth") )
                                # Use nominal tt 
                                templates["tt+tW"][syst][recoObs][m].Add(templates["TTbar"]["nominal"][recoObs][m])
                                if useSmoothing:
                                    originalTemplates["tt+tW"][syst][recoObs][m] = originalTemplates["ST_tW"][syst][recoObs][m].Clone(originalTemplates["ST_tW"][syst][recoObs][m].GetName().replace("tWactual","tttWactual") )
                                    originalTemplates["tt+tW"][syst][recoObs][m].Add(originalTemplates["TTbar"]["nominal"][recoObs][m])

                            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
                            templates["tt+tW"][syst][recoObs][m].SetTitle("%s %s %s smoothed  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t} + tW", massStr))
                            if useSmoothing:
                                originalTemplates["tt+tW"][syst][recoObs][m].SetTitle("%s %s %s  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t} + tW", massStr))
                    else:
                        for m in masses["ST_tW"]:
                            if systType not in tWOnlySysts:
                                templates["tt+tW"][syst][recoObs][m] = templates["TTbar"][syst][recoObs][m].Clone(templates["TTbar"][syst][recoObs][m].GetName().replace("ttactual","tttWactual") )
                                if systType not in ttOnlySysts:
                                    templates["tt+tW"][syst][recoObs][m].Add(templates["ST_tW"][syst][recoObs][m])
                                else:
                                    # Add nominal tW
                                    templates["tt+tW"][syst][recoObs][m].Add(templates["ST_tW"]["nominal"][recoObs][m])

                            else:
                                templates["tt+tW"][syst][recoObs][m] = templates["ST_tW"][syst][recoObs][m].Clone(templates["ST_tW"][syst][recoObs][m].GetName().replace("tWactual","tttWactual") )
                                # Use nominal tt 
                                templates["tt+tW"][syst][recoObs][m].Add(templates["TTbar"]["nominal"][recoObs][m])

                            exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
                            templates["tt+tW"][syst][recoObs][m].SetTitle("%s %s %s  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t} + tW", massStr))

                    # Morphed tt+tW templates (summed after morphing)
                    for m in morph_masses:
                        if systType not in tWOnlySysts:
                            tttW_morphed[syst][recoObs][m] = tt_morphed[syst][recoObs][m].Clone(tt_morphed[syst][recoObs][m].GetName().replace("tt","tttW"))
                            if systType not in ttOnlySysts:
                                tttW_morphed[syst][recoObs][m].Add(tW_morphed[syst][recoObs][m])
                            else:
                                tttW_morphed[syst][recoObs][m].Add(tW_morphed["nominal"][recoObs][m])
                        else:
                            tttW_morphed[syst][recoObs][m] = tW_morphed[syst][recoObs][m].Clone(tW_morphed[syst][recoObs][m].GetName().replace("tW","tttW"))
                            tttW_morphed[syst][recoObs][m].Add(tt_morphed["nominal"][recoObs][m])

                        exec('massStr = "%.' + str(precision) + 'f" % (m/decimalScaling)')
                        tttW_morphed[syst][recoObs][m].SetTitle("%s %s %s  m_{t} = %s GeV" % (reco, obsTitle[obs], "t#bar{t} + tW", massStr))
    
        print "done"


    if addBinStats and "TTbar" in signal and "ST_tW" in signal:
        # Add an additional nuisance parameter for each bin
        print "Adding binwise stat uncertainty nuisance parameters for each signal bin %s" % ("with normalization fixed" if binStatsFixNorm else "without fixing overall normalization")
        nSigBins = {}
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)
                
                nSigBins[recoObs] = tttW_morphed["nominal"][recoObs][int(decimalScaling*172.5)].GetNbinsX()
                
                for _b in xrange(1,nSigBins[recoObs]+1):
                    if "%s_bin%dUp" % (recoObs,_b) not in tttW_morphed:
                        tttW_morphed["%s_bin%dUp" % (recoObs,_b)] = {recoObs:{} }
                        tttW_morphed["%s_bin%dDown" % (recoObs,_b)] = {recoObs:{} }

                    for m in morph_masses:
                        tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m] = tttW_morphed["nominal"][recoObs][m].Clone(tttW_morphed["nominal"][recoObs][m].GetName() + "_bin%dUp" % _b)
                        tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m] = tttW_morphed["nominal"][recoObs][m].Clone(tttW_morphed["nominal"][recoObs][m].GetName() + "_bin%dDown" % _b)
                        
                        # Bin up = Nominal + bin_error
                        tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m].SetBinContent(_b, tttW_morphed["nominal"][recoObs][m].GetBinContent(_b) + tttW_morphed["nominal"][recoObs][m].GetBinError(_b))
                        if binStatsFixNorm:
                            tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m].Scale(tttW_morphed["nominal"][recoObs][m].Integral() / tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m].Integral())

                        # Bin down = Nominal - bin_error
                        tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m].SetBinContent(_b, tttW_morphed["nominal"][recoObs][m].GetBinContent(_b) - tttW_morphed["nominal"][recoObs][m].GetBinError(_b))
                        if binStatsFixNorm:
                            tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m].Scale(tttW_morphed["nominal"][recoObs][m].Integral() / tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m].Integral())


    if not isToy:
        # Fine to create data_obs distribution here if not in toy mode
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)
#            data_obs[recoObs] = tttW_morphed["nominal"][recoObs][int(decimalScaling*175.5)].Clone("data_obs")
                if useAsimov:
                    data_obs[recoObs] = tttW_morphed["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
                else:
                    data_obs[recoObs] = templates["TTbar"]["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
                    data_obs[recoObs].Add(templates["ST_tW"]["nominal"][recoObs][int(decimalScaling*172.5)])
                #data_obs[recoObs] = tt_morphed["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
                data_obs[recoObs].SetTitle("data_obs")
                #data_obs[recoObs].Add(templates["ST_tW"]["nominal"][recoObs][int(decimalScaling*172.5)])
#            data_obs[recoObs].Add(tW_morphed["nominal"][recoObs][int(decimalScaling*172.5)])
                for b in background:
                    # Omit WJets
                    if b != "WJets":
                        data_obs[recoObs].Add(templates[b][recoObs])


    if normalize:
        # Normalize disributions to unity
        print "Normalizing morphed templates to nominal rate"
        
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)
                
                for m in morph_masses:
                    tttW_morphed[syst][recoObs][m].Scale(tttW_morphed["nominal"][recoObs][int(decimalScaling*175.5)].Integral()/tttW_morphed["nominal"][recoObs][m].Integral()) 
                
#                for syst,systDir in systematics.iteritems():
#                    
#                    if syst.find("Up") >= 0:
#                        systType = syst[:syst.find("Up")]
#                    elif syst.find("Down") >= 0:
#                        systType = syst[:syst.find("Down")]
#                    else:
#                        systType = syst
#                    
#                    for m in morph_masses:
#                        tttW_morphed[syst][recoObs][m].Scale(tttW_morphed["nominal"][recoObs][m].Integral()/tttW_morphed[syst][recoObs][m].Integral()) 

    outFile = TFile.Open(outF, "recreate")
    for obs in observables:
        for reco in recoLvls:
            recoObs = "%s_%s" % (reco,obs)
            outFile.mkdir(recoObs)
            outFile.cd(recoObs)
            #outFile.mkdir("ttbins")
            #outFile.mkdir("tWbins")
            data_obs[recoObs].Write()

            if addBinStats:
                for _b in xrange(1,nSigBins[recoObs]+1):
                    for m in morph_masses:
                        tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m].Write(tttW_morphed["%s_bin%dUp" % (recoObs,_b)][recoObs][m].GetName()[len(recoObs)+1:] )
                        tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m].Write(tttW_morphed["%s_bin%dDown" % (recoObs,_b)][recoObs][m].GetName()[len(recoObs)+1:] )

            for syst in systematics.keys():
                if syst.find("Up") >= 0:
                    systType = syst[:syst.find("Up")]
                elif syst.find("Down") >= 0:
                    systType = syst[:syst.find("Down")]
                else:
                    systType = syst
                
                
                if includeGraphs:
                    try:
                        for b in xrange(1,len(tt_BinG[syst][recoObs])):
                            tt_BinG[syst][recoObs][b].Write(tt_BinG[syst][recoObs][b].GetName()[len(recoObs)+1:])
                            tt_BinMorphG[syst][recoObs][b].Write(tt_BinMorphG[syst][recoObs][b].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass

                    try:
                        for b in xrange(1,len(tW_BinG[syst][recoObs])):
                            tW_BinG[syst][recoObs][b].Write(tW_BinG[syst][recoObs][b].GetName()[len(recoObs)+1:])
                            tW_BinMorphG[syst][recoObs][b].Write(tW_BinMorphG[syst][recoObs][b].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass

                
                for m in morph_masses:
                    if PU_mt1725only and "pileup" in syst and m != int(decimalScaling*172.5): continue
                
                    # tt+tW templates
                    tttW_morphed[syst][recoObs][m].Write(tttW_morphed[syst][recoObs][m].GetName()[len(recoObs)+1:])
                    try:
                        tt_morphed[syst][recoObs][m].Write(tt_morphed[syst][recoObs][m].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass

                    try:
                        tW_morphed[syst][recoObs][m].Write(tW_morphed[syst][recoObs][m].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass
    
                if "TTbar" in signal and "ST_tW" in signal:
                    for m in masses["ST_tW"]:
                        templates["tt+tW"][syst][recoObs][m].Write()
                        if useSmoothing:
                            originalTemplates["tt+tW"][syst][recoObs][m].Write()


                for s in signal:
                    if s == "TTbar" and systType in tWOnlySysts: continue
                    if s == "ST_tW" and systType in ttOnlySysts: continue
                    for m in masses[s]:
                        templates[s][syst][recoObs][m].Write()
                        if useSmoothing:
                            originalTemplates[s][syst][recoObs][m].Write()
                            
#                        try:
#                            templates[s][syst][recoObs][m].Write()
#                        except KeyError:
#                            print "Error saving template %s %s %s %d" % (s, syst, recoObs, m)
#                            pass
#                        #templates[s][syst][recoObs][m].Write(templates[s][syst][recoObs][m].GetName()[len(recoObs):])
                


                if includeGraphs:
                    try:
                        tt_G2D[syst][recoObs].Write(tt_G2D[syst][recoObs].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass
                        
                    try:
                        tW_G2D[syst][recoObs].Write(tW_G2D[syst][recoObs].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass

            for b in background:
                templates[b][recoObs].Write()
                #templates[b][recoObs].Write(templates[b][recoObs].GetName()[len(recoObs):])
            #outFile.cd()
            #outFile.cd(recoObs)
            #gDirectory.cd(recoObs)
    
    outFile.Close()
    print "Output templates saved to %s" % outF

    if makePlots:
        if verbosity > 0: print "Making per-bin plots"
        c = TCanvas("c","c", 800, 600)
        # Make individual bin plots
        for syst,systDir in systematics.iteritems():
            for obs in observables:
                for reco in recoLvls:
                    recoObs = "%s_%s" % (reco,obs) 

                    try:
                        for b in xrange(1, tt_morphed[syst][recoObs][morph_masses[0]].GetNbinsX()+1):
                            l = TLegend(0.425, 0.8, 0.575, 0.9)
                            l.SetFillStyle(0)
                            l.SetBorderSize(0)
                            tt_BinMorphG[syst][recoObs][b].SetFillColor(0)
                            tt_BinG[syst][recoObs][b].SetFillColor(0)
                            tt_BinG[syst][recoObs][b].SetLineWidth(2)
                            #tt_BinMorphG[syst][recoObs][b].Draw("ALP")
                            if useSmoothing:
                                tt_originalBinG[syst][recoObs][b].SetLineColor(kBlue)
                                tt_originalBinG[syst][recoObs][b].SetFillColor(0)
                                tt_originalBinG[syst][recoObs][b].SetLineWidth(2)
                                #tt_originalBinG[syst][recoObs][b].Draw("LP SAME")
                                tt_originalBinG[syst][recoObs][b].Draw("ALP")
                                tt_originalBinG[syst][recoObs][b].SetTitle("t#bar{t} Morphed Bin %d" % b)
                                tt_originalBinG[syst][recoObs][b].GetXaxis().SetTitle(tt_BinG[syst][recoObs][b].GetXaxis().GetTitle())
                                tt_originalBinG[syst][recoObs][b].GetXaxis().SetTitleOffset(tt_BinG[syst][recoObs][b].GetXaxis().GetTitleOffset())
                                tt_originalBinG[syst][recoObs][b].GetYaxis().SetTitle(tt_BinG[syst][recoObs][b].GetYaxis().GetTitle())
                                tt_originalBinG[syst][recoObs][b].GetYaxis().SetTitleOffset(tt_BinG[syst][recoObs][b].GetYaxis().GetTitleOffset())

                                tt_BinMorphG[syst][recoObs][b].Draw("LP SAME")
                                tt_BinG[syst][recoObs][b].SetLineColor(kGreen+1)
                                tt_BinG[syst][recoObs][b].SetMarkerStyle(23)
                                tt_BinG[syst][recoObs][b].SetMarkerColor(kGreen+1)
                                tt_BinG[syst][recoObs][b].Draw("LP SAME")
                                l.AddEntry(tt_originalBinG[syst][recoObs][b], "Actual")
                                l.AddEntry(tt_BinG[syst][recoObs][b], "Smoothed")
                                l.AddEntry(tt_BinMorphG[syst][recoObs][b], "Morphed")
                            else:
                                tt_BinG[syst][recoObs][b].SetLineColor(kBlue)
                                #tt_BinG[syst][recoObs][b].Draw("LP SAME")
                                tt_BinG[syst][recoObs][b].Draw("ALP")
                                tt_BinMorphG[syst][recoObs][b].Draw("LP SAME")
                                l.AddEntry(tt_BinG[syst][recoObs][b], "Actual")
                                l.AddEntry(tt_BinMorphG[syst][recoObs][b], "Morphed")
                            fit = tt_BinG[syst][recoObs][b].Fit(interp, "S" if verbosity > 1 else "SQ")
                            l.Draw("SAME")
                            os.system("mkdir -p %s/%s/%s" % (plotDir,recoObs,syst))
                            c.SaveAs("%s/%s/%s/tt_morph_bin_%d.png" % (plotDir,recoObs,syst,b))
                    except KeyError:
                        pass

                    try:
                        for b in xrange(1, tW_morphed[syst][recoObs][morph_masses[0]].GetNbinsX()+1):
                            l = TLegend(0.425, 0.8, 0.575, 0.9)
                            l.SetFillStyle(0)
                            l.SetBorderSize(0)
                            tW_BinMorphG[syst][recoObs][b].SetFillColor(0)
                            tW_BinG[syst][recoObs][b].SetFillColor(0)
                            tW_BinG[syst][recoObs][b].SetLineWidth(2)
                            tW_BinMorphG[syst][recoObs][b].Draw("ALP")
                            if useSmoothing:
                                tW_originalBinG[syst][recoObs][b].SetLineColor(kBlue)
                                tW_originalBinG[syst][recoObs][b].SetFillColor(0)
                                tW_originalBinG[syst][recoObs][b].SetLineWidth(2)
                                tW_originalBinG[syst][recoObs][b].Draw("LP SAME")
                                tW_BinG[syst][recoObs][b].SetLineColor(kGreen+1)
                                tW_BinG[syst][recoObs][b].SetMarkerStyle(23)
                                tW_BinG[syst][recoObs][b].SetMarkerColor(kGreen+1)
                                tW_BinG[syst][recoObs][b].Draw("LP SAME")
                                l.AddEntry(tW_originalBinG[syst][recoObs][b], "Actual")
                                l.AddEntry(tW_BinG[syst][recoObs][b], "Smoothed")
                                l.AddEntry(tW_BinMorphG[syst][recoObs][b], "Morphed")
                            else:
                                tW_BinG[syst][recoObs][b].SetLineColor(kBlue)
                                tW_BinG[syst][recoObs][b].Draw("LP SAME")
                                l.AddEntry(tW_BinG[syst][recoObs][b], "Actual")
                                l.AddEntry(tW_BinMorphG[syst][recoObs][b], "Morphed")
                            fit = tW_BinG[syst][recoObs][b].Fit("pol1", "S" if verbosity > 1 else "SQ")
                            l.Draw("SAME")
                            os.system("mkdir -p %s/%s/%s" % (plotDir,recoObs,syst))
                            c.SaveAs("%s/%s/%s/tW_morph_bin_%d.png" % (plotDir,recoObs,syst,b))
                    except KeyError:
                        pass

        print "Plots saved to %s" % plotDir

    if debug:
        with gzip.open(debugOut, "wb") as f:
            pickle.dump(templates, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tt_G2D, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tt_GFit, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tt_GErrors, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tW_G2D, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tW_GFit, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(tW_GErrors, f, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(diffSyst, f, protocol=pickle.HIGHEST_PROTOCOL)
        print "Debugging info saved to %s\n" % debugOut

    sys.stdout.flush()
    #return templates,tt_morphed,tW_morphed,tt_G2D,tt_GFit,tt_GErrors,diffSyst
    print "Now clearing dictionaries"
    templates.clear()
    tt_morphed.clear()
    tW_morphed.clear()
    tt_G2D.clear()
    tt_GFit.clear()
    tt_GErrors.clear()
    diffSyst.clear()
    print "Done clearing."
    return
    

def main():
    parser = ArgumentParser()
    parser.add_argument("--topDir", default="/uscms_data/d3/msaunder/TopMass/CMSSW_8_0_26_patch1/src/TopNtuplizer/Plotting", help="path to Plotting directory")
    parser.add_argument("-i", dest="inDir", default="histograms", help="Input directory containing ttrees")
    parser.add_argument("-o", dest="outF", default="mtTemplatesForCH.root", help="Output template file")
    parser.add_argument("--debug", action="store_true", default=False, help="store debugging info")
    parser.add_argument("--debugOut", default="", help="output file to store pickled template info")
    parser.add_argument("-b", "--rebin", type=int, default=1, help="Integer rebin width in GeV")
    parser.add_argument("--systs", default="", nargs="*", choices=(allSystematics + ["none","None"]), help="ONLY plot these systematics")
    parser.add_argument("--toy", action="store_true", default=False, help="create toy templates by fluctuating bin contents according to MC uncertainty")
    parser.add_argument("--toyFunc", default="gaussian", choices=["poisson", "gaussian"], help="distribution to use for fluctuating bin contents when creating toys")
    parser.add_argument("--toySeed", type=int, default=0, help="seed for TRandom3 to calculate toy fluctuations")
    parser.add_argument("--cutMin", type=int, default=0, help="bins to cut from the left in GeV (before rebinning)")
    parser.add_argument("--cutMax", type=int, default=0, help="bins to cut from the right in GeV (before rebinning)")
    parser.add_argument("--minmt", type=float, default=166.5, help="minimum mass for morphing range")
    parser.add_argument("--maxmt", type=float, default=178.5, help="maximum mass for morphing range")
    parser.add_argument("--deltaMT", type=float, default=0.1, help="morphing mass increment (in GeV)") 
    parser.add_argument("--normalize", action="store_true", default=False, help="normalize to unity") 
    parser.add_argument("--newtoppt", action="store_true", default=False, help="use top pT rw as nominal wth 2 sided toppt systematic") 
    parser.add_argument("--newtW", action="store_true", default=False, help="use ST tW 0.5(DR+DS) as nominal, DR down, DS up") 
    parser.add_argument("--noMorphRates", action="store_true", default=False, help="interpolate non-normalized histograms")
    parser.add_argument("--noScalingToNominal", action="store_true", default=False, help="don't scale certain systematics to the nominal rate")
    parser.add_argument("--useNewMorphing", action="store_true", default=False, help="use new morphing method")
    parser.add_argument("--useMorphFile", action="store_true", default=False, help="use morphed templates from external file instead of doing morphing here")
    parser.add_argument("--extMorphFile", default="/uscms/homes/m/msaunder/work/DNN_TemplateMorphing/new/newplots/morphTF.root", help="external morph template file loaded when useMorphFile option is selected")
    parser.add_argument("-r", "--rateScaling", type=float, default=1., help="scale rates by a constant factor")
#    parser.add_argument("--bins", type=str, default="[0, 15, 30, 45, 60, 70, 80, 90, 100, 110, 120, 130, 140]", help="List of variable bin ranges. The last entry is the upper edge of the last bin. All other entries are the lower bin edges")
    parser.add_argument("--smooth", action="store_true", default=False, help="apply template smoothing")
    parser.add_argument("--binF", "--binFile", dest="binFile", default="", help="file of bins for rec/gen observables")
    parser.add_argument("--bins", type=str, default="", help="List of variable bin ranges. The last entry is the upper edge of the last bin. All other entries are the lower bin edges")
    parser.add_argument("--noBinStats", action="store_true", default=False, help="don't add binwise stat unc nps")
    parser.add_argument("--binStatsFixNorm", action="store_true", default=False, help="when computing binwise stats, adjust other bin contents to preserve the overall normalization")
    parser.add_argument("--newErrors", action="store_true", default=False, help="use new error method")
    parser.add_argument("--precision", type=int, default=1, help="number of decimal places to use for morphing")
    parser.add_argument("--obs", nargs="+", default=['ptll'], choices=(_observables + ["kin","all","diff"]), help="include these observables (or 'all')")
    parser.add_argument("--reco", nargs="+", default=["rec"], choices=["rec","gen"], help="reco lvl")
    parser.add_argument("--interp", default="pol3", help="ttbar interpolation function to use in ROOT Fit")
    parser.add_argument("--plots", action = "store_true", default=False, help="create bin plots")
    parser.add_argument("--plotDir", default="", help="directory to store bin plots if --plots is selected")
    parser.add_argument("--includeGraphs", action="store_true", default=False, help="Store per-bin graphs in output root file")
    parser.add_argument("-a", "--asimov", action="store_true", default=False, help="Use asimov set at nominal mass for data_obs")
    parser.add_argument("-v", "-V", "--verbosity", dest="verbosity", type=int, default=0, help="verbosity of output")
    args = parser.parse_args()

    isToy = args.toy
    toyFunc = args.toyFunc
    toySeed = args.toySeed

    print ""
    if isToy:
        print "Generating toy templates using %s smearing and seed = %d%s" % (toyFunc, toySeed, " (randomized)" if toySeed == 0 else "")

    newtoppt = args.newtoppt 
    newtW = args.newtW
    if newtoppt: 
        print "Using top pT reweighted samples for all ttbar masses"
    else:
        print "Using one-sided top pT reweighting systematic"
        oneSidedSysts.append("toppt")

    if newtW:
        print "Using 0.5(DR + DR) as nominal ST tW sample"
    else:
        print "Using ST tW DR as nominal sample with DS as one-sided sytematic"
        oneSidedSysts.append("DS")

    useAsimov = args.asimov

    if useAsimov:
        print "Using asimov dataset at mt=172.5 for data_obs"

    scaleToNominal = not args.noScalingToNominal
    if "all" in args.obs:
        # Include all observables
        print "All observables selected"
        args.obs = _observables
    elif "kin" in args.obs:
        print "Kinematic observables selected"
        args.obs = ["ptll","Mll","Epos","Eneg","ptpos","ptneg","Ep_Em","ptp_ptm"]
    elif "diff" in args.obs:
        print "Differential observables selected"
        args.obs = _diffDists

    if args.useNewMorphing:
        print "Using new morphing method"

    useSmoothing = args.smooth
    if useSmoothing:
        print "Template smoothing will be applied"

    if args.precision < 0:
        print "Cannot have %d decimal places! Defaulting to 1" % args.precision
        args.precision = 1

    elif args.deltaMT < 10**(-args.precision):
        # Increase precision to leading decimal value in deltaMT
        args.precision = int(TMath.Ceil(TMath.Log10(1./args.deltaMT)))
        print "Using precision: %d" % args.precision 

    if args.debugOut == "":
        if "mtTemplatesForCH.root" in args.outF:
            args.debugOut = args.outF.replace("mtTemplatesForCH.root", "debugTemplates.pklz")
        elif ".root" in args.outF:
            args.debugOut = args.outF.replace(".root", "debugTemplates.pklz")
        else:
            args.debugOut = "debugTemplates.pklz"


    if args.plots and args.plotDir == "":
        # Directory for bin plots
        if args.outF.find("mtTemplatesForCH.root") >= 0:
            args.plotDir = args.outF.replace("mtTemplatesForCH.root", "binPlots")
        elif args.outF[-5:] == ".root":
            args.plotDir = args.outF.replace(".root", "_binPlots")
        else:
            args.plotDir = "binPlots"

    if args.binFile != "":
        print "Using binning from %s" % args.binFile

    morphRates = not args.noMorphRates
    if args.noMorphRates:
        print "Templates will be morphed using normalized templates then scaled to nominal 172.5 template rate"


    global precision,decimalScaling
    precision = args.precision          # Number of decimal places 
    decimalScaling = 10**precision     # Scaling such that the specified precision can be represented by integers
    
    global useNewErrors,Neff
    useNewErrors = args.newErrors
    Neff = 329857.

    if args.verbosity <= 1:
        ROOT.gErrorIgnoreLevel = kWarning


    if args.topDir[-1] == "/": args.topDir = args.topDir[:-1]
    if args.inDir[-1] == "/": args.inDir = args.inDir[:-1]
    if args.rebin < 1:
        print "Invalid rebin value. Must be >= 1! Defaulting to 2 GeV"
        args.rebin = 2
    else:
        if args.verbosity > 0: print "Rebinning to %d GeV" % args.rebin

    
    # Create a set of templates for the given ttres and config 
    create_templates(inDir="%s/%s" % (args.topDir,args.inDir), includedSysts=args.systs, isToy=isToy, toyFunc=toyFunc, toySeed=toySeed, cutMin=args.cutMin, cutMax=args.cutMax, massMin=args.minmt, massMax=args.maxmt, deltaMT=args.deltaMT, rateScaling=args.rateScaling, observables=args.obs, recoLvls=args.reco, rebin=args.rebin, bins=args.bins, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir, debug=args.debug, debugOut=args.debugOut, includeGraphs=args.includeGraphs, morphRates=morphRates, useNewtoppt=newtoppt, useNewtW=newtW, useNewMorphing=args.useNewMorphing, useMorphFile=args.useMorphFile, extMorphFile=args.extMorphFile, useSmoothing=useSmoothing, addBinStats=(not args.noBinStats), binStatsFixNorm=args.binStatsFixNorm, scaleToNominal=scaleToNominal,normalize=args.normalize,useAsimov=useAsimov,binFile=args.binFile,verbosity=args.verbosity)

    #return

    print "Finished"
    return


if __name__ == "__main__":
    sys.exit(main())

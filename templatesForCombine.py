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

gROOT.SetBatch(True)
allSystematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS',"MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"]

systematicsToScale = {"TTbar":["Q2","MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"],
                      "ST_tW":["Q2"]}
_observables = ['ptll', 'Mll', 'ptpos', 'Epos', 'ptp_ptm', 'Ep_Em', "leadLepPt", "leadJetPt"] 
obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)", "leadLepPt":"Leading Lepton p_{T}", "leadJetPt":"Leading Jet p_{T}"}
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
ttOnlySysts = ['toppt','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp',"MEscale1","MEscale2","MEscale3","MEscale4","MEscale5","MEscale6"]
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


def morphTemplates(templates, morph_masses, name, title, precision=1, systematic="", variableBins=None, interp = "pol3", morphRates = False, verbosity=1):
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
                if variableBins is None:
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



def morphTemplates2D(templates, morph_masses, name, title, systematic="", variableBins = None, interp = "pol3", morphRates = False, verbosity=1):
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



def create_templates(inDir, includedSysts, rebin, cutMin, cutMax, massMin, massMax, deltaMT, rateScaling, observables, recoLvls, interp, outF, bins, makePlots, plotDir, debugOut, includeGraphs, morphRates, useMorphFile, extMorphFile, scaleToNominal, verbosity=1): 
    binning = None
    if bins != "":
        binning = eval(bins)
        print "Using binning:", binning
    # If includedSysts is not an empty string, remove all other systs from systematics dictionary
    if "none" in includedSysts or "None" in includedSysts:
        # No systematics
        #systematics = {"nominal":"hists"}
        systsToRemove = deepcopy(allSystematics)
        for removeSys in systsToRemove:
            systematics.pop(removeSys+"Up")
            systematics.pop(removeSys+"Down")
            
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
    
    if verbosity > 0: print "Will produce %d morphed templates between %s and %s at %s GeV increments\n" % (len(morph_masses), minstr, maxstr, str(round(deltaMT,precision)))


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

#    if len(obsList) == 0 or "all" in obsList:
#        # Use all observables
#        observables = _observables
#    else:
#        # Use only these observables
#        observables = obsList

    print "Observables:",
    for _obs in observables: print (" "+_obs),
    print "" 
    
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
                    templates[b][recoObs] = tmp.Rebin(len(binning)-1, b, array('d',binning))
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
                                rebin=1
                                break

                    if useVariableBinning:
                        _bins = templates[b][recoObs].GetXaxis().GetXbins()
                        variableBins = [_bins[_b] for _b in xrange(len(_bins))]
                        if verbosity > 0:
                            print "Variable binning detected:"
                            print variableBins
                        
#                if useVariableBinning and binning is None:
#                    templates[b][recoObs].Scale(1., "width")

                # Rebin
   #             templates[b][recoObs].Rebin(rebin)
        f.Close()

    for s in signal:
        sample = "ttactual" if s == "TTbar" else "tWactual"
        templates[s] = {}
        for syst,systDir in systematics.iteritems():
            
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
                            templates[s][syst][recoObs][int(decimalScaling*172.5)] = tmp.Rebin(len(binning)-1, _histname, array('d', binning))
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
#                        if useVariableBinning and binning is None:
#                            templates[s][syst][recoObs][int(decimalScaling*172.5)].Scale(1., "width")


                        # Rate scaling
                        templates[s][syst][recoObs][int(decimalScaling*172.5)].Scale(rateScaling)

            else:
                for m in masses[s]:
                    fName = s if m == int(decimalScaling*172.5) else "%s_mt%d" % (s,int(float(m)/decimalScaling*10) )
                    f = TFile.Open("%s/%s/%s.root" % (inDir, systDir, fName), "read")
                    for obs in observables:
                        for reco in recoLvls:
                            recoObs = "%s_%s" % (reco,obs)
                            if not recoObs in templates[s][syst]:
                                templates[s][syst][recoObs] = {}
                            
                            _histname = "%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst)
                            if binning is not None:
                                tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__"+_histname)
                                templates[s][syst][recoObs][m] = tmp.Rebin(len(binning)-1, _histname, array('d', binning))
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
#                            if useVariableBinning and binning is None:
#                                templates[s][syst][recoObs][m].Scale(1., "width")
                            
#                            templates[s][syst][recoObs][m] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst))
#                            templates[s][syst][recoObs][m].SetDirectory(0)
                            
                            
                            # Rate scaling
                            templates[s][syst][recoObs][m].Scale(rateScaling)
        
                    f.Close()

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
    if binning is None:
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
                
                #print "Scaling %s %s" % (s, syst)
                for obs in observables:
                    for reco in recoLvls:
                        recoObs = "%s_%s" % (reco,obs)
                        for m in masses[s]:
                            templates[s][syst][recoObs][m].Scale(templates[s]["nominal"][recoObs][m].Integral() / templates[s][syst][recoObs][m].Integral())

    
   

    global tt_BinG, tt_BinMorphG,tttW_morphed
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

    tttW_morphed = {}
    tttW_G2D = {}


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
                    try:
                        tt_morphed[syst][recoObs],tt_BinG[syst][recoObs],tt_BinMorphG[syst][recoObs] = morphTemplates(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, variableBins=variableBins, interp=interp, morphRates=morphRates, verbosity=verbosity)
                    except KeyError:
                        pass

                try:
                    tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs] = morphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), precision=precision, systematic = "" if syst == "nominal" else "_"+syst, variableBins=variableBins, interp="pol1", morphRates=morphRates, verbosity=verbosity)
                except KeyError:
                    pass

                try:
                    tt_G2D[syst][recoObs],tt_GFit[syst][recoObs],tt_GErrors[syst][recoObs] = morphTemplates2D(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, variableBins=variableBins, interp=interp, morphRates=morphRates, verbosity=verbosity)
                except KeyError:
                    pass
                
                try:
                    tW_G2D[syst][recoObs],tW_GFit[syst][recoObs],tW_GErrors[syst][recoObs] = morphTemplates2D(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, variableBins=variableBins, interp="pol1", morphRates=morphRates, verbosity=verbosity)
                except KeyError:
                    pass


    # Create sum tt+tW templates
    if "TTbar" in signal and "ST_tW" in signal:
        templates["tt+tW"] = {}
        for syst,systDir in systematics.iteritems():
            
            if syst.find("Up") >= 0:
                systType = syst[:syst.find("Up")]
            elif syst.find("Down") >= 0:
                systType = syst[:syst.find("Down")]
            else:
                systType = syst
           
            templates["tt+tW"][syst] = {}
            tttW_morphed[syst] = {}

            for obs in observables:
                for reco in recoLvls:
                    recoObs = "%s_%s" % (reco,obs)
                    if not recoObs in templates["tt+tW"][syst]:
                        templates["tt+tW"][syst][recoObs] = {}
                    if not recoObs in tttW_morphed[syst]:
                        tttW_morphed[syst][recoObs] = {}

                    # Actual tt+tW templates
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
    data_obs = {}
    for obs in observables:
        for reco in recoLvls:
            recoObs = "%s_%s" % (reco,obs)
            #data_obs[recoObs] = templates["TTbar"]["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
            data_obs[recoObs] = tt_morphed["nominal"][recoObs][int(decimalScaling*172.5)].Clone("data_obs")
            data_obs[recoObs].SetTitle("data_obs")
            #data_obs[recoObs].Add(templates["ST_tW"]["nominal"][recoObs][int(decimalScaling*172.5)])
#            data_obs[recoObs].Add(tW_morphed["nominal"][recoObs][int(decimalScaling*172.5)])
            for b in background:
                # Omit WJets
                if b != "WJets":
                    data_obs[recoObs].Add(templates[b][recoObs])

    outFile = TFile.Open(outF, "recreate")
    for obs in observables:
        for reco in recoLvls:
            recoObs = "%s_%s" % (reco,obs)
            outFile.mkdir(recoObs)
            outFile.cd(recoObs)
            #outFile.mkdir("ttbins")
            #outFile.mkdir("tWbins")
            data_obs[recoObs].Write()
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

                for s in signal:
                    if s == "TTbar" and systType in tWOnlySysts: continue
                    if s == "ST_tW" and systType in ttOnlySysts: continue
                    for m in masses[s]:
                        templates[s][syst][recoObs][m].Write()
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
                            tt_BinMorphG[syst][recoObs][b].SetFillColor(0)
                            tt_BinG[syst][recoObs][b].SetFillColor(0)
                            tt_BinG[syst][recoObs][b].SetLineWidth(2)
                            tt_BinMorphG[syst][recoObs][b].Draw("ALP")
                            tt_BinG[syst][recoObs][b].SetLineColor(kBlue)
                            tt_BinG[syst][recoObs][b].Draw("LP SAME")
                            fit = tt_BinG[syst][recoObs][b].Fit(interp, "S" if verbosity > 1 else "SQ")
                            l.AddEntry(tt_BinMorphG[syst][recoObs][b], "Morphed")
                            l.AddEntry(tt_BinG[syst][recoObs][b], "Actual")
                            l.SetFillStyle(0)
                            l.Draw("SAME")
                            os.system("mkdir -p %s/%s/%s" % (plotDir,recoObs,syst))
                            c.SaveAs("%s/%s/%s/tt_morph_bin_%d.png" % (plotDir,recoObs,syst,b))
                    except KeyError:
                        pass

                    try:
                        for b in xrange(1, tW_morphed[syst][recoObs][morph_masses[0]].GetNbinsX()+1):
                            l = TLegend(0.425, 0.8, 0.575, 0.9)
                            tW_BinMorphG[syst][recoObs][b].SetFillColor(0)
                            tW_BinG[syst][recoObs][b].SetFillColor(0)
                            tW_BinG[syst][recoObs][b].SetLineWidth(2)
                            tW_BinMorphG[syst][recoObs][b].Draw("ALP")
                            tW_BinG[syst][recoObs][b].SetLineColor(kBlue)
                            tW_BinG[syst][recoObs][b].Draw("LP SAME")
                            fit = tW_BinG[syst][recoObs][b].Fit("pol1", "S" if verbosity > 1 else "SQ")
                            l.AddEntry(tW_BinMorphG[syst][recoObs][b], "Morphed")
                            l.AddEntry(tW_BinG[syst][recoObs][b], "Actual")
                            l.SetFillStyle(0)
                            l.Draw("SAME")
                            os.system("mkdir -p %s/%s/%s" % (plotDir,recoObs,syst))
                            c.SaveAs("%s/%s/%s/tW_morph_bin_%d.png" % (plotDir,recoObs,syst,b))
                    except KeyError:
                        pass

        print "Plots saved to %s" % plotDir

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
    return templates,tt_morphed,tW_morphed,tt_G2D,tt_GFit,tt_GErrors,diffSyst

#    for sys,sys_weight in systematics.items():
#        # Signals
#        tt[sys] = makeHist("%s/mc_TT_mt" % inDir, name = "ttactual", title="t#bar{t}", config=config, masses=ttmasses, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        tW[sys] = makeHist("%s/mc_ST_tW_top_mt" % inDir, name = "tWactual", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)
#        tWantitop[sys] = makeHist("%s/mc_ST_tW_antitop_mt" % inDir, name = "tW", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)
#
#        # Combine tW top & antitop
#        for obs,masses in tW[sys].items():
#            for m in masses.keys():
#                tW[sys][obs][m].Add(tWantitop[sys][obs][m])
#
#
#        # Backgrounds
#        DY[sys] = makeHist("%s/mc_DYJetsToLL" % inDir, name="DY", title="DY", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        WW[sys] = makeHist("%s/mc_WWTo2L2Nu" % inDir, name="WW", title="WW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        WZ[sys] = makeHist("%s/mc_WZTo3LNu" % inDir, name="WZ", title="WZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        ZZ[sys] = makeHist("%s/mc_ZZTo2L2Nu" % inDir, name="ZZ", title="ZZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        TTW[sys] = makeHist("%s/mc_TTWJetsToLNu" % inDir, name="TTW", title="TTW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        TTZ[sys] = makeHist("%s/mc_TTZToLLNuNu" % inDir, name="TTZ", title="TTZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        WJets[sys] = makeHist("%s/mc_WJetsToLNu" % inDir, name="WJets", title="WJets", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        STs[sys] = makeHist("%s/mc_ST_s" % inDir, name="STs", title="ST s", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        STt[sys] = makeHist("%s/mc_ST_t_top" % inDir, name="STt", title="ST t", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#        STtbar[sys] = makeHist("%s/mc_ST_t_antitop" % inDir, name="tbar", title="ST t antitop", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
#
#        # Combine t top & antitop
#        for obs in STt[sys].keys():
#            STt[sys][obs].Add(STtbar[sys][obs])
#
#
#        tt_BinG[sys] = {}
#        tW_BinG[sys] = {}
#        ttmorphed[sys] = {}
#        tWmorphed[sys] = {}
#        ttG2D[sys] = {}
#        tWG2D[sys] = {}
#        ttGFit[sys] = {}
#        tWGFit[sys] = {}
#        ttGerrors[sys] = {}
#        tWGerrors[sys] = {}
#
#        for obs in config.keys():
#            # Returns morphed templates, bin graphs
#            ttmorphed[sys][obs],tt_BinG[sys][obs] = morphTemplates(tt[sys][obs], morph_masses, name=obs+"_tt", title = config[obs]["title"] + " t#bar{t}", systematic = "" if sys == "nominal" else "_"+sys, interp=interp)
#            tWmorphed[sys][obs],tW_BinG[sys][obs] = morphTemplates(tW[sys][obs], morph_masses, name=obs+"_tW", title = config[obs]["title"] + " tW", systematic = "" if sys == "nominal" else "_"+sys, interp="pol1")
#            
#            ttG2D[sys][obs],ttGFit[sys][obs],ttGerrors[sys][obs] = morphTemplates2D(tt[sys][obs], morph_masses, name=obs+"_tt", title = config[obs]["title"] + " t#bar{t}", systematic = "" if sys == "nominal" else "_"+sys, interp=interp)
#            tWG2D[sys][obs],tWGFit[sys][obs], tWGerrors[sys][obs] = morphTemplates2D(tW[sys][obs], morph_masses, name=obs+"_tW", title = config[obs]["title"] + " tW", systematic = "" if sys == "nominal" else "_"+sys, interp="pol1")
#
#    #pprint(ttmorphed)
#    #pprint(tWmorphed)
#   
#    #pprint(tt_BinG)
#    #pprint(tW_BinG)
#
#    
#    return ttG2D,ttGFit,ttGerrors


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--topDir", default="/uscms_data/d3/msaunder/TopMass/CMSSW_8_0_26_patch1/src/TopNtuplizer/Plotting", help="path to Plotting directory")
    parser.add_argument("-i", dest="inDir", default="histograms", help="Input directory containing ttrees")
    parser.add_argument("-o", dest="outF", default="mtTemplatesForCH.root", help="Output template file")
    parser.add_argument("--debugOut", default="", help="output file to store pickled template info")
    parser.add_argument("-b", "--rebin", type=int, default=1, help="Integer rebin width in GeV")
    parser.add_argument("--systs", default="", nargs="*", choices=(allSystematics + ["none","None"]), help="ONLY plot these systematics")
    parser.add_argument("--cutMin", type=int, default=0, help="bins to cut from the left in GeV (before rebinning)")
    parser.add_argument("--cutMax", type=int, default=0, help="bins to cut from the right in GeV (before rebinning)")
    parser.add_argument("--minmt", type=float, default=166.5, help="minimum mass for morphing range")
    parser.add_argument("--maxmt", type=float, default=178.5, help="maximum mass for morphing range")
    parser.add_argument("--deltaMT", type=float, default=0.1, help="morphing mass increment (in GeV)") 
    #parser.add_argument("--obs", nargs="+", default=[], help="include these observables")
    parser.add_argument("--morphRates", action="store_true", default=False, help="interpolate non-normalized histograms")
    parser.add_argument("--noScalingToNominal", action="store_true", default=False, help="don't scale certain systematics to the nominal rate")
    parser.add_argument("--useMorphFile", action="store_true", default=False, help="use morphed templates from external file instead of doing morphing here")
    parser.add_argument("--extMorphFile", default="/uscms/homes/m/msaunder/work/DNN_TemplateMorphing/new/newplots/morphTF.root", help="external morph template file loaded when useMorphFile option is selected")
    parser.add_argument("-r", "--rateScaling", type=float, default=1., help="scale rates by a constant factor")
#    parser.add_argument("--bins", type=str, default="[0, 15, 30, 45, 60, 70, 80, 90, 100, 110, 120, 130, 140]", help="List of variable bin ranges. The last entry is the upper edge of the last bin. All other entries are the lower bin edges")
    parser.add_argument("--bins", type=str, default="", help="List of variable bin ranges. The last entry is the upper edge of the last bin. All other entries are the lower bin edges")
    parser.add_argument("--newErrors", action="store_true", default=False, help="use new error method")
    parser.add_argument("--precision", type=int, default=1, help="number of decimal places to use for morphing")
    parser.add_argument("--obs", nargs="+", default=['ptll'], choices=(_observables + ["all"]), help="include these observables (or 'all' for all 6)")
    parser.add_argument("--reco", nargs="+", default=["rec"], choices=["rec","gen"], help="reco lvl")
    parser.add_argument("--interp", default="pol3", help="ttbar interpolation function to use in ROOT Fit")
    parser.add_argument("--plots", action = "store_true", help="create bin plots")
    parser.add_argument("--plotDir", default="morphed_bins", help="directory to store bin plots if --plots is selected")
    parser.add_argument("--includeGraphs", action="store_true", default=False, help="Store per-bin graphs in output root file")
    parser.add_argument("-v", "-V", "--verbosity", dest="verbosity", type=int, default=0, help="verbosity of output")
    args = parser.parse_args()

    scaleToNominal = not args.noScalingToNominal
    if "all" in args.obs:
        # Include all observables
        args.obs = _observables

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
        if args.verbosity > 1: print "Rebinning to %d GeV" % args.rebin

    
    # Create a set of templates for the given ttres and config 
    #ttG2D,ttGFit,ttGerrors = create_templates(inDir=args.inDir, rebin=args.rebin, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir)
    templates,tt_morphed,tW_morphed,tt_G2D,tt_GFit,tt_GErrors, diffSyst = create_templates(inDir="%s/%s" % (args.topDir,args.inDir), includedSysts=args.systs, cutMin=args.cutMin, cutMax=args.cutMax, massMin=args.minmt, massMax=args.maxmt, deltaMT=args.deltaMT, rateScaling=args.rateScaling, observables=args.obs, recoLvls=args.reco, rebin=args.rebin, bins=args.bins, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir, debugOut=args.debugOut, includeGraphs=args.includeGraphs, morphRates=args.morphRates, useMorphFile=args.useMorphFile, extMorphFile=args.extMorphFile, scaleToNominal=scaleToNominal,verbosity=args.verbosity)



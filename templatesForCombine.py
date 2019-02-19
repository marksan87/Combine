#!/usr/bin/env python
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

gROOT.SetBatch(True)

_observables = ['ptll', 'Mll', 'ptpos', 'Epos', 'ptp_ptm', 'Ep_Em'] 
obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}
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
               }

separateSystSamples = ['isr','fsr','DS','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
ttOnlySysts = ['toppt','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

PU_mt1725only = False 


def clamp(val, minV, maxV):
    # Clamps a value between min and max 
    return max(minV, min(val,maxV))


def systMorph(histMT, nominalMeanBin, diffHist):
    mtMeanBin = histMT.FindBin(histMT.GetMean())
    binShift = 0
    if mtMeanBin != nominalMeanBin:
        binShift = mtMeanBin - nominalMeanBin
#        print "Mean shifted by %s" % binShift

    Nbins = diffHist.GetNbinsX()
    morphed = histMT.Clone()
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
                hist[obs][m] = TH1F("%s%s%d%s" % (obs,name,m,systematic), "%s  %s  m_{t} = %.1f GeV" % (vals["title"], title, m/10.), vals["nbins"], vals["min"], vals["max"])
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


def morphTemplates(templates, morph_masses, name, title, systematic="", interp = "pol3", verbose=False):
    morph = {}
    fit = {}
    fitFunc = {}
    binG = {}
    binMorphG = {}
    rates = {}
    normalized = {}
    actual_masses = sorted(templates.keys())
    #print "morphTemplates"
#    print "name =", name
    print "Morphing  %s  %s" % (name,systematic[1:])
    obs = name[4:-3]
    #print "obs =", obs
    for m in actual_masses:
        #print "Now on %s  %s  mt = %.1f" % (name,systematic,m)
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone()
        normalized[m].SetDirectory(0)
        normalized[m].Scale(1.0 / rates[m])


    for b in xrange(1, normalized[actual_masses[0]].GetNbinsX()+1):
        binG[b] = TGraphErrors(len(actual_masses), array('d', [m/10.0 for m in actual_masses]), array('d', [normalized[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [normalized[m].GetBinError(b) for m in actual_masses]))
        binG[b].SetName("%s%sbin_%d" % (name, "" if systematic == "" else systematic+"_",b))
        binG[b].SetTitle("%s  Bin %d" % (obs,b))
        binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binG[b].GetYaxis().SetTitle("Entries")
        binG[b].GetYaxis().SetTitleOffset(1.3)
        binG[b].Draw()
        fit[b] = binG[b].Fit(interp, "S"+ ("" if verbose else "Q"))
        fitFunc[b] = binG[b].GetFunction(interp)
        errors = array('d', [0.] * len(morph_masses))
        fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/10.0 for m in morph_masses]), errors, 2./3., False)
        for i,m in enumerate(morph_masses):
            if m not in morph:
                morph[m] = TH1F("%s%d%s" % (name,m,systematic), "%s  m_{t} = %.1f GeV" % (title, m/10.), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
            morph[m].SetBinContent(b, fitFunc[b].Eval(m/10.0))
            morph[m].SetBinError(b, errors[i])
    
    for b in xrange(1, morph[morph_masses[0]].GetNbinsX()+1):
        binMorphG[b] = TGraphErrors(len(morph_masses), array('d', [m/10.0 for m in morph_masses]), array('d', [morph[m].GetBinContent(b) for m in morph_masses]), array('d', [0.] * len(morph_masses)), array('d', [morph[m].GetBinError(b) for m in morph_masses]))
        binMorphG[b].SetName("%s%smorphed_bin_%d" % (name,"" if systematic == "" else systematic+"_",b))
        binMorphG[b].SetTitle("Morphed Bin %d" % b)
        binMorphG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binMorphG[b].GetYaxis().SetTitle("Entries")
        binMorphG[b].GetYaxis().SetTitleOffset(1.3)

    
    #print "Done with %s  %s" % (name,systematic[1:])    
    for m in morph_masses:
        morph[m].Scale(rates[1725])


    return morph,binG,binMorphG



def morphTemplates2D(templates, morph_masses, name, title, systematic="", interp = "pol3", verbose=False):
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
        normalized[m].Scale(1.0 / rates[m])

    #print "morphTemplates2D"
    obs = name[4:-3]
    nbins_obs = templates[actual_masses[0]].GetNbinsX()
    nbins_mass = len(morph_masses)

    graph2D = TGraph2DErrors(nbins_obs*nbins_mass)
    for m in actual_masses:
        for b in xrange(1, nbins_obs):
            N = graph2D.GetN()+1
            graph2D.SetPoint(N, m / 10.0, normalized[m].GetBinCenter(b), normalized[m].GetBinContent(b))
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
    fit = graph2D.Fit(interp, "S" + ("" if verbose else "Q"))
    fitFunc = graph2D.FindObject(interp)
    #fitFunc = graph2D.GetFunction(interp)
    errors = array('d', [0.] * len(morph_masses))
    fit.GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/10.0 for m in morph_masses]), errors, 2./3., False)
    return graph2D, fit, errors



def create_templates(inDir, rebin, cutMin, cutMax, obsList, recoLvls, interp, outF, makePlots, plotDir, debugOut, includeGraphs, verbose=False): 
    masses = {}
    masses["TTbar"] = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
    masses["ST_tW"] = [1695, 1725, 1755]
    deltaM = 1
    morph_masses = range(masses["TTbar"][0], masses["TTbar"][-1] + deltaM, deltaM)
    print "Processing templates from %s" % inDir 
     

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

    if len(obsList) == 0:
        # Use all observables
        observables = _observables
    else:
        # Use only these observables
        observables = obsList

   
    
    diffSyst = {}   # Ratio of syst/nominal @ mt=172.5 for separate sample systematics


    # Backgrounds
    for b in background:
        templates[b] = {}
        f = TFile.Open("%s/%s/%s.root" % (inDir, systematics["nominal"], b), "read")
        for obs in observables:
            for reco in recoLvls:
                recoObs = "%s_%s" % (reco,obs)  # rec_ptll, etc..
                #templates[b][recoObs] = f.Get("%s_%s" % (recoObs, b)).Clone(b)
                #templates[b][recoObs].SetDirectory(0)
                if cutMin == 0 and cutMax == 0:
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

            if systType in separateSystSamples:
                # Samples only exist for mt = 172.5
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
                        if cutMin == 0 and cutMax == 0:
                            templates[s][syst][recoObs][1725] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,1725,"" if syst=="nominal" else "_%s" % syst))
                            templates[s][syst][recoObs][1725].SetDirectory(0)
                        else:
                            tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__%s%d%s" % (sample,1725,"" if syst=="nominal" else "_%s" % syst)) 
                            binW = tmp.GetBinCenter(2) - tmp.GetBinCenter(1)
                            templates[s][syst][recoObs][1725] = TH1F(b, tmp.GetTitle(), tmp.GetNbinsX()-cutMin-cutMax, tmp.GetBinCenter(1)+cutMin-binW/2., tmp.GetBinCenter(tmp.GetNbinsX())-cutMax+binW/2.)
                            templates[s][syst][recoObs][1725].SetDirectory(0)
                            for _bin in xrange(templates[s][syst][recoObs][1725].GetNbinsX()+1):
                                templates[s][syst][recoObs][1725].SetBinContent(_bin, tmp.GetBinContent(_bin+cutMin))
                                templates[s][syst][recoObs][1725].SetBinError(_bin, tmp.GetBinError(_bin+cutMin))
#                        templates[s][syst][recoObs][1725] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,1725,"" if syst=="nominal" else "_%s" % syst))
#                        templates[s][syst][recoObs][1725].SetDirectory(0)

                        templates[s][syst][recoObs][1725].SetTitle("%s %s %s  m_{t} = %.1f GeV" % (reco, obsTitle[obs], "t#bar{t}" if s == "TTbar" else "tW", 172.5)) 

            else:
                for m in masses[s]:
                    fName = s if m == 1725 else "%s_mt%d" % (s,m)
                    f = TFile.Open("%s/%s/%s.root" % (inDir, systDir, fName), "read")
                    for obs in observables:
                        for reco in recoLvls:
                            recoObs = "%s_%s" % (reco,obs)
                            if not recoObs in templates[s][syst]:
                                templates[s][syst][recoObs] = {}
                            
                            if cutMin == 0 and cutMax == 0:
                                templates[s][syst][recoObs][m] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst))
                                templates[s][syst][recoObs][m].SetDirectory(0)
                            else:
                                tmp = f.Get("%s_%s" % (recoObs,fName)).Clone("__%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst))
                                binW = tmp.GetBinCenter(2) - tmp.GetBinCenter(1)
                                templates[s][syst][recoObs][m] = TH1F(b, tmp.GetTitle(), tmp.GetNbinsX()-cutMin-cutMax, tmp.GetBinCenter(1)+cutMin-binW/2., tmp.GetBinCenter(tmp.GetNbinsX())-cutMax+binW/2.)
                                templates[s][syst][recoObs][m].SetDirectory(0)
                                for _bin in xrange(templates[s][syst][recoObs][m].GetNbinsX()+1):
                                    templates[s][syst][recoObs][m].SetBinContent(_bin, tmp.GetBinContent(_bin+cutMin))
                                    templates[s][syst][recoObs][m].SetBinError(_bin, tmp.GetBinError(_bin+cutMin))                            
                            
                            templates[s][syst][recoObs][m].SetTitle("%s %s %s  m_{t} = %.1f GeV" % (reco, obsTitle[obs], "t#bar{t}" if s == "TTbar" else "tW", m/10.)) 
                            
                            
#                            templates[s][syst][recoObs][m] = f.Get("%s_%s" % (recoObs,fName)).Clone("%s%d%s" % (sample,m,"" if syst=="nominal" else "_%s" % syst))
#                            templates[s][syst][recoObs][m].SetDirectory(0)

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

            if systType in separateSystSamples:
                print "Creating systematic morphed templates for %s %s" % (s,syst)
                diffSyst[s][syst] = {}
                for m in masses[s]:
                    if m == 1725: continue

                    for obs in observables:
                        for reco in recoLvls:
                            recoObs = "%s_%s" % (reco,obs)
                        
                            diffSyst[s][syst][recoObs] = templates[s][syst][recoObs][1725].Clone()
                            diffSyst[s][syst][recoObs].Divide(templates[s]['nominal'][recoObs][1725])

                            templates[s][syst][recoObs][m] = systMorph(templates[s]["nominal"][recoObs][m], templates[s]["nominal"][recoObs][1725].FindBin(templates[s]["nominal"][recoObs][1725].GetMean()), diffSyst[s][syst][recoObs])


    # Rebin prior to morphing
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

    global tt_BinG, tt_BinMorphG
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
                try:
                    tt_morphed[syst][recoObs],tt_BinG[syst][recoObs],tt_BinMorphG[syst][recoObs] = morphTemplates(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp=interp, verbose=verbose)
                except KeyError:
                    pass

                try:
                    tW_morphed[syst][recoObs],tW_BinG[syst][recoObs],tW_BinMorphG[syst][recoObs] = morphTemplates(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", verbose=verbose)
                except KeyError:
                    pass

                try:
                    tt_G2D[syst][recoObs],tt_GFit[syst][recoObs],tt_GErrors[syst][recoObs] = morphTemplates2D(templates["TTbar"][syst][recoObs], morph_masses, name=recoObs+"_tt", title = "%s %s t#bar{t}" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp=interp, verbose=verbose)
                except KeyError:
                    pass
                
                try:
                    tW_G2D[syst][recoObs],tW_GFit[syst][recoObs],tW_GErrors[syst][recoObs] = morphTemplates2D(templates["ST_tW"][syst][recoObs], morph_masses, name=recoObs+"_tW", title = "%s %s tW" % (reco, obsTitle[obs]), systematic = "" if syst == "nominal" else "_"+syst, interp="pol1", verbose=verbose)
                except KeyError:
                    pass


    data_obs = {}
    for obs in observables:
        for reco in recoLvls:
            recoObs = "%s_%s" % (reco,obs)
            #data_obs[recoObs] = templates["TTbar"]["nominal"][recoObs][1725].Clone("data_obs")
            data_obs[recoObs] = tt_morphed["nominal"][recoObs][1725].Clone("data_obs")
            data_obs[recoObs].SetTitle("data_obs")
            #data_obs[recoObs].Add(templates["ST_tW"]["nominal"][recoObs][1725])
            data_obs[recoObs].Add(tW_morphed["nominal"][recoObs][1725])
            for b in background:
                data_obs[recoObs].Add(templates[b][recoObs])

    print ""
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
                    if PU_mt1725only and "pileup" in syst and m != 1725: continue
                    
                    try:
                        tt_morphed[syst][recoObs][m].Write(tt_morphed[syst][recoObs][m].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass

                    try:
                        tW_morphed[syst][recoObs][m].Write(tW_morphed[syst][recoObs][m].GetName()[len(recoObs)+1:])
                    except KeyError:
                        pass
    
                for s in signal:
                    for m in masses[s]:
                        try:
                            templates[s][syst][recoObs][m].Write()
                        except KeyError:
                            pass
                        #templates[s][syst][recoObs][m].Write(templates[s][syst][recoObs][m].GetName()[len(recoObs):])
                
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
    print "Templates saved to %s" % outF

    if makePlots:
        print "Making per-bin plots"
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
                            fit = tt_BinG[syst][recoObs][b].Fit(interp, "S" if verbose else "SQ")
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
                            fit = tW_BinG[syst][recoObs][b].Fit("pol1", "S" if verbose else "SQ")
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
    print "Debugging info saved to %s" % debugOut 
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
    parser.add_argument("--debugOut", default="debug_templates.pklz", help="output file to store pickled template info")
    parser.add_argument("-b", "--rebin", type=int, default=2, help="Integer rebin width in GeV")
    parser.add_argument("--cutMin", type=int, default=0, help="bins to cut from the left in GeV (before rebinning)")
    parser.add_argument("--cutMax", type=int, default=0, help="bins to cut from the right in GeV (before rebinning)")
    parser.add_argument("--obs", nargs="+", default=[], help="include these observables")
    parser.add_argument("--reco", nargs="+", default=["rec"], choices=["rec","gen"], help="reco lvl")
    parser.add_argument("--interp", default="pol3", help="ttbar interpolation function to use in ROOT Fit")
    parser.add_argument("--plots", action = "store_true", help="create bin plots")
    parser.add_argument("--plotDir", default="bins", help="directory to store bin plots if --plots is selected")
    parser.add_argument("--includeGraphs", action="store_true", default=False, help="Store per-bin graphs in output root file")
    parser.add_argument("-v", "-V", "--verbose", dest="verbose", action="store_true", default=False, help="verbose output")
    args = parser.parse_args()

    if args.topDir[-1] == "/": args.topDir = args.topDir[:-1]
    if args.inDir[-1] == "/": args.inDir = args.inDir[:-1]
    if args.rebin < 1:
        print "Invalid rebin value. Must be >= 1! Defaulting to 2 GeV"
        args.rebin = 2
    else:
        print "Rebinning to %d GeV" % args.rebin

    
    # Create a set of templates for the given ttres and config 
    #ttG2D,ttGFit,ttGerrors = create_templates(inDir=args.inDir, rebin=args.rebin, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir)
    templates,tt_morphed,tW_morphed,tt_G2D,tt_GFit,tt_GErrors, diffSyst = create_templates(inDir="%s/%s" % (args.topDir,args.inDir), cutMin=args.cutMin, cutMax=args.cutMax, obsList=args.obs, recoLvls=args.reco, rebin=args.rebin, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir, debugOut=args.debugOut, includeGraphs=args.includeGraphs, verbose=args.verbose)



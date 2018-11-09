#!/usr/bin/env python
import ROOT
from ROOT import *
from argparse import ArgumentParser
from array import array
import os
from time import sleep
from pprint import pprint
import pickle

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

sysList = {"toppt", "toppt_Up"}

def makeHist(fName, name, title, config, masses = [], systematic = "", weight="weight"):
    hist = {}
    if len(masses) > 0:
        for m in masses:
            f = TFile.Open("%s%d.root" % (fName,m), "read")
            tree = f.Get("goodEvents")
            for obs,vals in config.items():
                if obs not in hist: hist[obs] = {}
                hist[obs][m] = TH1D("%s%s%d%s" % (obs,name,m,systematic), "%s  %s  m_{t} = %.1f GeV" % (vals["title"], title, m/10.), vals["nbins"], vals["min"], vals["max"])
                tree.Draw("%s>>%s" % (obs,hist[obs][m].GetName()), "%s * (%s > %f && %s < %f)" % (weight, obs, vals["min"], obs, vals["max"]))
                hist[obs][m].SetDirectory(0)
            
            f.Close() 
    else:
        f = TFile.Open("%s.root" % fName, "read")
        tree = f.Get("goodEvents")
        for obs,vals in config.items():
            hist[obs] = TH1D("%s%s%s" % (obs,name,systematic), "%s  %s" % (vals["title"], title), vals["nbins"], vals["min"], vals["max"])
            tree.Draw("%s>>%s" % (obs,hist[obs].GetName()), "%s * (%s > %f && %s < %f)" % (weight, obs, vals["min"], obs, vals["max"]))
            hist[obs].SetDirectory(0)
        
        f.Close() 
    return hist


def morphTemplates(templates, morph_masses, name, title, systematic="", interp = "pol3"):
    morph = {}
    fit = {}
    fitFunc = {}
    binG = {}
    rates = {}
    normalized = {}
    actual_masses = sorted(templates.keys())
    print "morphTemplates"
    print "name =", name
    obs = name[4:-3]
    print "obs =", obs
    for m in actual_masses:
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

        fit[b] = binG[b].Fit(interp, "S")
        fitFunc[b] = binG[b].GetFunction(interp)
        errors = array('d', [0.] * len(morph_masses))
        fit[b].GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/10.0 for m in morph_masses]), errors, 2./3., False)
        for i,m in enumerate(morph_masses):
            if m not in morph:
                morph[m] = TH1D("%s%d%s" % (name,m,systematic), "%s  m_{t} = %.1f GeV" % (title, m/10.), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
            morph[m].SetBinContent(b, fitFunc[b].Eval(m/10.0))
            morph[m].SetBinError(b, errors[i])
    
    for m in morph_masses:
        morph[m].Scale(rates[1725])

    return morph,binG



def morphTemplates2D(templates, morph_masses, name, title, systematic="", interp = "pol3"):
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

    print "morphTemplates2D"
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
    #print "About to fit 2D graph %s_%s" % (name,systematic)
    print "About to fit 2D graph %s" % graph2D.GetName()
    fit = graph2D.Fit(interp, "S")
    fitFunc = graph2D.FindObject(interp)
    #fitFunc = graph2D.GetFunction(interp)
    errors = array('d', [0.] * len(morph_masses))
    fit.GetConfidenceIntervals(len(morph_masses), 1, 1, array('d', [m/10.0 for m in morph_masses]), errors, 2./3., False)
    return graph2D, fit, errors


    for i,m in enumerate(morph_masses):
        if m not in morph:
            morph[m] = TH1D("%s%d%s" % (name,m,systematic), "%s  m_{t} = %.1f GeV" % (title, m/10.), normalized[actual_masses[0]].GetNbinsX(), normalized[actual_masses[0]].GetXaxis().GetXmin(), normalized[actual_masses[0]].GetXaxis().GetXmax())
        morph[m].SetBinContent(b, fitFunc[b].Eval(m/10.0))
        morph[m].SetBinError(b, errors[i])
    
    for m in morph_masses:
        morph[m].Scale(rates[1725])

    return morph,binG


def create_templates(inDir, configF, interp, outF, makePlots, plotDir): 
    # Load config file for distributions and bin min,max,count
    config = {}
    with open(configF, "r") as f:
        for i,line in enumerate(f):
            if i == 0 and line[0] == "#": continue
            l = line.split()
            config[l[0]] = {"nbins":int(l[1]), "min":float(l[2]), "max":float(l[3]), "title":obsTitle[l[0][4:]]}

    print "\nFound observables from %s:\n" % args.config
    for obs,vals in config.items():
        print "%s\t%d bins from %.1f-%.1f GeV  (bin width = %.1f GeV)" % (obs if len(obs) >= 8 else obs+"\t", vals["nbins"], vals["min"], vals["max"], (vals["max"] - vals["min"]) / float(vals["nbins"])) 
    print "\n"

    ttmasses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
    tWmasses = [1695, 1725, 1755]
    deltaM = 1
    morph_masses = range(ttmasses[0], ttmasses[-1] + deltaM, deltaM)
    print "Processing ttrees from %s" % inDir

#gROOT.SetBatch(True)

# Systematics
    systematics = {"nominal":"weight", "topptDown":"weight", "topptUp":"weight_TopPtReweight", "pileupUp":"weight_pileupUp", "pileupDown":"weight_pileupDown"}
    tt = {}
    tW = {}
    tWantitop = {}
    DY = {}
    WW = {}
    WZ = {}
    ZZ = {}
    TTW = {}
    TTZ = {}
    WJets = {}
    STs = {}
    STt = {}
    STtbar = {}


    tt_BinG = {}
    tW_BinG = {}
    ttmorphed = {}
    tWmorphed = {}

    ttG2D = {}
    tWG2D = {}
    ttGFit = {}
    tWGFit = {}
    ttGerrors = {}
    tWGerrors = {}

    for sys,sys_weight in systematics.items():
        # Signals
        tt[sys] = makeHist("%s/mc_TT_mt" % inDir, name = "ttactual", title="t#bar{t}", config=config, masses=ttmasses, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        tW[sys] = makeHist("%s/mc_ST_tW_top_mt" % inDir, name = "tWactual", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)
        tWantitop[sys] = makeHist("%s/mc_ST_tW_antitop_mt" % inDir, name = "tW", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)

        # Combine tW top & antitop
        for obs,masses in tW[sys].items():
            for m in masses.keys():
                tW[sys][obs][m].Add(tWantitop[sys][obs][m])


        # Backgrounds
        DY[sys] = makeHist("%s/mc_DYJetsToLL" % inDir, name="DY", title="DY", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        WW[sys] = makeHist("%s/mc_WWTo2L2Nu" % inDir, name="WW", title="WW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        WZ[sys] = makeHist("%s/mc_WZTo3LNu" % inDir, name="WZ", title="WZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        ZZ[sys] = makeHist("%s/mc_ZZTo2L2Nu" % inDir, name="ZZ", title="ZZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        TTW[sys] = makeHist("%s/mc_TTWJetsToLNu" % inDir, name="TTW", title="TTW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        TTZ[sys] = makeHist("%s/mc_TTZToLLNuNu" % inDir, name="TTZ", title="TTZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        WJets[sys] = makeHist("%s/mc_WJetsToLNu" % inDir, name="WJets", title="WJets", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        STs[sys] = makeHist("%s/mc_ST_s" % inDir, name="STs", title="ST s", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        STt[sys] = makeHist("%s/mc_ST_t_top" % inDir, name="STt", title="ST t", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
        STtbar[sys] = makeHist("%s/mc_ST_t_antitop" % inDir, name="tbar", title="ST t antitop", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)

        # Combine t top & antitop
        for obs in STt[sys].keys():
            STt[sys][obs].Add(STtbar[sys][obs])


        tt_BinG[sys] = {}
        tW_BinG[sys] = {}
        ttmorphed[sys] = {}
        tWmorphed[sys] = {}
        ttG2D[sys] = {}
        tWG2D[sys] = {}
        ttGFit[sys] = {}
        tWGFit[sys] = {}
        ttGerrors[sys] = {}
        tWGerrors[sys] = {}

        for obs in config.keys():
            # Returns morphed templates, bin graphs
            ttmorphed[sys][obs],tt_BinG[sys][obs] = morphTemplates(tt[sys][obs], morph_masses, name=obs+"_tt", title = config[obs]["title"] + " t#bar{t}", systematic = "" if sys == "nominal" else "_"+sys, interp=interp)
            tWmorphed[sys][obs],tW_BinG[sys][obs] = morphTemplates(tW[sys][obs], morph_masses, name=obs+"_tW", title = config[obs]["title"] + " tW", systematic = "" if sys == "nominal" else "_"+sys, interp="pol1")
            
            ttG2D[sys][obs],ttGFit[sys][obs],ttGerrors[sys][obs] = morphTemplates2D(tt[sys][obs], morph_masses, name=obs+"_tt", title = config[obs]["title"] + " t#bar{t}", systematic = "" if sys == "nominal" else "_"+sys, interp=interp)
            tWG2D[sys][obs],tWGFit[sys][obs], tWGerrors[sys][obs] = morphTemplates2D(tW[sys][obs], morph_masses, name=obs+"_tW", title = config[obs]["title"] + " tW", systematic = "" if sys == "nominal" else "_"+sys, interp="pol1")

    #pprint(ttmorphed)
    #pprint(tWmorphed)
   
    #pprint(tt_BinG)
    #pprint(tW_BinG)

    data_obs = {}
    for obs in config.keys():
        data_obs[obs] = tt["nominal"][obs][1725].Clone("data_obs")
        data_obs[obs].SetTitle("data_obs")
        data_obs[obs].Add(tW["nominal"][obs][1725])
        data_obs[obs].Add(DY["nominal"][obs])
        data_obs[obs].Add(WW["nominal"][obs])
        data_obs[obs].Add(WZ["nominal"][obs])
        data_obs[obs].Add(ZZ["nominal"][obs])
        data_obs[obs].Add(WJets["nominal"][obs])
        data_obs[obs].Add(STs["nominal"][obs])
        data_obs[obs].Add(STt["nominal"][obs])
        data_obs[obs].Add(TTZ["nominal"][obs])
        data_obs[obs].Add(TTW["nominal"][obs])


    outFile = TFile.Open(outF, "recreate")
    for obs in config.keys(): 
        #outFile.cd()
        #gDirectory.cd()
        outFile.mkdir(obs)
        outFile.cd(obs)
        #outFile.mkdir("ttbins")
        #outFile.mkdir("tWbins")
        data_obs[obs].Write()
        for sys in systematics.keys():
            #gDirectory.cd("%s/ttbins" % obs)
            for b in xrange(1,len(tt_BinG[sys][obs])):
                tt_BinG[sys][obs][b].Write(tt_BinG[sys][obs][b].GetName()[len(obs)+1:])
            #gDirectory.cd("../tWbins")
            for b in xrange(1,len(tW_BinG[sys][obs])):
                tW_BinG[sys][obs][b].Write(tW_BinG[sys][obs][b].GetName()[len(obs)+1:])
            #gDirectory("../")
            for m in morph_masses:
                ttmorphed[sys][obs][m].Write(ttmorphed[sys][obs][m].GetName()[len(obs)+1:])
                tWmorphed[sys][obs][m].Write(tWmorphed[sys][obs][m].GetName()[len(obs)+1:])
            for m in ttmasses:
                tt[sys][obs][m].Write(tt[sys][obs][m].GetName()[len(obs):])
            for m in tWmasses:
                tW[sys][obs][m].Write(tW[sys][obs][m].GetName()[len(obs):])
            #outFile.cd()
            #outFile.cd(obs)
            #gDirectory.cd(obs)
            ttG2D[sys][obs].Write()
            tWG2D[sys][obs].Write()
            DY[sys][obs].Write(DY[sys][obs].GetName()[len(obs):])
            WW[sys][obs].Write(WW[sys][obs].GetName()[len(obs):])
            WZ[sys][obs].Write(WZ[sys][obs].GetName()[len(obs):])
            ZZ[sys][obs].Write(ZZ[sys][obs].GetName()[len(obs):])
            TTW[sys][obs].Write(TTW[sys][obs].GetName()[len(obs):])
            TTZ[sys][obs].Write(TTZ[sys][obs].GetName()[len(obs):])
            WJets[sys][obs].Write(WJets[sys][obs].GetName()[len(obs):])
            STs[sys][obs].Write(STs[sys][obs].GetName()[len(obs):])
            STt[sys][obs].Write(STt[sys][obs].GetName()[len(obs):])
    outFile.Close()

    print "Templates saved to %s" % outF
    with open("debug_templates.pkl", "wb") as f:
        pickle.dump(ttG2D, f)
        pickle.dump(ttGFit, f)
        pickle.dump(ttGerrors, f)
        pickle.dump(tWG2D, f)
        pickle.dump(tWGFit, f)
        pickle.dump(tWGerrors, f)
    print "Debugging info saved to debug_templates.pkl" 
    
    return ttG2D,ttGFit,ttGerrors
    """ 
    if makePlots:
        # Save plots of individual bins and differences for actual mass points
        gROOT.SetBatch(True)
        os.system("mkdir -p %s" % plotDir)
        c = TCanvas("foo","bar", 2000, 1200)
        line = TF1("line", "0", -99999, 99999)
        gStyle.SetOptStat(0)
        #for m in masses:
        #    diff[m].GetYaxis().SetRangeUser(-0.3, 0.3)
        #    diff[m].Draw()
        #    c.Update()
        #    tline = TLine(c.GetUxmin(), 0, c.GetUxmax(), 0)
        #    tline.Draw("SAME")
        #    c.SaveAs("%s/diff_mt%d.png" % (plotDir, m))

        for obs,vals in config.items():
            # Create seperate bin plots for each observable
            os.system("mkdir -p %s/%s" % (plotDir, obs)
            for b in xrange(1, vals["nbins"]+1):
                l = TLegend(0.425, 0.8, 0.575, 0.9)
                binMorphG[b].SetFillColor(0)
                binG[b].SetFillColor(0)
                binG[b].SetLineWidth(2)
                binMorphG[b].Draw("ALP")
                binG[b].SetLineColor(kBlue)
                binG[b].Draw("LP SAME")
                fit = binG[b].Fit(args.func, "S")
                l.AddEntry(binMorphG[b], "Morphed")
                l.AddEntry(binG[b], "Actual")
                l.SetFillStyle(0)
                #l.AddEntry(fit.Get(), "Fit")
                l.Draw("SAME")
                c.SaveAs("%s/%s/morph_bin_%d.png" % (plotDir,b))
    """


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", dest="inDir", default="plots2018", help="Input directory containing ttrees")
    parser.add_argument("-o", dest="outF", default="mtTemplatesForCH.root", help="Output template file")
    parser.add_argument("-c", dest="config", default="config.txt", help="Config file specifying observables and bin ranges/sizes")
    parser.add_argument("--interp", default="pol3", help="ttbar interpolation function to use in ROOT Fit")
    parser.add_argument("--plots", action = "store_true", help="create bin plots")
    parser.add_argument("--plotDir", default="bins", help="directory to store bin plots if --plots is selected")
    args = parser.parse_args()

    if args.inDir[-1] == "/": args.inDir = args.inDir[:-1]
    
    # Create a set of templates for the given ttres and config 
    ttG2D,ttGFit,ttGerrors = create_templates(inDir=args.inDir, configF=args.config, interp=args.interp, outF=args.outF, makePlots=args.plots, plotDir=args.plotDir)



#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
from array import array
import os
from time import sleep

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}

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
    for m in actual_masses:
        rates[m] = templates[m].Integral()
        normalized[m] = templates[m].Clone()
        normalized[m].SetDirectory(0)
        normalized[m].Scale(1.0 / rates[m])


    for b in xrange(1, normalized[actual_masses[0]].GetNbinsX()+1):
        binG[b] = TGraphErrors(len(actual_masses), array('d', [m/10.0 for m in actual_masses]), array('d', [normalized[m].GetBinContent(b) for m in actual_masses]), array('d', [0.] * len(actual_masses)), array('d', [normalized[m].GetBinError(b) for m in actual_masses]))
        binG[b].SetName("bin_%d" % b)
        binG[b].SetTitle("Bin %d" % b)
        binG[b].GetXaxis().SetTitle("m_{t} [GeV]")
        binG[b].GetYaxis().SetTitle("Entries")
        binG[b].GetYaxis().SetTitleOffset(1.3)

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

    return morph


parser = ArgumentParser()
parser.add_argument("-i", dest="inDir", default="plots2018", help="Input directory containing ttrees")
parser.add_argument("-o", dest="outF", default="mtTemplatesForCH.root", help="Output template file")
parser.add_argument("-c", dest="config", default="config.txt", help="Config file specifying observables and bin ranges/sizes")
args = parser.parse_args()

if args.inDir[-1] == "/": args.inDir = args.inDir[:-1]

# Load config file
config = {}

with open(args.config, "r") as f:
    for i,line in enumerate(f):
        if i == 0 and line[0] == "#": continue
        l = line.split()
        config[l[0]] = {"nbins":int(l[1]), "min":float(l[2]), "max":float(l[3]), "title":obsTitle[l[0][4:]]}

print "\nFound observables from %s:\n" % args.config
for obs,vals in config.items():
    print "%s\t%d bins from %.1f-%.1f GeV  (bin width = %.1f GeV)" % (obs, vals["nbins"], vals["min"], vals["max"], (vals["max"] - vals["min"]) / float(vals["nbins"])) 
print "\n"

ttmasses = [1665, 1695, 1715, 1725, 1735, 1755, 1785]
tWmasses = [1695, 1725, 1755]
deltaM = 1
morph_masses = range(ttmasses[0], ttmasses[-1] + deltaM, deltaM)
print "Processing ttrees from %s" % args.inDir

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


ttmorphed = {}
tWmorphed = {}

for sys,sys_weight in systematics.items():
    # Signals
    tt[sys] = makeHist("%s/mc_TT_mt" % args.inDir, name = "ttactual", title="t#bar{t}", config=config, masses=ttmasses, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    tW[sys] = makeHist("%s/mc_ST_tW_top_mt" % args.inDir, name = "tWactual", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)
    tWantitop[sys] = makeHist("%s/mc_ST_tW_antitop_mt" % args.inDir, name = "tW", title = "ST tW", config=config, masses=tWmasses, systematic = "" if sys == "nominal" else "_"+sys, weight=sys_weight)

    # Combine tW top & antitop
    for obs,masses in tW[sys].items():
        for m in masses.keys():
            tW[sys][obs][m].Add(tWantitop[sys][obs][m])


    # Backgrounds
    DY[sys] = makeHist("%s/mc_DYJetsToLL" % args.inDir, name="DY", title="DY", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    WW[sys] = makeHist("%s/mc_WWTo2L2Nu" % args.inDir, name="WW", title="WW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    WZ[sys] = makeHist("%s/mc_WZTo3LNu" % args.inDir, name="WZ", title="WZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    ZZ[sys] = makeHist("%s/mc_ZZTo2L2Nu" % args.inDir, name="ZZ", title="ZZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    TTW[sys] = makeHist("%s/mc_TTWJetsToLNu" % args.inDir, name="TTW", title="TTW", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    TTZ[sys] = makeHist("%s/mc_TTZToLLNuNu" % args.inDir, name="TTZ", title="TTZ", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    WJets[sys] = makeHist("%s/mc_WJetsToLNu" % args.inDir, name="WJets", title="WJets", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    STs[sys] = makeHist("%s/mc_ST_s" % args.inDir, name="STs", title="ST s", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    STt[sys] = makeHist("%s/mc_ST_t_top" % args.inDir, name="STt", title="ST t", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)
    STtbar[sys] = makeHist("%s/mc_ST_t_antitop" % args.inDir, name="tbar", title="ST t antitop", config=config, systematic = "" if sys == "nominal" else "_"+sys, weight = sys_weight)

    # Combine t top & antitop
    for obs in STt[sys].keys():
        STt[sys][obs].Add(STtbar[sys][obs])



    ttmorphed[sys] = {}
    tWmorphed[sys] = {}
    for obs in config.keys():
        ttmorphed[sys][obs] = morphTemplates(tt[sys][obs], morph_masses, name=obs+"tt", title = config[obs]["title"] + " t#bar{t}", systematic = "" if sys == "nominal" else "_"+sys)
        tWmorphed[sys][obs] = morphTemplates(tW[sys][obs], morph_masses, name=obs+"tW", title = config[obs]["title"] + " tW", systematic = "" if sys == "nominal" else "_"+sys, interp="pol1")

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


outF = TFile.Open(args.outF, "recreate")
for obs in config.keys(): 
    outF.mkdir(obs)
    outF.cd(obs)
    data_obs[obs].Write()
    for sys in systematics.keys():
        for m in morph_masses:
            ttmorphed[sys][obs][m].Write(ttmorphed[sys][obs][m].GetName()[len(obs):])
            tWmorphed[sys][obs][m].Write(tWmorphed[sys][obs][m].GetName()[len(obs):])
        for m in ttmasses:
            tt[sys][obs][m].Write(tt[sys][obs][m].GetName()[len(obs):])
        for m in tWmasses:
            tW[sys][obs][m].Write(tW[sys][obs][m].GetName()[len(obs):])
        DY[sys][obs].Write(DY[sys][obs].GetName()[len(obs):])
        WW[sys][obs].Write(WW[sys][obs].GetName()[len(obs):])
        WZ[sys][obs].Write(WZ[sys][obs].GetName()[len(obs):])
        ZZ[sys][obs].Write(ZZ[sys][obs].GetName()[len(obs):])
        TTW[sys][obs].Write(TTW[sys][obs].GetName()[len(obs):])
        TTZ[sys][obs].Write(TTZ[sys][obs].GetName()[len(obs):])
        WJets[sys][obs].Write(WJets[sys][obs].GetName()[len(obs):])
        STs[sys][obs].Write(STs[sys][obs].GetName()[len(obs):])
        STt[sys][obs].Write(STt[sys][obs].GetName()[len(obs):])
outF.Close()

print "Templates saved to %s" % args.outF

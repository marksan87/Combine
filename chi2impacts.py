#!/usr/bin/env python
from ROOT import *
from array import array
from pprint import pprint
import os
import sys
from argparse import ArgumentParser
import json

parser=ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input template file")
parser.add_argument("-o", "--outF", default="chi2.json", help="chi2 impacts")
args = parser.parse_args()


separateSystSamples = ['isr','fsr','DS','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
oneSidedSysts = ["toppt", "CRerdON", "CRGluon", "CRQCD", "DS", "amcanlo", "madgraph", "herwigpp" ]

ttOnlySysts = ['toppt','hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']

observables = ["rec_ptll"]
signal = ["tt"]

tWactual = [1695,1725,1755]
mtactual = [1665,1695,1715,1725,1735,1755,1785]
mtmorphed = [1665 + i for i in xrange(121)]
masses = {}
masses["tt"] = {"actual":mtactual, "morph":mtmorphed}
masses["tW"] = {"actual":tWactual, "morph":mtmorphed}

#systematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
systematics = ["pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD']


f = TFile.Open(args.inF, 'read')
morphH = {}
systH = {}

for m in mtmorphed:
    morphH[m] = f.Get("rec_ptll/tt%d" % m).Clone()
    morphH[m].SetDirectory(0)

for syst in systematics:
    systH[syst] = {"Up":f.Get("rec_ptll/ttactual1725_%sUp" % syst).Clone()}
    systH[syst]["Up"].SetDirectory(0)
    
    if syst not in oneSidedSysts:
        systH[syst]["Down"] = f.Get("rec_ptll/ttactual1725_%sDown" % syst).Clone()
        systH[syst]["Down"].SetDirectory(0)


chi2 = {}
mt_minChi2 = {}
g = {}

N = len(mtmorphed)

for syst in systematics:
    chi2[syst] = {"Up":[]}
    mt_minChi2[syst] = {}
    g[syst] = {}
    if syst not in oneSidedSysts:
        chi2[syst]["Down"] = []

    for m in mtmorphed:
        chi2[syst]["Up"].append(morphH[m].Chi2Test(systH[syst]["Up"], "WW CHI2/NDF"))
        if syst not in oneSidedSysts:
            chi2[syst]["Down"].append(morphH[m].Chi2Test(systH[syst]["Down"], "WW CHI2/NDF"))
    
    mt_minChi2[syst]["Up"] = mtmorphed[chi2[syst]["Up"].index(min(chi2[syst]["Up"]))]
    if syst not in oneSidedSysts:
        mt_minChi2[syst]["Down"] = mtmorphed[chi2[syst]["Down"].index(min(chi2[syst]["Down"]))]

    
    # Create chi2 graphs

    g[syst]["Up"] = TGraph(N, array('d', [m/10. for m in mtmorphed]), array('d', chi2[syst]["Up"]))
    g[syst]["Up"].SetMarkerStyle(22)
    #g[syst]["Up"].Fit("pol2")
    g[syst]["Up"].SetName("%s_Up" % syst)
    g[syst]["Up"].SetTitle("%s%s" % (syst,"" if syst in oneSidedSysts else "_Up"))

    if syst not in oneSidedSysts:
        g[syst]["Down"] = TGraph(N, array('d', [m/10. for m in mtmorphed]), array('d', chi2[syst]["Down"]))
        #g[syst]["Down"].Fit("pol2")
        g[syst]["Down"].SetMarkerStyle(22)
        g[syst]["Down"].SetName("%s_Down" % syst)
        g[syst]["Down"].SetTitle("%s Down" % syst)


# Fill json file
jsondata = {"POIs":[ {"fit":[0.,0.,0.], "name":"MT"} ] }
jsondata["params"] = []

for syst in systematics:
    info = {}
    up = mt_minChi2[syst]["Up"]
    if syst not in oneSidedSysts:
        down = mt_minChi2[syst]["Down"]
    else:
        down = 1725.

    info["MT"] = [down, 1725., up]
    info["impact_MT"] = max(abs(1725 - down), abs(1725 - up))
    info["name"] = syst
    info["groups"] = []
    info["type"] = "Gaussian"
    info["fit"] = [0.,0.,0.]
    info["prefit"] = [-1., 0., 1.]

    jsondata["params"].append(info)

with open(args.outF, "w") as f:
    json.dump(jsondata, f, indent=4, separators=(',', ': '))

print "Output written to %s" % args.outF


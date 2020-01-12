#!/usr/bin/env python
from ROOT import *
from array import array
import os
import sys
import json
from argparse import ArgumentParser
from pprint import pprint
from glob import glob

systematics = ["BkgNorm", "pileup", "Lumi", "EleIDEff", "EleRecoEff", "EleScale", "EleSmear", "MuIDEff", "MuIsoEff", "MuTrackEff", "MuScale", "TrigEff", "BTagSF", "JEC", "JER", "toppt", "Q2", "Pdf", "isr", "fsr", 'hdamp','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp','DS']

parser = ArgumentParser()
parser.add_argument("-i", "--inDir", default="bin10_ptll_withbinstat_impacts")
parser.add_argument("-t", "--toyDir", default="submissionScripts/20bins_MCstatToyResults")
parser.add_argument("-o", "--outDir", default="pullPlots")
args = parser.parse_args()

inDir = args.inDir
toyDir = args.toyDir
outDir = args.outDir
if inDir[-1] == "/": inDir = inDir[:-1]
if toyDir[-1] == "/": toyDir = toyDir[:-1]
if outDir[-1] == "/": outDir = outDir[:-1]

paramsExpected = 28     # Number of impacts expected (includes stat and MCbinstats values, not used here)

impactF = "%s/%s.json" % (inDir, os.path.basename(inDir))
toyF = "%s/vals_%s.py" % (toyDir, os.path.basename(toyDir))
print "Using impact results from %s" % impactF
print "Using toy results from %s" % toyF

with open(toyF, "r") as f:
    exec(f.read())

#print "mean =", mean
#print "sigma =", sigma
#sys.exit()


with open(impactF, 'r') as f:
    jsondata = json.load(f)



toy = {}
mtvals = {"nominal":[]}
incomplete = [] 
    
# Find the MT POI
for j,p in enumerate(jsondata["POIs"]):
    if p["name"] == "MT":
        MTpoi = j
        break

    
mt = jsondata['POIs'][MTpoi]['fit'][1] / 10.
toy = {"nominal":mt}
pullCombine = {} 
pullToy = {}
# Loop over systematics



for syst in jsondata["params"]:
    name = syst["name"]
    if name == "stat" or name == "MCbinStats" or (len(name) > 2 and name[:3] == "bin"): continue
    mt = syst["MT"][1] / 10.
    toy[name] = mt
    
    pre = syst['prefit']
    fit = syst['fit']
    
    pre_err_hi = (pre[2] - pre[1])
    pre_err_lo = (pre[1] - pre[0])
    pull = fit[1] - pre[1]
    pull = (pull/pre_err_hi) if pull >= 0 else (pull/pre_err_lo)
    pull_hi = fit[2] - pre[1]
    pull_hi = (pull_hi/pre_err_hi) if pull_hi >= 0 else (pull_hi/pre_err_lo)
    pull_hi = pull_hi - pull
    pull_lo = fit[0] - pre[1]
    pull_lo = (pull_lo/pre_err_hi) if pull_lo >= 0 else (pull_lo/pre_err_lo)
    pull_lo =  pull - pull_lo
    
    pullCombine[name] = [-pull_lo, pull, pull_hi]

    impactUp = abs(fit[2] - fit[1])
    impactDown = abs(fit[1] - fit[0])

    up = (impactUp**2 + sigma[name]**2)**0.5
    dn = (impactDown**2 + sigma[name]**2)**0.5

    pullToy[name] = [-dn, toy[name] - mean[name], up]
# Pull from toys

#    pull = 
#    pullToy[syst] = 



sys.exit()

os.system("mkdir -p %s" % outDir)
gROOT.SetBatch(True)
c = TCanvas("c","c",1600,1200)

h = {}
mean = {}
sigma = {}
for syst in (["nominal"] + systematics):
    if syst not in mtvals:
        print "%s not found" % syst
        continue
    h[syst] = TH1F("MCstatToy_%s" % syst, "MC stat toys  %s" % syst, 50, 171, 174)
    N = len(mtvals[syst])
    h[syst].FillN(N, array('d', mtvals[syst]), array('d', [1]*N))
    h[syst].GetXaxis().SetTitle("m_{t} [GeV]")
    h[syst].GetXaxis().SetTitleOffset(1.3)
    h[syst].GetYaxis().SetTitle("Entries")
    h[syst].GetYaxis().SetTitleOffset(1.3)

    mean[syst] = h[syst].GetMean()
    sigma[syst] = h[syst].GetStdDev()

    h[syst].Draw("hist")
    gPad.Update()
    stats = h[syst].FindObject("stats")
    stats.SetX1NDC(0.7)
    stats.SetX2NDC(0.9)
    stats.SetY1NDC(0.7)
    stats.SetY2NDC(0.9)

    c.Modified()
    c.Update()
    c.SaveAs("%s/%s_mt.png" % (outDir,syst))

print "MC stat unc: %.2f GeV" % sigma["nominal"]

with open("%s/vals_%s.py" % (outDir,outDir), "w") as f:
    f.write("mean = ")
    pprint(mean,f)
    f.write("\nsigma = ")
    pprint(sigma,f)
    f.write("\n")


outF = TFile.Open("%s/%s.root" % (outDir,outDir), "recreate")
for syst in (["nominal"] + systematics):
    try:
        h[syst].Write()
    except KeyError:
        continue

outF.Close()




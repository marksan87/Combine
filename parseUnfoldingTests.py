#!/usr/bin/env python
from ROOT import *
import os
import sys
from argparse import ArgumentParser
from glob import glob
from pprint import pprint

observables = ["ptll","Mll","Epos","Eneg","ptpos","ptneg","Ep_Em","ptp_ptm"]

parser = ArgumentParser()
parser.add_argument("-i", "--inDir", default="unfoldingBinningTests", help="directory of unfolding impact plot data")
parser.add_argument("--obs", choices=(observables+['all']), default="ptll")
parser.add_argument("-o", "--outDir", default="resultsUnfoldingTests")
args = parser.parse_args()



if args.inDir[-1] == "/":
    args.inDir = args.inDir[:-1]

if args.outDir[-1] == "/":
    args.outDir = args.outDir[:-1]

inDir = args.inDir
outDir = args.outDir

print "Parsing results from %s" % inDir

if outDir != "":
    os.system("mkdir -p %s" % outDir)

if args.obs == "all":
    obsList = observables
else:
    obsList = [args.obs]

print "obsList =", obsList

fileList = {}
params = {}
missingParams = {}

nbins = {}
numNuisances = {}
numMCBinStats = {}
MCUnc = {}
statUnc = {}
systUnc = {}
fitUnc = {}
values = {}
for obs in obsList:
    fileList[obs] = glob("%s/%s/*/*.debug" % (inDir,obs))
    #pprint(fileList[obs])
    
    params[obs] = []
    missingParams[obs] = []
    
    nbins[obs] = {}
    numNuisances[obs] = {}
    numMCBinStats[obs] = {}
    MCUnc[obs] = {}
    statUnc[obs] = {}
    systUnc[obs] = {}
    fitUnc[obs] = {}
    values[obs] = {}

    for fName in fileList[obs]:
        p = int(os.path.basename(fName).replace("gen_unfbin","").replace("_impacts.debug",""))
        params[obs].append(p)
        with open(fName, "r") as f:
            line = f.readline()
            # Skip header line
            line = f.readline()
            values[obs][p] = line.split()

        #print values[obs][p]
    
    params[obs] = sorted(params[obs])
#    print "params =", params[obs]
#    print ""

    refParams = range(10,26)
    missingParams[obs] = range(10,26)
    for p in params[obs]:
        missingParams[obs].remove(p)


for p in params[obs]:
    # Get number of bins from the root files
    f = TFile.Open("gen_unfbin%d_mtTemplatesForCH.root" % p, "read")
    for obs in obsList:
#        print "Now getting bins for observables %s  param %d" % (obs,p)
        h = f.Get("gen_%s/ttactual1725" % obs)
        nbins[obs][p] = h.GetNbinsX()

        numNuisances[obs][p] = int(values[obs][p][0])
        numMCBinStats[obs][p] = int(values[obs][p][1])
        MCUnc[obs][p] = float(values[obs][p][2])
        statUnc[obs][p] = float(values[obs][p][3])
        systUnc[obs][p] = float(values[obs][p][4])
        fitUnc[obs][p] = float(values[obs][p][5])

    f.Close()

#pprint(nbins)



for obs in obsList:
    with open("%s/summary_%s.txt" % (outDir,obs),"w") as f:
        print "Observable: %s\n" % obs
        f.write("Observable: %s\n" % obs)
        print "param\tNbins\t27-#nuis\tnbins-#MCBins\tMCUnc\tstatUnc\tsystUnc\tfitUnc\tPassed Bin Checks"
        f.write("param\tNbins\t27-#nuis\tnbins-#MCBins\tMCUnc\tstatUnc\tsystUnc\tfitUnc\tPassed Bin Checks\n")
        print "-"*100
        f.write("-"*100 + "\n")
        for p in params[obs]:
            deltaNuisance = 27 - numNuisances[obs][p]
            deltaMCbins = nbins[obs][p] - numMCBinStats[obs][p]
            passed = "Yes" if (not deltaNuisance and not deltaMCbins) else "No"
            print "%d\t%d\t%d\t\t%d\t\t%.3f\t%.3f\t%.3f\t%.3f\t\t%s" % (p, nbins[obs][p], 27 - numNuisances[obs][p], nbins[obs][p]-numMCBinStats[obs][p], MCUnc[obs][p], statUnc[obs][p], systUnc[obs][p], fitUnc[obs][p], passed)
            f.write("%d\t%d\t%d\t\t%d\t\t%.3f\t%.3f\t%.3f\t%.3f\t\t%s\n" % (p, nbins[obs][p], 27 - numNuisances[obs][p], nbins[obs][p]-numMCBinStats[obs][p], MCUnc[obs][p], statUnc[obs][p], systUnc[obs][p], fitUnc[obs][p], passed))
        print "-"*100,"\n"
        f.write("-"*100 + "\n\n")

        if len(missingParams[obs]) > 0:
            missing = ""
            for p in missingParams[obs]:
                missing += "%d " % p
            print "Missing params:", missing 
            print "\n"

print "Results saved to %s\n" % outDir


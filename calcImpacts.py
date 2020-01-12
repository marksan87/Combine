#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
import sys
from glob import glob
from pprint import pprint
from array import array

gROOT.SetBatch(True)

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="", help="CH datacard")
parser.add_argument("-o", "--outDir", default="")#"grid_scans")
parser.add_argument("-m", "--mt", type=float, default=172.5, help="initial mt for generating asimov dataset")
parser.add_argument("--minmt", type=float, default=171.0)
parser.add_argument("--maxmt", type=float, default=174.0)
parser.add_argument("--unc", "--uncDisplay", default="both", choices=["total","syststat","both"])
parser.add_argument("--log", default="", help ="log file with output of combine commands")
parser.add_argument("--noFreeze", action="store_true", default=False, help="float other nuisances when scanning one np")
parser.add_argument("-d", "--data_obs", action="store_true", default=False, help="use data_obs instead of asimov dataset")
parser.add_argument("-v", "--verbosity", type=int, default=0, help="verbosity level")
parser.add_argument("-q", "--quiet", action="store_true", default=False, help="no terminal output (will still log with specified verbosity level")
args = parser.parse_args()

freezePars = not args.noFreeze  # Controls whether to freeze other nps when scanning one np
if not args.quiet:
    if freezePars:
        print "Other nuisances fixed to best-fit values when scanning a np"
    else:
        print "All other nuisances floated while scanning a np"


useAsimov = not args.data_obs
if not args.quiet:
    if useAsimov:
        print "Using Asimov dataset with mt = %.1f GeV" % args.mt
    else:
        print "Using data_obs"

if args.outDir == "":
    args.outDir = args.inF.replace(".txt", "").replace(".root","")
    args.outDir = args.outDir.replace("_cardMorph", "_impacts")
else:
    if args.outDir[-1] == "/":
        args.outDir = args.outDir[-1]
    #args.outDir += "/" + args.inF.replace(".txt", "").replace(".root","").replace("_cardMorph","_impacts")
    #args.outDir += os.path.basename(args.outDir.replace("_cardMorph", "_impacts"))

if args.log == "":
    args.log = os.path.basename(args.outDir).replace("_impacts",".log")
else:
    args.log = impactscan.log 

outDir = args.outDir
jsonFile = os.path.basename(outDir) + ".json" 
print "args.inF =", args.inF
print "outDir =", outDir
print "args.log =", args.log
print "jsonFile =", jsonFile
command = "mkdir -p %s; " % outDir

cardRoot = args.inF.replace(".txt",".root")
if args.inF.find(".txt") >= 0:
    # Create root workspace
    command += "text2workspace.py %s; " % args.inF
    
command += "mv %s %s/; rm -f %s/%s" % (cardRoot, outDir, outDir, args.log)
os.system(command)

if args.inF.find(".txt") >= 0:
    command = "mv %s %s %s/" % (args.inF, args.inF.replace("cardMorph.txt", "outputfileMorph.root"), outDir)
    os.system(command)
else:
    command = "mv %s %s %s/" % (args.inF.replace(".root",".txt"), args.inF.replace("cardMorph.root", "outputfileMorph.root"), outDir)
    print command
    os.system(command)


#loggingCmd = " |& tee -a %s" % args.log if not args.quiet else " &>> %s" % args.log
loggingCmd =    " &>> %s" % args.log if args.quiet else " |& tee -a %s" % args.log
silenceOutput = " &>> %s" % args.log if args.quiet else ""


f = TFile.Open("%s/%s" % (outDir, cardRoot))
w = f.Get("w")
MTvar = w.var("MT")
precision = len(str(int(MTvar.getVal()))) - 3      # Decimal places of precision

mt = int((10**precision)*args.mt)
minmt = int((10**precision)*args.minmt)
maxmt = int((10**precision)*args.maxmt)


nset = w.set("nuisances")
nuisances = []
niter = nset.createIterator()

for i in range(len(nset)):
    nuisances.append(niter.Next().getPlotLabel())

#nuisances = ['bin2']
if not args.quiet: print "Found np's:", nuisances
f.Close()
#pprint(nuisances)

# saveSpecifiedNuis causes seg fault
#options = "--saveSpecifiedNuis all --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"
#options = "--setRobustFitTolerance 0.2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"
options = "--floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --expectSignal 1"

# Initial fit
combineCmd = "combine -M MultiDimFit -n _initialFit_Test --algo singles --redefineSignalPOIs MT,r %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d -d %s -v %d --saveWorkspace %s" % (options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, cardRoot, args.verbosity, loggingCmd) 
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, silenceOutput, combineCmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

# Stat only
combineCmd = "combine -M MultiDimFit -n _paramFit_Test_stat --algo impact --redefineSignalPOIs MT,r %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d higgsCombine_initialFit_Test.MultiDimFit.mH125.root -v %d --snapshotName MultiDimFit --freezeNuisanceGroups all %s" % (options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, args.verbosity, loggingCmd) 
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, silenceOutput, combineCmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

# MC Bin stats only
combineCmd = "combine -M MultiDimFit -n _paramFit_Test_MCbinStats --algo impact --redefineSignalPOIs MT,r %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d higgsCombine_initialFit_Test.MultiDimFit.mH125.root -v %d --snapshotName MultiDimFit --freezeNuisanceGroups exp,theory %s" % (options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, args.verbosity, loggingCmd) 
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, silenceOutput, combineCmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

# Run impact scans
for n in nuisances:
    freezeArgs = ""
    if freezePars:
        # Comma separated list of all additional nps
        #freezeArgs = " --freezeParameters "
        freezeArgs = " --freezeParameters "
        for freezeMe in nuisances:
            if freezeMe != n:
                freezeArgs += "%s," % freezeMe

        # Remove trailing comma
        freezeArgs = freezeArgs[:-1]

    combineCmd = "combine -M MultiDimFit -n _paramFit_Test_%s --algo impact%s --redefineSignalPOIs MT,r -P %s %s --setParameterRanges MT=%d,%d:r=0.5,1.5 %s -m 125 --setParameters MT=%d higgsCombine_initialFit_Test.MultiDimFit.mH125.root -v %d --snapshotName MultiDimFit %s" % (n, freezeArgs, n, options, minmt, maxmt, "%s" % ("-t -1" if useAsimov else ""), mt, args.verbosity, loggingCmd)
    os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))


# Create impact json file
#jsonFile = cardRoot.replace("_cardMorph.root", "_impacts.json")
#paramList = "--named stat,"
#for n in nuisances:
#    paramList += "%s," % n
#paramList = paramList[:-1]  # Remove trailing comma
#combineCmd = "combineTool.py -M Impacts -d %s -m 125 --redefineSignalPOIs MT,r -o %s %s %s" % (cardRoot, jsonFile, paramList, loggingCmd)
combineCmd = "combineTool.py -M Impacts -d %s -m 125 --redefineSignalPOIs MT,r -o %s %s" % (cardRoot, jsonFile, loggingCmd)
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, silenceOutput, combineCmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, combineCmd, loggingCmd, combineCmd))

cmd = "%s/addStatImpact.py -j %s%s" % (os.getcwd(), jsonFile, " -q" if args.quiet else "")
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, cmd, silenceOutput, cmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, cmd, loggingCmd, cmd))


# Plot impacts
#cmd = "../createImpactPlot.py -i %s -o %s %s" % (jsonFile, jsonFile.replace(".json",""), loggingCmd)
cmd = "%s/createImpactPlot.py -i %s -o %s --mergeBinStats --unc %s %s" % (os.getcwd(), jsonFile, os.path.basename(outDir), args.unc, loggingCmd)
#os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, cmd, silenceOutput, cmd))
os.system("pushd %s%s; echo '%s'%s; %s" % (outDir, silenceOutput, cmd, loggingCmd, cmd))

print "\nImpact plot saved to %s/%s\n" % (outDir, os.path.basename(outDir)+".pdf")

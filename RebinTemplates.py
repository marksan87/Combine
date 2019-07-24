#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import os
import sys
from array import array
from pprint import pprint

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="bin1_cut220_mtTemplatesForCH.root", help="input template root file binned at 1 GeV")
parser.add_argument("--obs", default="ptll")
parser.add_argument("-b", "--bins", type=str, default="", help="List of variable bin ranges. The last entry is the upper edge of the last bin. All other entries are the lower bin edges")
parser.add_argument("-q", "--quantiles", type=int, default=0, help="rebin to quantiles")
parser.add_argument("-o", "--outF", default="rebinned_mtTemplatesForCH.root", help="output template root file")
args = parser.parse_args()

if args.bins != "":
    bins = eval(args.bins)
inF = TFile.Open(args.inF, "read")
d = inF.Get("rec_%s" % args.obs)
#bin_scaling = 1 
#binsScaled = [bin_scaling*b for b in bins]

nom = d.Get("ttactual1725").Clone("nominal")
rate = nom.Integral()

if args.quantiles > 0:
    print "Rebinning to %d quantiles" % args.quantiles
    xq = array('d', [0.] * (args.quantiles))  # Input range from [0 ... 1]
    yq = array('d', [0.] * (args.quantiles))  # Output quantiles
    for i in range(args.quantiles): xq[i] = float(i+1)/float(args.quantiles)
    nom.GetQuantiles(args.quantiles, yq, xq)

    bins = [0] + [i for i in yq]

print "Bins:", bins 
print "\nRounded:", [int(b) for b in bins]
sys.exit()

hists = []
# Loop over histograms in input file and save a list of rebinned templates
for k in d.GetListOfKeys():
    h = k.ReadObj()
    if h.ClassName() == "TH1F":
        tmp = h.Clone("_" + h.GetName())
        name = tmp.GetName()

        if name.find("rec_ptll_") >= 0:
            # Remove rec_ptll_ from the beginning of the histogram name
            name = name[len("rec_ptll_"):]
  
        if args.quantiles > 0:
            rebin = tmp.Rebin(args.quantiles, "_"+name, array('d', bins))
        else:
            #rebin = tmp.Rebin(len(bins)-1, "rebin_"+name, array('d',[10*_b for _b in bins]))
            _rebinned = tmp.Rebin(len(bins)-1, "rebin_"+name, array('d', bins))
            #rebin = TH1F(name[1:], tmp.GetTitle(), len(bins)-1, array('d', binsScaled))
            rebin = TH1F(name[1:], tmp.GetTitle(), len(bins)-1, array('d', [ 0.5 + i for i in range(len(bins)+1)])) 
            for i in xrange(1,len(bins)):
                rebin.SetBinContent(i, _rebinned.GetBinContent(i))
                rebin.SetBinError(i, _rebinned.GetBinError(i))

        rebin.SetDirectory(0)

        # Scale by bin width
        #rebin.Scale(1., "width")
        
        # Renormalize to nominal rate
       # rebin.Scale(rate/rebin.Integral())

        hists.append(rebin)

inF.Close()

outF = TFile.Open(args.outF, "recreate")
outF.mkdir("rec_ptll")
outF.cd("rec_ptll")

for h in hists:
    #h.SetName(h.GetName()[len("rebin_"):])
    if args.quantiles > 0:
        h.SetName(h.GetName()[2:])
#        h.Write(h.GetName()[2:])
#    else:
#        h.Write()
    h.Write()

outF.Close()


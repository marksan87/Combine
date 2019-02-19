#!/usr/bin/env python
from ROOT import *
import numpy as np
from argparse import ArgumentParser
import gzip
import pickle


parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="mtTemplatesForCH.root", help="input template root file")
parser.add_argument("-o", "--outF", default="tensorflowTemplates.pklz", help="output datacard and root file name")
parser.add_argument("--obs", default="ptll", help="observable")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reco level")
args = parser.parse_args()

f = TFile.Open(args.inF, 'read')

tWmasses = [1695, 1725, 1755]
ttmasses = [1665,1695,1715,1725,1735,1755,1785]
mtmorph = [1665 + i for i in range(121)]

masses = {}
masses["tt"] = {"actual":ttmasses, "morph":mtmorph}
masses["tW"] = {"actual":tWmasses, "morph":mtmorph}
signal = ["tt","tW"]

templates = {}

print "Loading templates from %s" % args.inF

for s in signal:
    templates[s] = {"actual":{}, "morph":{}}
    for mode in ["actual", "morph"]:
        for m in masses[s][mode]:
            h = f.Get("rec_ptll/%s%s%d" % (s,"actual" if mode == "actual" else "",m)).Clone()
            h.SetDirectory(0)
            h.Scale(1./h.Integral())
            nbins = h.GetNbinsX()
            mass = np.array([m/10. for i in range(nbins)], dtype=np.float64)

            binCenters = []  # Ptll
            binContents = [] # Event count
            binErrors = []   # Error in events
            
            for b in range(1,nbins+1):
                binCenters.append(h.GetBinCenter(b))
                binContents.append(h.GetBinContent(b))
                binErrors.append(h.GetBinError(b))

            ptll = np.array(binCenters, dtype=np.float64)
            events = np.array(binContents, dtype=np.float64)
            errors = np.array(binErrors, dtype=np.float64)

            templates[s][mode][m] = {u"ptll":ptll, u"events":events, u"errors":errors}



with gzip.open(args.outF, "wb") as f:
    pickle.dump(templates["tt"], f, pickle.HIGHEST_PROTOCOL)
    pickle.dump(templates["tW"], f, pickle.HIGHEST_PROTOCOL)

print "Template info saved to %s" % args.outF

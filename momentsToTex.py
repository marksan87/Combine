#!/usr/bin/env python
import os
import sys
import pickle
import gzip
import json
from argparse import ArgumentParser
from pprint import pprint

obsTitle = {"ptll":"\\ptll",
            "ptpos":"\\ptpos",
            "Epos":"\\Epos",
            "ptp_ptm":"\\ptpptm",
            "Ep_Em":"\\EpEm",
            "Mll":"\\Mll",
    }

allobservables = ["ptll","Mll","ptpos","Epos","ptp_ptm","Ep_Em"]

parser = ArgumentParser()
parser.add_argument("-i", "--inDir", default="all_rbf_plots")
parser.add_argument("-s", "--sig", dest="signal", default="tttW", choices=["tt","tW","tttW"])
parser.add_argument("--reco", default="rec", choices=["rec","gen"])
parser.add_argument("--obs", default=["all"], nargs="+", choices=(["all"]+allobservables))
parser.add_argument("-o", "--outF", default="momentTable.tex")
args = parser.parse_args()

signal = args.signal 
momentRange = [1,2]

inDir = args.inDir
if inDir[-1] == "/": inDir = inDir[:-1]
outF = args.outF

reco = args.reco

if "all" in args.obs:
    observables = allobservables
else:
    observables = args.obs

# Load moment pickle file
with gzip.open("%s/moments.pklz" % inDir) as f:
    momentData = pickle.load(f)

moments = {}
momentErrors = {}

for obs in observables:
    moments[obs] = {}
    momentErrors[obs] = {}
    for mom in momentRange:
        centralN = (len(momentData[signal]["actual"]["%s_%s"%(reco,obs)]["nominal"]["m%d"%mom])-1)/2
        moments[obs][mom] = momentData[signal]["actual"]["%s_%s"%(reco,obs)]["nominal"]["m%d"%mom][centralN]
        momentErrors[obs][mom] = momentData[signal]["actual"]["%s_%s"%(reco,obs)]["nominal"]["m%derr"%mom][centralN]


mtFit = {}
mtFitError = {}

for obs in observables:
    mtFit[obs] = {}
    mtFitError[obs] = {}
    for mom in [1,2]:
    # Load json file with moment masses
        jsonFileName = "%s/%s_%s/%s_%s_%s_m%d_uncalibrated_impacts.json" % (inDir,reco,obs,signal,reco,obs,mom)
        with open(jsonFileName) as jsonF:
            data = json.load(jsonF)


        POIs = [ele['name'] for ele in data['POIs']]
        selectedPOI = -1
        for i,p in enumerate(POIs):
            if str(p) == str("MT"):
                selectedPOI = i
                break
        if selectedPOI < 0:
            # Couldn't find the POI
            print "Couldn't find POI %s in %s! Defaulting to first POI" % ("MT",jsonFileName)
            selectedPOI = 0


        # Decimal points of precision, determined by number of whole number digits after the first 3
        precision = len(str(int(data[u'POIs'][selectedPOI][u'fit'][0]))) - 3
        decimalScaling = 10**precision
        
        try:
            for p in range(len(data[u'params'])):
                #print "Now on %s" % data[u'params'][p]['name']
                if str(data[u'params'][p]['name']) == "stat":
                    mtFit[obs][mom] = data[u'params'][p]["MT"][1]/10.
                    mtFitError[obs][mom] = (data[u'params'][p]["MT"][1] - data[u'params'][p]["MT"][0])/10.
                break
            if mom not in mtFit[obs] or mom not in mtFitError[obs]:
                raise KeyError

        except KeyError:
            print "Can't find stat uncertainty in %s!" % jsonFileName


# Write table
f = open(outF, "w")

f.write(\
"""\\begin{table*}[h]
\t\\begin{center}
\t\t\\footnotesize
\t\t\\caption{First two moments and calibrated masses of reco-level observables fit to nominal \\ttbar+\\tW simulation at \\mt = 172.5 GeV. Uncertainties in extracted masses are statistical.}
\t\t\\label{table:rec_moments}
\t\t\\begin{tabular}{l | c c | c c }
\t\t\t\\hline
\t\t\tObservable & Moment 1 [GeV] & \\mt [GeV] & Moment 2 [GeV\\textsuperscript{2}] & \\mt [GeV] \\\\
\t\t\t\\hline
""")

for obs in observables:
    f.write("\t\t\t%s & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ \\\\\n" % (obsTitle[obs], moments[obs][1], momentErrors[obs][1], mtFit[obs][1], mtFitError[obs][1], moments[obs][2], momentErrors[obs][2], mtFit[obs][2], mtFitError[obs][2]) ) 

f.write("""\t\t\t\\hline\n""")

f.write(
"""\t\t\\end{tabular}
\t\\end{center}
\\end{table*}
"""
    )

f.close()

print "Impact table written to %s" % outF



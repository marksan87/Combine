#!/usr/bin/env python
# Adds stat uncertainty to json file output from combineTool.py -M Impacts
import json
from argparse import ArgumentParser
from ROOT import TFile 

parser = ArgumentParser()
parser.add_argument("-j", "--json", help="input json file")
parser.add_argument("-s", "--statF", default="higgsCombine_paramFit_Test_stat.MultiDimFit.mH125.root")
parser.add_argument("--addBinStats", action="store_true", default=False, help="add barlow-beeston lite binwise stat unc nps")
parser.add_argument("-b", "--binStatF", default="higgsCombine_paramFit_Test_MCbinStats.MultiDimFit.mH125.root")
parser.add_argument("-q", "--quiet", action="store_true", default=False, help="run silently")
parser.add_argument("-o", "--outF", default="", help="store results in this output file (or overwrite original json file if left empty")
args = parser.parse_args()

addBinStats = args.addBinStats
with open(args.json) as jsonfile:
    data = json.load(jsonfile)

f = TFile.Open(args.statF)
t = f.Get("limit")

mtvals = [t.MT for evt in t]
mt = [mtvals[1], mtvals[0], mtvals[2]]  # [ -sigma, nominal, +sigma ]
impact = max(abs(mt[2] - mt[1]), abs(mt[0] - mt[1]))


#print "Stat uncertainty: %.2f GeV" % (impact/10.)
statPosition = -1  # Position of stat parameter (if it already exists) 
for p in xrange(len(data[u'params'])):
    if data[u'params'][p][u'name'] == u'stat':
        statPosition = p

statVals = \
{
    u'name': u'stat',
    u'MT': mt,
    u'impact_MT': impact,
    u'impact_r': 0.0,
    u'prefit': [-1.0, 0.0, 1.0],
    u'fit': [ -1.0, 0.0, 1.0],
    u'groups': [],
    u'r': [1.0, 1.0, 1.0],
    u'type': "Gaussian",
}

if statPosition >= 0:
#    if not args.quiet: print "Replacing stat uncertainty values in json file"
    data[u'params'][statPosition] = statVals
else:
#    if not args.quiet: print "Adding stat uncertainty values to json file"
    data[u'params'].append(statVals)

if addBinStats:
    f = TFile.Open(args.binStatF)
    t = f.Get("limit")

    MCmtvals = [t.MT for evt in t]
    MCmt = [MCmtvals[1], MCmtvals[0], MCmtvals[2]]  # [ -sigma, nominal, +sigma ]
    MCimpact = max(abs(MCmt[2] - MCmt[1]), abs(MCmt[0] - MCmt[1]))

    MCstatPosition = -1  # Position of MC stat parameter (if it already exists) 
    for p in xrange(len(data[u'params'])):
        if data[u'params'][p][u'name'] == u'MCbinStats':
            statPosition = p

    MCstatVals = \
    {
        u'name': u'MCbinStats',
        u'MT': MCmt,
        u'impact_MT': MCimpact,
        u'impact_r': 0.0,
        u'prefit': [-1.0, 0.0, 1.0],
        u'fit': [ -1.0, 0.0, 1.0],
        u'groups': [],
        u'r': [1.0, 1.0, 1.0],
        u'type': "Gaussian",
    }

    if MCstatPosition >= 0:
#    if not args.quiet: print "Replacing stat uncertainty values in json file"
        data[u'params'][MCstatPosition] = MCstatVals
    else:
#    if not args.quiet: print "Adding stat uncertainty values to json file"
        data[u'params'].append(MCstatVals)



jsondata = json.dumps(data, sort_keys=True, indent=2, separators=(',', ': '))

outF = args.json if args.outF == "" else args.outF

with open(outF, "w") as jsonfile:
    jsonfile.write(jsondata)

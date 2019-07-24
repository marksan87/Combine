#!/usr/bin/env python
from ROOT import *
from argparse import ArgumentParser
import json
import sys

obsTitle = {"ptll":"p_{T}(ll)",
            "ptpos":"p_{T}(l^{+})",
            "Epos":"E(l^{+})",
            "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})",
            "Ep_Em":"E(l^{+}) + E(l^{-})",
            "Mll":"M(ll)",
            "leadLepPt":"Leading Lepton p_{T}",
            "leadJetPt":"Leading Jet p_{T}",
            "ptll_M0_E0":"p_{T}(ll)  M(ll) < 80  E(l^{+}) < 60",
            "ptll_M0_E1":"p_{T}(ll)  M(ll) < 80  60 <= E(l^{+}) < 100",
            "ptll_M0_E2":"p_{T}(ll)  M(ll) < 80  100 <= E(l^{+}) < 600",
            "ptll_M1_E0":"p_{T}(ll)  80 <= M(ll) < 130  E(l^{+}) < 60",
            "ptll_M1_E1":"p_{T}(ll)  80 <= M(ll) < 130  60 <= E(l^{+}) < 100",
            "ptll_M1_E2":"p_{T}(ll)  80 <= M(ll) < 130  100 <= E(l^{+}) < 600",
            "ptll_M2_E0":"p_{T}(ll)  130 <= M(ll) < 600  E(l^{+}) < 60",
            "ptll_M2_E1":"p_{T}(ll)  130 <= M(ll) < 600  60 <= E(l^{+}) < 100",
            "ptll_M2_E2":"p_{T}(ll)  130 <= M(ll) < 600  100 <= E(l^{+}) < 600",
            }


parser = ArgumentParser()
parser.add_argument("-1", "--inF1", help="input json file with bin stats impacts")
parser.add_argument("-2", "--inF2", help="input json file with bin stats impacts")
parser.add_argument("--obs", default="ptll")
parser.add_argument("--reco", default="rec")
parser.add_argument("--allSysts", action="store_true", default=False, help="all systematics included in the fits (for table text only)")
parser.add_argument("-o", "--outF", default="", help="output plot file")
args = parser.parse_args()

recoObs = "%s_%s" % (args.reco, args.obs)
if args.outF == "":
    args.outF = "MC_stats_compare_impacts.tex"

print "Comparing impacts from %s and %s" % (args.inF1, args.inF2)

with open(args.inF1) as jf:
    data1 = json.load(jf)

with open(args.inF2) as jf:
    data2 = json.load(jf)

binning = [0, 17, 25, 32, 37, 42, 47, 52, 56, 60, 64, 68, 72, 76, 81, 86, 92, 98, 105, 114, 126]
bins1 = {}
bins2 = {}

bins1Total = bins2Total = 0.

for p in xrange(len(data1[u'params'])):
    if "bin" in data1[u'params'][p]['name']:
        b = int( data1[u'params'][p]['name'][ data1[u'params'][p]['name'].find("bin")+3:] )
        bins1[b] = data1[u'params'][p][u'impact_MT'] / 10.
        bins1Total += bins1[b]**2
bins1Total = bins1Total**0.5

for p in xrange(len(data2[u'params'])):
    if "bin" in data2[u'params'][p]['name']:
        b = int(data2[u'params'][p]['name'][ data2[u'params'][p]['name'].find("bin")+3:] )
        bins2[b] = data2[u'params'][p][u'impact_MT'] / 10.
        bins2Total += bins2[b]**2
bins2Total = bins2Total**0.5

print "File 1 has %d bins and total MC stat unc: %.2f" % (len(bins1),bins1Total)
print "File 2 has %d bins and total MC stat unc: %.2f" % (len(bins2),bins2Total)

f = open(args.outF, "w")
f.write( \
"""\\documentclass{article}

\\begin{document}

\\begin{table*}[h]
\t\\begin{center}
\t\t\\caption{Comparison of MC stat impacts for %s $%s$ with 10 and 20 bins %s}
\t\t\\label{table:cutflow}
\t\t\\begin{tabular}{ c | c c | c }
\t\t\t\\hline
\t\t\t$%s$ [GeV] & 20 Bins & 20 Bins Total & 10 Bins \\\\
\t\t\t\\hline
""" % (args.reco, obsTitle[args.obs], "(all systematics included in fits)" if args.allSysts else "(no other systematics included)", obsTitle[args.obs])
)

for b in xrange(1, len(bins2)+1):
    f.write("\t\t\t%s & %s & %s & %s \\\\\n" % ("%d - %d"%(binning[b-1],binning[b]), "%.2f" % bins2[b], "" if b % 2 == 1 else "%.2f" % ((bins2[b-1]**2+bins2[b]**2)**0.5), "" if b % 2 == 1 else "%.2f" % bins1[b/2] ) )
    if b % 2 == 0:
        f.write("\t\t\t\\hline\n")

f.write(
"""\t\t\tTotal MC unc & & %.2f & %.2f \\\\
\t\t\t\\hline
\t\t\\end{tabular}
\t\\end{center}
\\end{table*}

\\end{document}
""" % (bins2Total, bins1Total)
)

f.close()

print "Output saved to %s" % args.outF

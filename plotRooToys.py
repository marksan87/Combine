#!/usr/bin/env python
import ROOT
from ROOT import *
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--inF", nargs="+", help="input toy root file(s)")
parser.add_argument("-t", "--toy", type=int, default=1, help="plot this toy", metavar="N")

args = parser.parse_args()

files = {}
h = {}

#gROOT.ProcessLine("""
##include <TFile.h>
#TFile* f = new TFile('%s');
#f->ls();
#
#""" % ("%s/%s" % (os.getcwd(),args.inF[0]))
#)

gROOT.ProcessLine("""
TFile* f1 = TFile::Open("/uscms_data/d3/msaunder/combine/CMSSW_8_1_0/src/UserCode/mtMorphing/submissionScripts/toy1_higgsCombineTest.GenerateOnly.mH1725.1.root");
f1->cd("toys");
TH1F* h1 = toy_1->createHistogram("CMS_th1x");
h1->SetDirectory(0);
f1->Close();

TFile* f2 = TFile::Open("/uscms_data/d3/msaunder/combine/CMSSW_8_1_0/src/UserCode/mtMorphing/submissionScripts/toy2_higgsCombineTest.GenerateOnly.mH1725.1.root");
f2->cd("toys");
TH1F* h2 = toy_1->createHistogram("CMS_th1x");
h2->SetDirectory(0);
f2->Close();

"""
)
from ROOT import h1,h2
h1.Draw("hist")

diff=h1.Clone()
diff.Add(h2,-1)
diff.Draw("hist")

#for i,inF in enumerate(args.inF):
#    f = TFile.Open(inF, "read")
#    h[i] = f.Get("toys/toy_%d" % args.toy).createHistogram("CMS_th1x")
#    h[i].SetDirectory(0)
#    f.Close()
#
#
#h[0].SetLineColor(kRed)
#h[1].SetLineColor(kBlue)
#
#diff = h[0].Clone()
#diff.Add(h[1], -1)
#
#h[0].Draw("hist")
#h[1].Draw("hist same")


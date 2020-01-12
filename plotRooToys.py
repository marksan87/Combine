#!/usr/bin/env python
import ROOT
from ROOT import *
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--inF", default="submissionScripts/interactive/higgsCombine_paramFit_Test_bin1.MultiDimFit.mH125.123456.root", help="input toy root file(s)")
parser.add_argument("-t", "--toy", type=int, default=1, help="plot this toy", metavar="N")
parser.add_argument("--orig", default="asimov_ptll_q20_mtTemplatesForCH.root")
args = parser.parse_args()

files = {}
h = {}

origF = TFile.Open(args.orig)
horig = origF.Get("rec_ptll/data_obs")

#gROOT.ProcessLine("""
##include <TFile.h>
#TFile* f = new TFile('%s');
#f->ls();
#
#""" % ("%s/%s" % (os.getcwd(),args.inF[0]))
#)

gROOT.ProcessLine("""
TFile* f = TFile::Open("%s");
f->cd("toys");
TH1F* h1 = (TH1F*)toy_1->createHistogram("CMS_th1x");
h1->SetDirectory(0);
TH1F* h2 = (TH1F*)toy_2->createHistogram("CMS_th1x");
h2->SetDirectory(0);
f->Close();

""" % args.inF
)
from ROOT import h1,h2
#h1.Draw("hist")

h1r = horig.Clone("toy1")
h1r.Reset()
for b in xrange(1, h1r.GetNbinsX()+1):
    h1r.SetBinContent(b, h1.GetBinContent(b))
    h1r.SetBinError(b, h1.GetBinError(b))

#diff=h1.Clone("h2_diff")
#diff.Add(h2,-1)
diff=h1r.Clone("h2_diff")
diff.Add(horig,-1)
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


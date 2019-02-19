#!/usr/bin/env python
from __future__ import print_function
from ROOT import *
import numpy as np
import os
import sys
import gzip
import pickle
import array
from pprint import pprint
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument("-i", "--inF", default="eventPredictionsFromTF.pklz", help="input prediction gziped pickle file")
parser.add_argument("-o", "--outDir", default="morphTF", help="output file")
args = parser.parse_args()

with gzip.open(args.inF, "rb") as f:
    events = pickle.load(f)
    masses = pickle.load(f)
    ptll = pickle.load(f)

# Reconstruct templates
nbins = len(ptll)
binW = ptll[1]-ptll[0]
binMin = ptll[ 0] - binW / 2.
binMax = ptll[-1] + binW / 2.

masses10 = [int(m*10) for m in masses]
contents = {}
for i,m in enumerate(masses10):
    contents[m] = events[nbins*i : nbins*(i+1)]

h = {}
for m in masses10:
    h[m] = TH1F("morph_tt%d" % m, "t#bar{t} DNN Morph  m_{t} = %.1f" % (m/10.), nbins, binMin, binMax)
    h[m].SetDirectory(0)

    h[m].FillN(nbins, array.array('d',ptll), array.array('d',contents[m]))
    
h[1725].Draw("hist")

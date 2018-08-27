#!/usr/bin/env python
from ROOT import TFile, TH1F

masses = [1665,
          1695,
          1715,
          1725,
          1735,
          1755,
          1785]
morphed_masses = range(1665, 1786, 1)
#obs = "Ep_Em"
obs = "pt_ll"
tt_h = {}

f = TFile.Open("morphed_templates.root", "read")
for m in morphed_masses:
    tt_h[m] = f.Get("morph_mt%d" % m).Clone("signal%d" % m)
    tt_h[m].SetDirectory(0)
f.Close()
#for m in masses:
#    ttF = TFile.Open("RootFiles/mc_TT_mt%i.root"%m,"read")
#    tt_h[m] = ttF.Get(obs).Clone("signal%d" % m)
#    tt_h[m].SetDirectory(0)
#    ttF.Close()


# Use tW for background, for now
tWFile = TFile("RootFiles/mc_ST_tW_top_mt1725.root","read")
#tWFile = TFile("RootFiles/mc_DY.root","read")
bkg = tWFile.Get(obs).Clone("bkg")


data_signal_comp = 1.1
data_bkg1_comp = 0.9

nBins = bkg.GetNbinsX()

#data_Hist = TH1F("data","data",nBins,0,nBins)
#data_Hist.Add(tt_h[1735], 0.5)
data_Hist = tt_h[1735].Clone("data_obs")
#data_Hist = tt_h[1725].Clone("data_obs")
#data_Hist.Scale(1.0/3.0)
#data_Hist.Add(tt_h[1695], 1.0/3.0)
#data_Hist.Add(tt_h[1735], 1.0/3.0)

#data_Hist.Scale(0.5)
#data_Hist.Add(tt_h[1755], 0.5)
data_Hist.Add(bkg, 1.)


outF = "mtTemplatesForCH"
print "Creating output file %s.root..."%outF,
outputFile = TFile("%s.root"%outF,"recreate")
outputFile.mkdir("emu")
outputFile.cd("emu")
data_Hist.Write("data_obs")

for m in morphed_masses:
    tt_h[m].Write()



bkg.Write("bkg")

outputFile.Close()

print "done."

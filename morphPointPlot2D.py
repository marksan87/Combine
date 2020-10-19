#!/usr/bin/env python
import sys
from ROOT import * 
from argparse import ArgumentParser
from array import array

obsTitle = {"ptll":"p_{T}(ll)", "ptpos":"p_{T}(l^{+})", "Epos":"E(l^{+})", "ptp_ptm":"p_{T}(l^{+}) + p_{T}(l^{-})", "Ep_Em":"E(l^{+}) + E(l^{-})", "Mll":"M(ll)"}
signalTitle = {"tt":"t#bar{t}", "tW":"tW", "tttW":"t#bar{t} + tW"}

parser = ArgumentParser()
parser.add_argument("-i", dest="inF", default="mtTemplatesForCH.root", help="Input debug root file from CombineHarvester with morphed templates")
parser.add_argument("-a", "--actual", action="store_true", default=False, help="use actual templates instead of morphed")
#parser.add_argument("-s", "--sig", nargs="+", default=["tt", "tW"], choices=["tt","tW", "tttW"], help="signal sample")
parser.add_argument("-s", "--sig", nargs="+", default=["tt"], choices=["tt","tW", "tttW"], help="signal sample")
parser.add_argument("--syst", default="", help="systematic")
parser.add_argument("--reco", default="rec", choices=["rec","gen"], help="reconstruction level")
parser.add_argument("--obs", dest="obs", default="ptll", help="Kinematic observable")
parser.add_argument("--ratio", action="store_true", default=False, help="plot ratio wrt the nominal template instead")
parser.add_argument("-n", "--norm", action="store_true", default=False, help="normalize all templates to unity")
parser.add_argument("--includeDebug", action="store_true", default=False, help="store templates used in 2D hist/graph in output root file")
parser.add_argument("-o", "--outF", default="", help="Output root file")
args = parser.parse_args()

print "Using %s %s templates" % ("actual" if args.actual else "morphed", args.sig)

useRatio = args.ratio
if useRatio: print "Plotting ratio to nominal"
normalize = args.norm

if args.obs not in obsTitle.keys():
    print "Invalid observable:", obs
    print "Choose from:"
    print obsTitle.keys()
    sys.exit()

if args.outF == "":
    if args.inF.find("mtTemplatesForCH.root") >= 0:
        args.outF = args.inF.replace("mtTemplatesForCH.root", "morph2D.root")    
    else:
        args.outF = "morph2D.root"

f = TFile.Open(args.inF)
deltaM = 1

if args.actual:
    masses = { \
        "tt":   [1665, 1695, 1715, 1725, 1735, 1755, 1785],
        "tW":   [1695, 1725, 1755],
        "tttW": [1695, 1725, 1755]
    }
    mtbins = { \
        "tt":     [165.0, 168.0, 171.0, 172.0, 173.0, 174.0, 177.0, 180.0],
        "tW":     [168.0, 171.0, 174.0, 177.0],
        "tttW":   [168.0, 171.0, 174.0, 177.0],
    }
else:
    morphMasses = range(1665, 1786, deltaM)
    masses = {"tt":morphMasses, "tW":morphMasses, "tttW":morphMasses}
    massBins = [(m/10. - deltaM/20.0) for m in morphMasses]
    massBins.append(morphMasses[-1]/10. + deltaM/20.0)
    massBins = [round(m, 5) for m in massBins]

#signal = ["tt","tW"]
signal = args.sig 
#systematics = ["nominal", "pileupUp", "pileupDown", "Q2Up", "Q2Down", "PdfUp", "PdfDown"]
#systematics = ["nominal", "bin1Up", "bin1Down", "bin2Up", "bin2Down"]

if args.syst == "nominal":
    systematics = ["nominal"]
else:
    systematics = ["nominal", args.syst]

systematics = ["nominal", 'BTagSFDown', 'BTagSFUp', 'CRGluonDown', 'CRGluonUp', 'CRQCDDown', 'CRQCDUp', 'CRerdONDown', 'CRerdONUp', 'DSDown', 'DSUp', 'EleIDEffDown', 'EleIDEffUp', 'EleRecoEffDown', 'EleRecoEffUp', 'EleScaleDown', 'EleScaleUp', 'EleSmearDown', 'EleSmearUp', 'JECDown', 'JECUp', 'JERDown', 'JERUp', 'LumiDown', 'LumiUp', 'MuIDEffDown', 'MuIDEffUp', 'MuIsoEffDown', 'MuIsoEffUp', 'MuScaleDown', 'MuScaleUp', 'MuTrackEffDown', 'MuTrackEffUp', 'PdfDown', 'PdfUp', 'Q2Down', 'Q2Up', 'TrigEffDown', 'TrigEffUp', 'UEDown', 'UEUp', 'amcanloDown', 'amcanloUp', 'fsrDown', 'fsrUp', 'hdampDown', 'hdampUp', 'herwigppDown', 'herwigppUp', 'isrDown', 'isrUp', 'madgraphDown', 'madgraphUp', 'pileupDown', 'pileupUp', 'topptDown', 'topptUp']

ttOnlySysts = ['toppt','Pdf','UE','CRerdON','CRGluon','CRQCD','amcanlo','madgraph','herwigpp']
tWOnlySysts = ['DS']


recoObs = "%s_%s" % (args.reco, args.obs)
print "Creating 2D plots for %s%s" % ("normalized " if normalize else "", recoObs)
print "Signals: %s" % signal
print "Systematics: %s" % systematics


h = {}
g2D = {}
morph = {}
nominalMorph = {}   # Clone of nominal templates before scaling for ratio
for s in signal:
    h[s] = {}
    g2D[s] = {}
    morph[s] = {}
    nominalMorph[s] = {}
    m10 = [m/10. for m in masses[s]]
    mbinavg = []
    for i in xrange(len(m10)-1):
        mbinavg.append((m10[i] + m10[i+1])/2.)

    if args.actual:
        massBins = mtbins[s]
    #print "%s mbinavg:" % s, mbinavg
    #massBins = [2*m10bins[0] - mbinavg[0]] + mbinavg + [2*m10bins[-1] - mbinavg[-1]]
    print "%s: " % s, massBins
   
    for syst in systematics:
        if syst[-2:] == "Up":
            systType = syst[:-2]
        elif syst[-4:] == "Down":
            systType = syst[:-4]
        else:
            # Nominal
            systType = syst
        if (s == "tt" and systType in tWOnlySysts) or (s == "tW" and systType in ttOnlySysts): continue

        #print "Now on %s" % syst
        morph[s][syst] = {}
        for m in masses[s]:
            morph[s][syst][m] = f.Get("%s/%s%s%d%s" % (recoObs, s, "actual" if args.actual else "", m, "" if syst=="nominal" else "_" + syst))
            morph[s][syst][m].SetDirectory(0)
            morph[s][syst][m].GetXaxis().SetTitle("%s [GeV]" % obsTitle[args.obs])
            morph[s][syst][m].GetXaxis().SetTitleOffset(1.2)
            morph[s][syst][m].GetYaxis().SetTitle("Entries")
            morph[s][syst][m].GetYaxis().SetTitleOffset(1.2)
           
            if normalize:
                morph[s][syst][m].Scale(1/morph[s][syst][m].Integral())

            if useRatio and syst == "nominal":
                nominalMorph[s][m] = morph[s]["nominal"][m].Clone(morph[s]["nominal"][m].GetName()+"_nomrate")

        if useRatio:
            if syst == "nominal":
                #nom = morph[s]["nominal"][1725].Clone("nom")
                nom = nominalMorph[s][m] 
            for m in masses[s]:
                if syst != "nominal":
                    nom = nominalMorph[s][m] 

                morph[s][syst][m].Divide(nom)
                morph[s][syst][m].SetTitle(morph[s][syst][m].GetTitle() + " ratio to nominal")

        uniformBinning = True
        obsBins = morph[s][syst][masses[s][0]].GetXaxis().GetXbins()
        obsBins = [p for p in obsBins]
        obsBinCenters = array('d', [0]*(len(obsBins)-1))
        morph[s][syst][masses[s][0]].GetXaxis().GetCenter(obsBinCenters)
        obsBinCenters = [i for i in obsBinCenters]
        diff = []
        for i in xrange(len(obsBins)-1):
            diff.append(obsBins[i+1] - obsBins[i])
       
        bin1 = diff[0]
        for i in xrange(1,len(diff)):
            if abs(diff[i] - bin1) > 1e-4:
                # Variable binning detected
                uniformBinning = False
                break

        if uniformBinning: 
            # Uniform binning
            print "Uniform binning detected"
            uniformBinWidth = morph[s][syst][masses[s][0]].GetBinLowEdge(2) - morph[s][syst][masses[s][0]].GetBinLowEdge(1)
            obsBins = []
            for _b in xrange(1, morph[s][syst][masses[s][0]].GetNbinsX()+2):
                obsBins.append(morph[s][syst][masses[s][0]].GetBinLowEdge(_b))

        else:
            print "Variable binning detected"

       # massPoints2D = []
       # for m in m10:
       #     massPoints2D += [m]*len(obsBinCenters)
        #for m in masses[s]:


        g2D[s][syst] = TGraph2DErrors()
        g2D[s][syst].SetName("%s_%s%s_graph2D" % (s,args.obs, "" if syst=="nominal" else "_" + syst))

        h[s][syst] = TH2F("%s_%s%s_2D" % (s,args.obs,"" if syst=="nominal" else "_" + syst), "%s%s %s%s Templates" % (signalTitle[s], "" if args.actual else " Morphed", obsTitle[args.obs],"" if syst == "nominal" else "_" + syst), len(massBins)-1, array('d',massBins), morph[s][syst][masses[s][0]].GetNbinsX(), array('d', obsBins) )
        h[s][syst].GetXaxis().SetTitle("m_{t} [GeV]")
        h[s][syst].GetYaxis().SetTitle("%s [GeV]" % obsTitle[args.obs])
        h[s][syst].GetZaxis().SetTitle("Events / %sGeV" % ("%s " % uniformBinWidth if uniformBinning else ""))

        h[s][syst].GetXaxis().SetTitleOffset(1.6)
        h[s][syst].GetYaxis().SetTitleOffset(1.8)
        h[s][syst].GetZaxis().SetTitleOffset(1.5)

        for m in masses[s]:
            for b in xrange(1, morph[s][syst][m].GetNbinsX() + 1):
                if uniformBinning:
                    h[s][syst].Fill(m / 10.0, morph[s][syst][m].GetBinCenter(b), morph[s][syst][m].GetBinContent(b))
                    point = g2D[s][syst].GetN()
                    g2D[s][syst].SetPoint(point, m/10., morph[s][syst][m].GetBinCenter(b),  morph[s][syst][m].GetBinContent(b))
                    g2D[s][syst].SetPointError(point, 0., 0., morph[s][syst][m].GetBinError(b))
                else:
                    h[s][syst].Fill(m / 10.0, morph[s][syst][m].GetBinCenter(b), morph[s][syst][m].GetBinContent(b) / morph[s][syst][m].GetXaxis().GetBinWidth(b))
                    point = g2D[s][syst].GetN()
                    g2D[s][syst].SetPoint(point, m/10., morph[s][syst][m].GetBinCenter(b),  morph[s][syst][m].GetBinContent(b) / morph[s][syst][m].GetXaxis().GetBinWidth(b))
                    g2D[s][syst].SetPointError(point, 0., 0., morph[s][syst][m].GetBinError(b) / morph[s][syst][m].GetXaxis().GetBinWidth(b))


outF = TFile.Open(args.outF, "recreate")
for s in signal:
    for syst in systematics:
        if syst[-2:] == "Up":
            systType = syst[:-2]
        elif syst[-4:] == "Down":
            systType = syst[:-4]
        else:
            # Nominal
            systType = syst
        if (s == "tt" and systType in tWOnlySysts) or (s == "tW" and systType in ttOnlySysts): continue
        h[s][syst].Write()
        g2D[s][syst].Write()
        if args.includeDebug:
            for m in masses[s]:
                #morph[s][syst][m].Write()
                morph[s][syst][m].Write("%s_%s_%d%s" % (s,args.obs,m, "" if syst == "nominal" else "_%s" % syst))

outF.Close()

print "Results saved in %s" % args.outF

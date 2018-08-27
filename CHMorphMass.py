#!/usr/bin/env python
import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT

cb = ch.CombineHarvester()
cb.SetVerbosity(3)

"""
masses = [  \
          '1665',
          '1695',
          '1715',
          '1725',
          '1735',
          '1755',
          '1785' \
          ]
"""
deltaM = 1
masses = [str(1665 + i * deltaM) for i in xrange(121)]
#masses = [str(1695 + i * deltaM) for i in xrange(61)]

systematics = {"pileup":"shape"}#, "toppt":"shape"}
signals = ['tt','tW']
cats = [(0,'rec_ptll'), \
        (1,'rec_ptpos'), \
        (2,'rec_ptp_ptm'), \
        (3,'rec_Ep_Em'), \
        (4,'rec_Epos'), \
        (5,'rec_Mll')]
cb.AddObservations(['*'],['tt'],['13TeV'],['mtop'],cats)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['DY'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['STs'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['STt'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['TTZ'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['TTW'],  cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['WJets'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['WW'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['WZ'],   cats,False)
cb.AddProcesses(['*'],   ['tt'],['13TeV'],['mtop'],['ZZ'],   cats,False)
cb.AddProcesses(masses,  ['tt'],['13TeV'],['mtop'],['tt'],cats,True)
cb.AddProcesses(masses,  ['tt'],['13TeV'],['mtop'],['tW'],cats,True)


# Add systematics later
#cb.cp().signals().AddSyst(cb,"width","emu",ch.SystMap('mass')(['4','5','6','7','8'], 1.))
#cb.cp().process(['bkg']).AddSyst(cb,"bkg_norm","rateParam",ch.SystMap('era')(['13TeV'], 1.))

#cb.cp().backgrounds().AddSyst(cb,"pileup","shape",ch.SystMap('era')(['13TeV'], 1.))

for sys,systype in systematics.items():
    cb.cp().signals().AddSyst(cb,sys,systype,ch.SystMap('era')(['13TeV'], 1.))

cb.cp().backgrounds().ExtractShapes("mtTemplatesForCH.root","$BIN/$PROCESS","$BIN/$PROCESS_$SYSTEMATIC")
cb.cp().signals().ExtractShapes("mtTemplatesForCH.root","$BIN/$PROCESS$MASS","$BIN/$PROCESS$MASS_$SYSTEMATIC")


#cb.cp().backgrounds().ExtractShapes("mtTemplatesForCH.root","$BIN/$PROCESS","")
#cb.cp().signals().ExtractShapes("mtTemplatesForCH.root","$BIN/$PROCESS$MASS","")

#bbb = ch.BinByBinFactory()
#bbb.SetAddThreshold(0.1).SetMergeThreshold(0.5).SetFixNorm(True)
#bbb.MergeBinErrors(cb.cp().backgrounds())
#bbb.AddBinByBin(cb.cp().backgrounds(), cb)
#bbb.MergeBinErrors(cb.cp().signals())
#bbb.AddBinByBin(cb.cp().signals(), cb)


w = ROOT.RooWorkspace('morph', 'morph')

mT = ROOT.RooRealVar("MT","",1730,1650,1800)
mT.setConstant(True)


debug = ROOT.TFile('debug.root', 'RECREATE')

for b in cb.cp().channel(['mtop']).bin_set():   # rec_ptll, rec_ptpos
    for sig in signals: # tt, tW
        morphing.BuildRooMorphing(w, cb, b, sig, mT,force_template_limit=False,file=debug, interp_mode="CSPLINE")

cb.AddWorkspace(w,True)

cb.cp().signals().ExtractPdfs(cb, 'morph', '$BIN_$PROCESS_morph', '')
# cb.cp().process(['signal']).ExtractPdfs(cb, 'morph', '$BIN/$PROCESS_morph', '')


cb.PrintAll()
cb.WriteDatacard("cardMorph.txt", "outputfileMorph.root")

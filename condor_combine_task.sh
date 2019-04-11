#!/bin/sh
ulimit -s unlimited
set -e

export USER=$(whoami)
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    echo "Running In Batch"
    (>&2 echo "Starting job on " `date`) # Date/time of start of job
    (>&2 echo "Running on: `uname -a`") # Condor job is running on this node
    (>&2 echo "System software: `cat /etc/redhat-release`") # Operating System on that node

    cd ${_CONDOR_SCRATCH_DIR}
    echo ${_CONDOR_SCRATCH_DIR}

    echo "Current directory"
    pwd

# copy tarred cmssw area over from eos (should be excluding .SCRAM area)
    echo "xrdcp root://cmseos.fnal.gov//store/user/"${USER}"/condorFiles/combine_CMSSW_8_1_0.tgz ."
    xrdcp root://cmseos.fnal.gov//store/user/${USER}/condorFiles/combine_CMSSW_8_1_0.tgz .

    export SCRAM_ARCH=slc6_amd64_gcc530
    source /cvmfs/cms.cern.ch/cmsset_default.sh

    eval `scramv1 project CMSSW CMSSW_8_1_0`

    echo "tar -xvf combine_CMSSW_8_1_0.tgz"
    tar -xzf combine_CMSSW_8_1_0.tgz 
    rm combine_CMSSW_8_1_0.tgz
    cd CMSSW_8_1_0/src/UserCode/mtMorphing
    eval `scramv1 runtime -sh`

    ls -alhtr
fi

if [ $1 -eq 0 ]; then
  combine -M MultiDimFit -n _paramFit_Test_BkgNorm --algo impact --redefineSignalPOIs r,MT -P BkgNorm --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 1 ]; then
  combine -M MultiDimFit -n _paramFit_Test_Lumi --algo impact --redefineSignalPOIs r,MT -P Lumi --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 2 ]; then
  combine -M MultiDimFit -n _paramFit_Test_BTagSF --algo impact --redefineSignalPOIs r,MT -P BTagSF --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 3 ]; then
  combine -M MultiDimFit -n _paramFit_Test_CRGluon --algo impact --redefineSignalPOIs r,MT -P CRGluon --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 4 ]; then
  combine -M MultiDimFit -n _paramFit_Test_CRQCD --algo impact --redefineSignalPOIs r,MT -P CRQCD --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 5 ]; then
  combine -M MultiDimFit -n _paramFit_Test_CRerdON --algo impact --redefineSignalPOIs r,MT -P CRerdON --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 6 ]; then
  combine -M MultiDimFit -n _paramFit_Test_EleIDEff --algo impact --redefineSignalPOIs r,MT -P EleIDEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 7 ]; then
  combine -M MultiDimFit -n _paramFit_Test_EleRecoEff --algo impact --redefineSignalPOIs r,MT -P EleRecoEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 8 ]; then
  combine -M MultiDimFit -n _paramFit_Test_EleScale --algo impact --redefineSignalPOIs r,MT -P EleScale --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 9 ]; then
  combine -M MultiDimFit -n _paramFit_Test_EleSmear --algo impact --redefineSignalPOIs r,MT -P EleSmear --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 10 ]; then
  combine -M MultiDimFit -n _paramFit_Test_JEC --algo impact --redefineSignalPOIs r,MT -P JEC --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 11 ]; then
  combine -M MultiDimFit -n _paramFit_Test_JER --algo impact --redefineSignalPOIs r,MT -P JER --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 12 ]; then
  combine -M MultiDimFit -n _paramFit_Test_MuIDEff --algo impact --redefineSignalPOIs r,MT -P MuIDEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 13 ]; then
  combine -M MultiDimFit -n _paramFit_Test_MuIsoEff --algo impact --redefineSignalPOIs r,MT -P MuIsoEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 14 ]; then
  combine -M MultiDimFit -n _paramFit_Test_MuScale --algo impact --redefineSignalPOIs r,MT -P MuScale --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 15 ]; then
  combine -M MultiDimFit -n _paramFit_Test_MuTrackEff --algo impact --redefineSignalPOIs r,MT -P MuTrackEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 16 ]; then
  combine -M MultiDimFit -n _paramFit_Test_Pdf --algo impact --redefineSignalPOIs r,MT -P Pdf --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 17 ]; then
  combine -M MultiDimFit -n _paramFit_Test_Q2 --algo impact --redefineSignalPOIs r,MT -P Q2 --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 18 ]; then
  combine -M MultiDimFit -n _paramFit_Test_TrigEff --algo impact --redefineSignalPOIs r,MT -P TrigEff --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 19 ]; then
  combine -M MultiDimFit -n _paramFit_Test_UE --algo impact --redefineSignalPOIs r,MT -P UE --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 20 ]; then
  combine -M MultiDimFit -n _paramFit_Test_fsr --algo impact --redefineSignalPOIs r,MT -P fsr --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 21 ]; then
  combine -M MultiDimFit -n _paramFit_Test_hdamp --algo impact --redefineSignalPOIs r,MT -P hdamp --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 22 ]; then
  combine -M MultiDimFit -n _paramFit_Test_isr --algo impact --redefineSignalPOIs r,MT -P isr --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 23 ]; then
  combine -M MultiDimFit -n _paramFit_Test_pileup --algo impact --redefineSignalPOIs r,MT -P pileup --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi
if [ $1 -eq 24 ]; then
  combine -M MultiDimFit -n _paramFit_Test_toppt --algo impact --redefineSignalPOIs r,MT -P toppt --floatOtherPOIs 1 --saveInactivePOI 1 --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --setParameterRanges MT=1665,1785 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t 10 --saveToys -m 125 -d cardMorph.root
fi

ls -alhtr

if [ ! -z ${_CONDOR_SCRATCH_DIR} ] ; then
    cp robustHesse*.root higgsCombine*.root  ${_CONDOR_SCRATCH_DIR}/
    ls -alhtr ${_CONDOR_SCRATCH_DIR}
fi

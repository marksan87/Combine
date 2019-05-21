#!/bin/bash
ulimit -s unlimited
set -e

export USER=$(whoami)
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
    mkdir -p interactive
    pushd interactive
else
    echo "Running In Batch"
    echo "ulimit -n  " `ulimit -n`
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

    #eval `scramv1 project CMSSW CMSSW_8_1_0`
    scramv1 project CMSSW CMSSW_8_1_0
    echo "tar -xvf combine_CMSSW_8_1_0.tgz"
    tar -xzf combine_CMSSW_8_1_0.tgz 
    rm combine_CMSSW_8_1_0.tgz
    cd CMSSW_8_1_0/src/UserCode/mtMorphing
    eval `scramv1 runtime -sh`

    ls -alhtr
fi

if [ ! -z "$2" ] ; then
  toys="-t $2"
else
  toys="-t -1"
fi

parameterArgs="--setParameters MT=1725 --setParameterRanges MT=1665,1785:r=0,2 -S 0"

fitArgs="--redefineSignalPOIs MT,r --floatOtherPOIs 1 --saveInactivePOI 1 --expectSignal 1 --robustFit 1"
#fitArgs="--redefineSignalPOIs MT,r --saveInactivePOI 1 --expectSignal 1 --robustFit 1 --robustHesse 1"



#fitArgs="--redefineSignalPOIs MT,r --floatOtherPOIs 1 --saveInactivePOI 1 --expectSignal 1 --robustFit 1 --robustHesse 1 --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1"
#cardroot="diff_mt1725_cardMorph.root"
card="morphRates_ptll_varbins_cardMorph"
cardroot="${card}.root"

echo "text2workspace.py ${card}.txt"
text2workspace.py ${card}.txt

name=""

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    cardroot="../../$cardroot"
fi

if [ $1 -eq 0 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_BkgNorm --algo impact -P BkgNorm $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_BkgNorm --algo impact -P BkgNorm $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_BkgNorm"
fi
if [ $1 -eq 1 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_Lumi --algo impact -P Lumi $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_Lumi --algo impact -P Lumi $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_Lumi"
fi
if [ $1 -eq 2 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_BTagSF --algo impact -P BTagSF $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_BTagSF --algo impact -P BTagSF $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_BTagSF"
fi
if [ $1 -eq 3 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_CRGluon --algo impact -P CRGluon $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_CRGluon --algo impact -P CRGluon $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_CRGluon"
fi
if [ $1 -eq 4 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_CRQCD --algo impact -P CRQCD $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_CRQCD --algo impact -P CRQCD $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_CRQCD"
fi
if [ $1 -eq 5 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_CRerdON --algo impact -P CRerdON $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_CRerdON --algo impact -P CRerdON $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_CRerdON"
fi
if [ $1 -eq 6 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_EleIDEff --algo impact -P EleIDEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_EleIDEff --algo impact -P EleIDEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_EleIDEff"
fi
if [ $1 -eq 7 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_EleRecoEff --algo impact -P EleRecoEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_EleRecoEff --algo impact -P EleRecoEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_EleRecoEff"
fi
if [ $1 -eq 8 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_EleScale --algo impact -P EleScale $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_EleScale --algo impact -P EleScale $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_EleScale"
fi
if [ $1 -eq 9 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_EleSmear --algo impact -P EleSmear $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_EleSmear --algo impact -P EleSmear $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_EleSmear"
fi
if [ $1 -eq 10 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_JEC --algo impact -P JEC $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_JEC --algo impact -P JEC $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_JEC"
fi
if [ $1 -eq 11 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_JER --algo impact -P JER $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_JER --algo impact -P JER $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_JER"
fi
if [ $1 -eq 12 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_MuIDEff --algo impact -P MuIDEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_MuIDEff --algo impact -P MuIDEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_MuIDEff"
fi
if [ $1 -eq 13 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_MuIsoEff --algo impact -P MuIsoEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_MuIsoEff --algo impact -P MuIsoEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_MuIsoEff"
fi
if [ $1 -eq 14 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_MuScale --algo impact -P MuScale $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_MuScale --algo impact -P MuScale $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_MuScale"
fi
if [ $1 -eq 15 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_MuTrackEff --algo impact -P MuTrackEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_MuTrackEff --algo impact -P MuTrackEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_MuTrackEff"
fi
if [ $1 -eq 16 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_Pdf --algo impact -P Pdf $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_Pdf --algo impact -P Pdf $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_Pdf"
fi
if [ $1 -eq 17 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_Q2 --algo impact -P Q2 $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_Q2 --algo impact -P Q2 $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_Q2"
fi
if [ $1 -eq 18 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_TrigEff --algo impact -P TrigEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_TrigEff --algo impact -P TrigEff $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_TrigEff"
fi
if [ $1 -eq 19 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_UE --algo impact -P UE $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_UE --algo impact -P UE $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_UE"
fi
if [ $1 -eq 20 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_fsr --algo impact -P fsr $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_fsr --algo impact -P fsr $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_fsr"
fi
if [ $1 -eq 21 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_hdamp --algo impact -P hdamp $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_hdamp --algo impact -P hdamp $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_hdamp"
fi
if [ $1 -eq 22 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_isr --algo impact -P isr $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_isr --algo impact -P isr $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_isr"
fi
if [ $1 -eq 23 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_pileup --algo impact -P pileup $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_pileup --algo impact -P pileup $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_pileup"
fi
if [ $1 -eq 24 ]; then
  echo "combine -M MultiDimFit -n _paramFit_Test_toppt --algo impact -P toppt $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _paramFit_Test_toppt --algo impact -P toppt $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_paramFit_Test_toppt"
fi
if [ $1 -eq 25 ]; then
  echo "combine -M MultiDimFit -n _initialFit_Test --algo singles $fitArgs $parameterArgs $toys -m 1725 -d $cardroot"
  combine -M MultiDimFit -n _initialFit_Test --algo singles $fitArgs $parameterArgs $toys -m 1725 -d $cardroot
  name="_initialFit_Test"
fi


if [ ! -z ${_CONDOR_SCRATCH_DIR} ] ; then
    ls -alhtr
    # Only copy hesse root file if it exists
    [ -f robustHesse${name}.root ] && cp robustHesse${name}.root ${_CONDOR_SCRATCH_DIR}/
    
    cp higgsCombine${name}.MultiDimFit.mH1725*.root  ${_CONDOR_SCRATCH_DIR}/
    rm -rf ${_CONDOR_SCRATCH_DIR}/docker_stderror
    ls -alhtr ${_CONDOR_SCRATCH_DIR}
else
    popd
fi

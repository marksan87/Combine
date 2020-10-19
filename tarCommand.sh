pushd /uscms_data/d3/msaunder/combine 
echo "Current directory: "
pwd
tar --exclude='*~' --exclude='*.json' --exclude='*.pyc' --exclude='*.gif' --exclude='*.pdf' --exclude='*.o' --exclude='CMSSW_10_2_13/src/UserCode/mtMorphing/RootFiles' --exclude='*_bias.root' --exclude='robustHesse_param*.root' --exclude='higgsCombine*.root' --exclude='CMSSW_10_2_13/src/UserCode/mtMorphing/toys' --exclude-vcs  --exclude='debug.root' --exclude='plots.root' --exclude='CMSSW_10_2_13/.SCRAM' --exclude='CMSSW_10_2_13/src/UserCode/old' --exclude 'CMSSW_10_2_13/src/UserCode/mtMorphing/.git' --exclude 'CMSSW_10_2_13/src/UserCode/mtMorphing/plots2018' --exclude 'CMSSW_10_2_13/src/UserCode/mtMorphing/pol*moments' --exclude 'CMSSW_10_2_13/src/UserCode/mtMorphing/bins' --exclude='*_mtTemplatesForCH*.root' --exclude='CMSSW_10_2_13/src/auxiliaries' --exclude='*.git*' --exclude='*tgz' --exclude='*pkl' --exclude '*pklz' --exclude='*.png' --exclude='*.out' --exclude='*.err' --exclude='*.log' --exclude='*moments' --exclude='CMSSW_10_2_13/src/CombineHarvester/*.root' --exclude='CMSSW_10_2_13/src/HiggsAnalysis/*.root' -zcf combine_CMSSW_10_2_13.tgz CMSSW_10_2_13

xrdcp -f combine_CMSSW_10_2_13.tgz root://cmseos.fnal.gov//store/user/${USER}/condorFiles/.

rm combine_CMSSW_10_2_13.tgz 
popd


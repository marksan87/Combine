#!/bin/bash
outputName="impacts"

minMT=1665
maxMT=1785

#toyArg=""
nosyst=0
if [ "$1" == "-i" ] ; then
    card="$2"
    if [ "$3" == "-o" ] ; then
        outputName="$4"
        if [ "$5" == "--nosyst" ] ; then
            nosyst=1
        fi
    else
        outputName="scan"
    fi
elif [ -z "$1" ] ; then
    # Default to cardMorph
    echo "Defaulting to cardMorph"
    card="cardMorph"
fi

precision=$(( ${#minMT}-3 ))    # Number of decimal places
if [ $precision -le 0  ] ; then
    mtScaled="MT"
else
    mtScaled="MT/$(( 10**precision ))"
fi

initialMT=$(( minMT + (maxMT - minMT)/2 ))
#parameterArgs="--redefineSignalPOIs MT,r --floatOtherPOIs 1 --expectSignal 1 --setParameters MT=${initialMT} --setParameterRanges MT=${minMT},${maxMT}:r=0.75,1.25"
#parameterArgs="--redefineSignalPOIs MT,r --floatOtherPOIs 1 --expectSignal 1 --setParameters MT=${initialMT} --setParameterRanges MT=${minMT},${maxMT}:r=0.75,1.25"
parameterArgs="--redefineSignalPOIs MT,r --floatOtherPOIs 1 --expectSignal 1 --setParameters MT=${initialMT} --setParameterRanges MT=${minMT},${maxMT} --rMin 0 --rMax 2"

#seed=1
#toyArg="-t -1 --seed $seed"
toyArg="-t -1"

cardRoot="../${card}.root"

echo "text2workspace.py ${card}.txt"
text2workspace.py ${card}.txt

echo "mkdir -p ${outputName}"
mkdir -p ${outputName}

echo "pushd ${outputName}"
pushd ${outputName}

# Full uncertainty
echo "combine -M MultiDimFit -m 125 ${cardRoot} -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n nominal"
combine -M MultiDimFit -m 125 ${cardRoot} -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n nominal

# Save workspace
echo "combine -M MultiDimFit ${cardRoot} -m 125 -P MT --robustFit 1 --algo none $parameterArgs ${toyArg} -n bestfit --saveWorkspace"
combine -M MultiDimFit ${cardRoot} -m 125 -P MT --robustFit 1 --algo none $parameterArgs ${toyArg} -n bestfit --saveWorkspace

if [ $nosyst == 1 ] ; then 
    echo "./plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --output ${outputName}"
    ../plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --output ${outputName}
else

    # Freeze all systematics
    echo "combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n stat --snapshotName MultiDimFit --freezeNuisanceGroups all"
    combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n stat --snapshotName MultiDimFit --freezeNuisanceGroups all

    # Freeze theory
    echo "combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n theory --snapshotName MultiDimFit --freezeNuisanceGroups theory"
    combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n theory --snapshotName MultiDimFit --freezeNuisanceGroups theory 


    # Freeze exp
    echo "combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n exp --snapshotName MultiDimFit --freezeNuisanceGroups exp"
    combine -M MultiDimFit -m 125 higgsCombinebestfit.MultiDimFit.mH125.root -P MT --robustFit 1 --algo grid --points 100 $parameterArgs ${toyArg} -n exp --snapshotName MultiDimFit --freezeNuisanceGroups exp 


    #echo "./plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --output ${outputName}"
    #./plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --output ${outputName}

    echo "../plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --others 'higgsCombinestat.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat --output ${outputName}"
    ../plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --POI ${mtScaled} --others 'higgsCombinestat.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat --output ${outputName}

    #plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --others 'higgsCombinetheory.MultiDimFit.mH125.root:Freeze theory:4' 'higgsCombineexp.MultiDimFit.mH125.root:Freeze exp:3' 'higgsCombinestat.MultiDimFit.mH125.root:Freeze systs:2' --breakdown theory,exp,syst,stat --POI MT/10
fi

popd

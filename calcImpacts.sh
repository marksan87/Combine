#!/bin/bash
rebin=$1
cut=$2
obs=$3
reco=$4
outputName="impacts"

if [ "$1" == "-i" ] ; then
    card="$2"
    if [ "$3" == "-o" ] ; then
        outputName="$4"
    else
        outputName="$2_impacts"
    fi
elif [ -z "$1" ] ; then
    # Default to cardMorph
    echo "Defaulting to cardMorph"
    card="cardMorph"
else
    if [ -z "$2" ] ; then
        card="rebin_${rebin}_cardMorph"
        outputName="rebin_${rebin}_impacts"
    else
        if [ -z "$3" ] ; then
            echo "Defaulting to rec_ptll"
            obs="ptll"
            reco="rec"
        elif [ -z "$4" ] ; then
            echo "Defaulting to rec"
            reco="rec"
        fi
        card="${reco}_${obs}_rebin_${rebin}_cut_${cut}_cardMorph"
        outputName="${reco}_${obs}_rebin_${rebin}_cut_${cut}_impacts"
    fi
fi

minMT=1695
maxMT=1755

cardRoot="${card}.root"

echo "text2workspace.py ${card}.txt"
text2workspace.py ${card}.txt


echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doInitialFit --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1"
combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doInitialFit --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1 

echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doFits --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1 -v 2 2>&1 | tee log_${card}.log"
combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doFits --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1 -v 2 2>&1 | tee log_${card}.log


#echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doFits --parallel 8 --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1"
#combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 --robustFit 1 --robustHesse 1 --doFits --parallel 8 --setParameterRanges MT=${minMT},${maxMT} --maxFailedSteps 9999999 --stepSize 0.005 --setRobustFitTolerance 0.1 -t -1


echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 -o ${outputName}.json"
combineTool.py -M Impacts -d ${cardRoot} -m 125 --redefineSignalPOIs r,MT --setParameters MT=1725 --expectSignal 1 -o ${outputName}.json

echo "./createImpactPlot.py -i ${outputName}.json -o ${outputName}"
./createImpactPlot.py -i ${outputName}.json -o ${outputName}


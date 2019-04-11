#!/bin/bash
rebin=$1
cut=$2
obs=$3
reco=$4
outputName="impacts"

minMT=1665
maxMT=1785

if [ "$1" == "-i" ] ; then
    card="$2"
    if [ "$3" == "-o" ] ; then
        outputName="$4"
        if [ ! -z "$5" ] ; then
            if [ "$5" == "--nosyst" ] ; then
                nosyst=1
            else
                minMT=$5
                if [ ! -z "$6" ] ; then
                    maxMT=$6
                fi
            fi
        fi
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

initialMT=$(( minMT + (maxMT - minMT)/2 ))
parameterArgs="--redefineSignalPOIs MT,r --setParameters MT=${initialMT} --setParameterRanges MT=${minMT},${maxMT}:r=0,2"
#poiArgs="--redefineSignalPOIs r,MT --floatOtherPOIs 1 --saveInactivePOI 1 --expectSignal 1"
#fitArgs="--robustFit 1 --robustHesse 1"

fitArgs="--robustFit 1"
#toyArgs="-t -1 --seed 1"
toyArgs="-t -1"
#toyArgs=""

#parallel="--parallel 8"
parallel=""

#debug="-v 2 2>&1 | tee log_${card}.log"
debug=""



if [ "$nosyst" == "1" ] ; then
    echo "Ignoring systematics"
    fitArgs="$fitArgs -S 0"
fi


cardRoot="../${card}.root"

echo "text2workspace.py ${card}.txt"
text2workspace.py ${card}.txt

echo "mkdir -p ${outputName}"
mkdir -p ${outputName}

echo "pushd ${outputName}"
pushd ${outputName}


echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --expectSignal 1 --floatOtherPOIs 1 $fitArgs --doInitialFit $parameterArgs $toyArgs"
combineTool.py -M Impacts -d ${cardRoot} -m 125 --expectSignal 1 --floatOtherPOIs 1 $fitArgs --doInitialFit $parameterArgs $toyArgs 

if [ -z "$nosyst" ] ; then
    echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 --expectSignal 1 $fitArgs --doFits $parallel $parameterArgs $toyArgs $debug"
    combineTool.py -M Impacts -d ${cardRoot} -m 125 --expectSignal 1 $fitArgs --doFits $parallel $parameterArgs $toyArgs $debug 

    echo "combineTool.py -M Impacts -d ${cardRoot} -m 125 $parameterArgs --expectSignal 1 -o ${outputName}.json"
    combineTool.py -M Impacts -d ${cardRoot} -m 125 $parameterArgs --expectSignal 1 -o ${outputName}.json

    echo "../createImpactPlot.py -i ${outputName}.json -o ${outputName}"
    ../createImpactPlot.py -i ${outputName}.json -o ${outputName}
fi

popd


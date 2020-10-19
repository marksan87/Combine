#!/bin/bash

inDir=$1
outDir=/uscms/homes/m/msaunder/private/documentation/AN-17-249/assets
templateDir=$outDir/templates
momentDir=$outDir/moments


if [ "$inDir" == "" ]; then
    echo "Usage: ./copyToAN.sh plotDir"
    exit 1
fi

if [ "${inDir:${#inDir}-1:1}" == "/" ] ; then
    inDir=${inDir:0:${#inDir}-1}
fi

declare -a obs=(
    "rec_ptll" \
    "rec_Mll" \
    "rec_ptpos" \
    "rec_Epos" \
    "rec_ptp_ptm" \
    "rec_Ep_Em" \
)

declare -a systematics=(
    "herwigpp" \
    "madgraph" \
    "amcanlo" \
    "CRQCD" \
    "CRGluon" \
    "CRerdON" \
    "UE" \
    "hdamp" \
    "fsr" \
    "isr" \
    "Pdf" \
    "Q2" \
    "toppt" \
    "JER" \
    "JEC" \
    "BTagSF" \
    "TrigEff" \
    "MuScale" \
    "MuTrackEff" \
    "MuIsoEff" \
    "MuIDEff" \
    "EleSmear" \
    "EleScale" \
    "EleRecoEff" \
    "EleIDEff" \
    "Lumi" \
    "pileup" \
    "DS" \
)

declare -a ttOnlySysts=(
        "Pdf" \
        "toppt" \
        "UE" \
        "CRerdON" \ 
        "CRGluon" \
        "CRQCD" \
        "amcanlo" \
        "madgraph" \
        "herwigpp" \
        )


declare -a tWOnlySysts=(
        "DS" \
        )


declare -a signal=(
    "TTbar" \
    "ST_tW" \
)

declare -a plots=(
    "histplots/mtscan/mtscan_TTbar_actual_nominal" \
    "histplots/mtscan/mtscan_TTbar_morph_nominal" \
    "histplots/mtscan/mtscan_ST_tW_actual_nominal" \
    "histplots/mtscan/mtscan_ST_tW_morph_nominal" \
)

declare -a systPlots=(
    "histplots/TTbar_actual_hist" \
    "histplots/TTbar_morph_hist" \
    "histplots/ST_tW_actual_hist" \
    "histplots/ST_tW_morph_hist" \
)


echo "Copying templates from $inDir to $templateDir"

for (( o=0; o < ${#obs[@]}; o++ )); do
    mkdir -p $templateDir/${obs[$o]}
    mkdir -p $momentDir/${obs[$o]}
    
    # Mass scans
    for (( pl=0; pl < ${#plots[@]}; pl++ )); do
        cp $inDir/${obs[$o]}/${plots[$pl]}_${obs[$o]}.pdf $templateDir/${obs[$o]}/
    done

    # Systematics
    for (( sys=0; sys < ${#systematics[@]}; sys++ )); do
        # Copy moments
        cp $inDir/${obs[$o]}/*${obs[$o]}_${systematics[$sys]}_m*.pdf $momentDir/${obs[$o]}/
        
        for (( s=0; s < ${#signal[@]}; s++ )); do
            if [ "${signal[$s]}" == "TTbar" ] && [[ "${tWOnlySysts[@]}" =~ "${systematics[$sys]}" ]] ; then
                continue
            elif [ "${signal[$s]}" == "ST_tW" ] && [[ "${ttOnlySysts[@]}" =~ "${systematics[$sys]}" ]] ; then
                continue
            fi
                
            # Copy templates
            cp $inDir/${obs[$o]}/histplots/${signal[$s]}_*_hist_${obs[$o]}_${systematics[$sys]}.pdf $templateDir/${obs[$o]}/

        done
    done

done


echo "Copying moments from $inDir to $momentDir"


#!/bin/bash

inDir=$1
outDir=/uscms/homes/m/msaunder/private/documentation/AN-17-249/assets

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
        "Q2" \
        "Pdf" \
        "toppt" \
        "hdamp" \
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
    "histplots/mtscan/mtscan_TTbar_actual" \
    "histplots/mtscan/mtscan_TTbar_morph" \
    "histplots/mtscan/mtscan_ST_tW_actual" \
    "histplots/mtscan/mtscan_ST_tW_morph" \
)

declare -a systPlots=(
    "histplots/TTbar_actual_hist" \
    "histplots/TTbar_morph_hist" \
    "histplots/ST_tW_actual_hist" \
    "histplots/ST_tW_morph_hist" \
)

mkdir -p $outDir
echo "Copying plots from $inDir to $outDir"

for (( o=0; o < ${#obs[@]}; o++ )); do
    mkdir -p $outDir/${obs[$o]}
    
    # Mass scans
    for (( pl=0; pl < ${#plots[@]}; pl++ )); do
        cp $inDir/${obs[$o]}/${plots[$pl]}_${obs[$o]}.pdf $outDir/${obs[$o]}/
    done

    # Systematics
    for (( s=0; s < ${#signal[@]}; s++ )); do
        for (( sys=0; sys < ${#systematics[@]}; sys++ )); do
            if [ "${signal[$s]}" == "TTbar" ] && [[ "${tWOnlySysts[@]}" =~ "${systematics[$sys]}" ]] ; then
                continue
            elif [ "${signal[$s]}" == "ST_tW" ] && [[ "${ttOnlySysts[@]}" =~ "${systematics[$sys]}" ]] ; then
                continue
            fi

            cp $inDir/${obs[$o]}/histplots/${signal[$s]}_actual_hist_${obs[$o]}_${systematics[$sys]}.pdf $outDir/${obs[$o]}/
        done
    done


#    for (( p=0; p < ${#systPlots[@]}; p++ )); do
#        for (( sys=0; sys < ${#systematics[@]}; sys++ )); do
#            cp $inDir/${obs[$o]}/${systPlots[$p]}_${obs[$o]}_${systematics[$sys]}.pdf $outDir/${obs[$o]}/
#        done
#    done

done



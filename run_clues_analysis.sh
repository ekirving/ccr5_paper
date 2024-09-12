#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# CCR5 analysis - 2024-07-22

# fetch the CLUES code
[ ! -d "bin/clues/" ] && git clone git@github.com:35ajstern/clues.git --branch aaron/palm-integration bin/clues/

# CCR5 and the 7 controls
DELETIONS=( rs333 rs61231801 rs66552573 rs67580019 rs143241023 rs150628438 rs369842709 rs556322139 )

# the different models to run
MODELS=( No_filter Permissive_filter Strict_filter p_Data )

# the different ascertainment modes
MODES=( no_freq mod_freq )

# make the output folders
parallel "mkdir -p clues/{}/" ::: "${DELETIONS[@]}"

# make the CLUES inputs
Rscript scripts/make_clues_inputs.R

# make the single epoch time bins
printf "0\n529\n" > clues/one-epoch.bins

for rsid in "${DELETIONS[@]}"; do

    # frequency of the deletion in FIN, GBR and TSI
    if [ "${rsid}" == "rs333" ]; then
        freq=0.1103
    elif [ "${rsid}" == "rs61231801" ]; then
        freq=0.1034
    elif [ "${rsid}" == "rs66552573" ]; then
        freq=0.1223
    elif [ "${rsid}" == "rs67580019" ]; then
        freq=0.1054
    elif [ "${rsid}" == "rs143241023" ]; then
        freq=0.0944
    elif [ "${rsid}" == "rs150628438" ]; then
        freq=0.0974
    elif [ "${rsid}" == "rs369842709" ]; then
        freq=0.1203
    elif [ "${rsid}" == "rs556322139" ]; then
        freq=0.1004
    else
        freq=-1
    fi

    for model in "${MODELS[@]}"; do

        for mode in "${MODES[@]}"; do

            echo "Running ${rsid} w/ ${model} and ${mode}"

            if [ "${mode}" == "mod_freq" ]; then
                daf="--popFreq ${freq}"
            else
                daf=""
            fi

            # run CLUES
            python bin/clues/inference.py \
                 --lik \
                 ${daf} \
                 --coal "relate/1000G_phase3-FIN_GBR_TSI-popsize.coal" \
                 --ancientSamps "clues/${rsid}/${rsid}-${model}.ancient" \
                 --timeBins "clues/one-epoch.bins" \
                 --betaParam 0.5 \
                 --out "clues/${rsid}/${rsid}-${model}-${mode}" \
                     &> "clues/${rsid}/${rsid}-${model}-${mode}.log"

            # extract the results
            python scripts/clues_parse_log.py \
              --rsid "${rsid}-${model}-${mode}" \
              --ancestry ALL \
              --use-freq "${mode}" \
              --mod-freq "${freq}" \
              --log "clues/${rsid}/${rsid}-${model}-${mode}.log" \
              --out "clues/${rsid}/${rsid}-${model}-${mode}.json"

            num_samples=$(wc -l < "clues/${rsid}/${rsid}-${model}.ancient" | sed 's/ *//')

            if [ "${mode}" == "mod_freq" ]; then
                mod_label="Modern DAF = ${freq}"
            else
                mod_label="No modern DAF"
            fi

            # make the label
            echo "{\"title\":\"${rsid} | ${model} | ${mod_label} | n=${num_samples}\", \"gwascat\": []}" \
                     > "clues/${rsid}/${rsid}-${model}-${mode}-label.json"

            # plot the trajectory
            python scripts/clues_plot_trajectory.py \
              --gen-time 28 \
              --params "clues/${rsid}/${rsid}-${model}-${mode}.json" \
              --label "clues/${rsid}/${rsid}-${model}-${mode}-label.json" \
              --ancestry ALL \
              --ext png \
              "clues/${rsid}/${rsid}-${model}-${mode}" \
              "clues/${rsid}/${rsid}-${model}-${mode}" 2> /dev/null
        done;
    done;
done;

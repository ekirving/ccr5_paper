#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# CCR5 analysis - 2024-07-22

# fetch the CLUES code
[ ! -d "bin/clues/" ] && git clone git@github.com:35ajstern/clues.git --branch aaron/palm-integration bin/clues/

# CCR5 and the 7 controls
DELETIONS=( rs333 rs61231801 rs66552573 rs67580019 rs143241023 rs150628438 rs369842709 rs556322139 Le_proxy_SNP )

# the different models to run
MODELS=( No_filter Permissive_filter Strict_filter p_Data )

# the different ascertainment modes
MODES=( no_freq mod_freq )

# make the output folders
parallel "mkdir -p clues/{}/" ::: "${DELETIONS[@]}"

# make the CLUES inputs
Rscript scripts/make_clues_inputs.R

for rsid in "${DELETIONS[@]}"; do

    # frequency of the each deletion 1kGP EUR
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
    elif [ "${rsid}" == "Le_proxy_SNP" ]; then
        freq=0.1362  # rs73833033
    fi

    for model in "${MODELS[@]}"; do

        if [ "${rsid}" == "Le_proxy_SNP" ] && { [ "${model}" == "Permissive_filter" ] || [ "${model}" == "Strict_filter" ]; }; then
            continue
        fi

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
                 --timeBins "data/one-epoch.bins" \
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

            # make the label
            echo "{\"title\":\"${rsid} | ${model} | n=${num_samples} samples\", \"gwascat\": []}" \
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

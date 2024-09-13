#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# CCR5 analysis - 2024-07-22

# fetch the CLUES code
[ ! -d "bin/clues/" ] && git clone git@github.com:35ajstern/clues.git --branch aaron/palm-integration bin/clues/

# the ancestral paths
ANCESTRIES=( ALL ANA CHG WHG EHG )

# the different models to run
MODELS=( HAPI_samples West_Eurasian_samples )

# the different ascertainment modes
MODES=( no_freq mod_freq )

# make the output folders
parallel "mkdir -p clues/ccr5_tags-{}/" ::: "${MODELS[@]}"

# make the CLUES inputs
Rscript scripts/make_ancestry_inputs.R

for ancestry in "${ANCESTRIES[@]}"; do

    # average frequency of the 4 tag SNPs (i.e., rs113341849, rs113010081, rs79815064, rs11574435) in each ancestry path
    if [ "${ancestry}" == "ALL" ]; then
        freq=0.1237
    elif [ "${ancestry}" == "ANA" ]; then
        freq=0.0
    elif [ "${ancestry}" == "CHG" ]; then
        freq=0.30355
    elif [ "${ancestry}" == "WHG" ]; then
        freq=0.0284
    elif [ "${ancestry}" == "EHG" ]; then
        freq=0.185575
    fi

    for model in "${MODELS[@]}"; do

        for mode in "${MODES[@]}"; do

            echo "Running CCR5 tags for ${ancestry} w/ ${model} and ${mode}"

            if [ "${mode}" == "mod_freq" ]; then
                daf="--popFreq ${freq}"
            else
                daf=""
            fi

            if [ "${ancestry}" == "ALL" ]; then
                (
                    set -x
                    # run CLUES in diploid mode
                    python bin/clues/inference.py \
                        --lik \
                        ${daf} \
                        --coal "relate/1000G_phase3-FIN_GBR_TSI-popsize.coal" \
                        --ancientSamps "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}.ancient" \
                        --timeBins "clues/one-epoch.bins" \
                        --betaParam 0.5 \
                        --out "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}" \
                        &> "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}.log"
                )
            else
                (
                    set -x
                    # run CLUES in haploid mode
                    python bin/clues/inference.py \
                        --lik \
                        ${daf} \
                        --coal "relate/1000G_phase3-FIN_GBR_TSI-popsize.coal" \
                        --ancientHaps "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}.ancient" \
                        --timeBins "clues/one-epoch.bins" \
                        --betaParam 0.5 \
                        --out "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}" \
                        &> "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}.log"
                )
            fi

            # extract the results
            python scripts/clues_parse_log.py \
              --rsid "${ancestry}-${model}-${mode}" \
              --ancestry ALL \
              --use-freq "${mode}" \
              --mod-freq "${freq}" \
              --log "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}.log" \
              --out "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}.json"

            num_samples=$(wc -l < "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}.ancient" | sed 's/ *//')

            # make the label
            echo "{\"title\":\"${ancestry} | ${model} | n=${num_samples} samples\", \"gwascat\": []}" \
                     > "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}-label.json"

            # plot the trajectory
            python scripts/clues_plot_trajectory.py \
              --gen-time 28 \
              --params "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}.json" \
              --label "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}-label.json" \
              --ancestry ALL \
              --ext png \
              "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}" \
              "clues/ccr5_tags-${model}/ccr5_tags-${model}-${ancestry}-${mode}" 2> /dev/null
        done;
    done;
done;

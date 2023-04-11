#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# CCR5 analysis - 4th version - 2 selection epochs (14.8 kBP - 3 kBP - 0) - 2022-03-30

# fetch the CLUES code
[ ! -d "bin/clues/" ] && git clone https://github.com/standard-aaron/clues.git bin/clues/

# the different models to run
MODELS=( p_Data Predicted_by_model Artifacts_filter Minus_haplo_filter minus_haplotype )

# make the output folders
parallel "mkdir -p clues/{}/" ::: "${MODELS[@]}"

# make the CLUES input
Rscript scripts/ccr5_clues_inputs.R

# make the multi-epoch time bins
printf "0\n529\n" > clues/one-epoch.bins
printf "0\n107\n529\n" > clues/two-epochs.bins

EPOCHS=( one-epoch two-epochs )

# run all the models with the modern ascertainment
for model in "${MODELS[@]}"; do

  for epoch in "${EPOCHS[@]}"; do

	echo "Running ${model} for ${epoch}, with the modern ascertainment"

    if [ "${model}" == "minus_haplotype" ]
    then
       # frequency of the minus haplotype in FIN, GBR and TSI
       mod=0.0151
    else
       # frequency of the deletion in FIN, GBR and TSI
       mod=0.1237
    fi

    cd bin/clues/ && python inference.py \
      --popFreq ${mod} \
      --betaParam 0.5 \
      --coal ../../relate/1000G_phase3-FIN_GBR_TSI-popsize.coal \
      --ancientSamps ../../clues/${model}/ccr5-${model}.ancient \
      --timeBins ../../clues/${epoch}.bins \
      --out ../../clues/${model}/ccr5-${model}-${epoch}-mod \
      &> ../../clues/${model}/ccr5-${model}-${epoch}-mod.epochs.log
	
	cd ../..
	
    # extract the results
    python scripts/clues_parse_log.py \
      --rsid ccr5-${model}-mod \
      --mode ancient \
      --ancestry ALL \
      --sex any \
      --log clues/${model}/ccr5-${model}-${epoch}-mod.epochs.log \
      --out clues/${model}/ccr5-${model}-${epoch}-mod.json

    num_samples=$(cat clues/${model}/ccr5-${model}.ancient | wc -l)

    # make the label
    echo "{\"title\":\"CCR5 ${model} - Modern frequency ${mod} - ${epoch} (n=${num_samples})\", \"gwascat\": []}" \
      > clues/${model}/ccr5-${model}-${epoch}-mod-label.json

    # plot the trajectory
    python scripts/clues_plot_trajectory.py \
      --gen-time 28 \
      --params clues/${model}/ccr5-${model}-${epoch}-mod.json \
      --label  clues/${model}/ccr5-${model}-${epoch}-mod-label.json \
      --ancestry ALL \
      --sex any \
      --ext png \
      clues/${model}/ccr5-${model}-${epoch}-mod \
      clues/${model}/ccr5-${model}-${epoch}-mod 2> /dev/null

  done;
done;

# re-run the models without the modern ascertainment
for model in "${MODELS[@]}"; do

  for epoch in "${EPOCHS[@]}"; do

	echo "Running ${model} for ${epoch}, without the modern ascertainment"

    cd bin/clues/ && python inference.py \
      --betaParam 0.5 \
      --coal ../../relate/1000G_phase3-FIN_GBR_TSI-popsize.coal \
      --ancientSamps ../../clues/${model}/ccr5-${model}.ancient \
      --timeBins ../../clues/${epoch}.bins \
      --out ../../clues/${model}/ccr5-${model}-${epoch}-no_mod \
      &> ../../clues/${model}/ccr5-${model}-${epoch}-no_mod.epochs.log

	cd ../..

    # extract the results
    python scripts/clues_parse_log.py \
      --rsid ccr5-${model}-no_mod \
      --mode ancient \
      --ancestry ALL \
      --sex any \
      --log clues/${model}/ccr5-${model}-${epoch}-no_mod.epochs.log \
      --out clues/${model}/ccr5-${model}-${epoch}-no_mod.json

    num_samples=$(cat clues/${model}/ccr5-${model}.ancient | wc -l)

    # make the label
    echo "{\"title\":\"CCR5 ${model} - No modern ascertainment - ${epoch} (n=${num_samples})\", \"gwascat\": []}" \
      > clues/${model}/ccr5-${model}-${epoch}-no_mod-label.json

    # plot the trajectory
    python scripts/clues_plot_trajectory.py \
      --gen-time 28 \
      --params clues/${model}/ccr5-${model}-${epoch}-no_mod.json \
      --label  clues/${model}/ccr5-${model}-${epoch}-no_mod-label.json \
      --ancestry ALL \
      --sex any \
      --ext png \
      clues/${model}/ccr5-${model}-${epoch}-no_mod \
      clues/${model}/ccr5-${model}-${epoch}-no_mod 2> /dev/null

  done;
done;

# extract the age of the mutation
Rscript scripts/ccr5_mutation_age.R
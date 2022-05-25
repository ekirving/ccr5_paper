#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))
quiet(library(tibble))
quiet(library(stringr))
quiet(library(ggplot2))
quiet(library(RcppCNPy))

# the lowest discretisations of frequency are:
# 6.168503e-11
# 1.112985e-04
# 4.448131e-04
MIN_FREQ <- 1.112985e-04 # second lowest possible value!

clues_trajectory <- function(model_path) {

  # load the model data
  epochs <- npyLoad(paste0(model_path, ".epochs.npy"))
  freqs <- npyLoad(paste0(model_path, ".freqs.npy"))
  logpost <- npyLoad(paste0(model_path, ".post.npy"), dotranspose = F)

  # add column names
  colnames(logpost) <- paste0("V", seq(ncol(logpost)))

  model <- as_tibble(logpost) %>%

    # convert posterior densities from log-likelihoods
    exp() %>%

    # add the frequency labels
    add_column(freq = freqs, .before = 2) %>%

    # add the title heights (and a little padding to make sure there are no gaps)
    # tile heights are not equal because there is higher sampling density near 0 and 1
    add_column(height = diff(c(0, freqs)) + 1e-4, .before = 3) %>%

    # pivot the columns into long format
    pivot_longer(-c(freq, height), names_to = "epoch", values_to = "density") %>%

    # convert the column names into epochs (and switch the direction of time)
    mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))]) %>%

    # sort by epoch age
    arrange(epoch)

  model
}

models <- c("p_Data", "Predicted_by_model", "Artifacts_filter", "Minus_haplo_filter", "minus_haplotype")
epochs <- c("one-epoch", "two-epochs")
modes <- c("mod", "no_mod")

# load the trajectories
traj <- bind_rows(
  lapply(models, function(model) {
    bind_rows(
      lapply(epochs, function(epoch_type) {
        lapply(modes, function(mode) {
          clues_trajectory(paste0("clues/", model, "/ccr5-", model, "-", epoch_type, "-", mode)) %>%
            mutate(model = paste0(model, "-", epoch_type, "-", mode))
        })
      })
    )
  })
)

age <- traj %>%
  filter(freq <= MIN_FREQ) %>%
  group_by(model, epoch) %>%
  summarise(density = sum(density), freq = max(freq)) %>%
  filter(density >= 0.5) %>%
  slice(which.max(epoch)) %>%
  mutate(age_bp = -epoch * 28) %>%
  select(model, freq, epoch, age_bp, density)

write_tsv(age, "clues/ccr5_mutation_ages.tsv")

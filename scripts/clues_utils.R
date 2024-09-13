#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet(library(RcppCNPy))
quiet(library(zoo))

clues_trajectory <- function(model_path, smooth = 10) {

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

    # extract the maximum posterior trajectory
    traj <- model %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch)

    # apply a little smoothing to the jagged steps in the trajectory
    if (smooth) {
        traj$freq <- rollapply(c(traj$freq, rep_len(NA, smooth - 1)), width = smooth, by = 1, FUN = mean, na.rm = TRUE, align = "left")
    }

    traj
}

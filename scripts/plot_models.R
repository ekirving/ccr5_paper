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
quiet(library(directlabels))
quiet(library(zoo))
quiet(library(jsonlite))

# the average generation time in years
gen_time <- 28

# frequency of the deletion in FIN, GBR and TSI
modern_freq <- 0.1237

models <- c(
    # "Predicted_by_model", # not included in the main text
    "Artifacts_filter", # Permissive filter
    "Minus_haplo_filter" # Strict filter
)
epochs <- c("one-epoch")
modes <- c("mod", "no_mod")

# load the estimate allele ages
ages <- read_tsv("clues/ccr5_mutation_ages.tsv", col_types = cols()) %>%
    separate(col = model, into = c("model", "epochs", "_", "mode"), sep = "-") %>%
    filter(model %in% models & epochs == "one" & mode %in% modes)


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

# load the trajectories
traj <- bind_rows(
    lapply(models, function(model) {
        bind_rows(
            lapply(epochs, function(epoch_type) {
                lapply(modes, function(mode) {
                    prefix <- paste0("clues/", model, "/ccr5-", model, "-", epoch_type, "-", mode)
                    data <- fromJSON(paste0(prefix, ".json"))
                    pval <- formatC(pchisq(2 * data$logLR, df = 1, lower.tail = FALSE), format = "e", digits = 1)
                    clues_trajectory(prefix) %>%
                        mutate(
                            model = model,
                            name = case_when(
                                model == "Artifacts_filter" ~ "Permissive filter",
                                model == "Minus_haplo_filter" ~ "Strict filter",
                                TRUE ~ model
                            ),
                            label = paste0(name, " (p=", pval, ")"),
                            epoch_type = epoch_type,
                            mode = mode
                        )
                })
            })
        )
    })
)

# constrain the extent of the plotting
xmin <- min(traj$epoch)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
xlabels <- round(xbreaks * gen_time / 1000)

# set the factor order and labels
traj$mode <- factor(traj$mode, levels = c("mod", "no_mod"), labels = c("With modern ascertainment", "Without modern ascertainment"))
ages$mode <- factor(ages$mode, levels = c("mod", "no_mod"), labels = c("With modern ascertainment", "Without modern ascertainment"))

plt <- traj %>%
    # plot the trajectories
    ggplot(aes(x = epoch, y = freq, color = model)) +

    # split by mod and no_mod
    facet_grid(~mode) +

    # show the modern frequency
    # geom_point(x=0, y=modern_freq, color="red", shape=21, cex=3) +

    # plot the maximum posterior trajectory
    geom_line() +

    # show the allele ages
    geom_point(aes(y = 0), ages) +

    # set the model colours
    # scale_color_manual(values = snp_colors) +

    # print the labels
    geom_dl(aes(label = label), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.6), na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, .25), breaks = seq(0, 1, .05), position = "left") +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 450))) +

    # set the labels
    labs(title = "CCR5-delta 32") +
    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

ggsave("figure/ccr5_delta32_trajectory.png", plt, width = 9, height = 4)

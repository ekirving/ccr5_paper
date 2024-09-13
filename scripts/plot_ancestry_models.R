#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse))
quiet(library(directlabels))
quiet(library(jsonlite))

source("scripts/clues_utils.R")

# the average generation time in years
gen_time <- 28

# ancestry paths
ancestries <- c(
    "ALL",
    "ANA",
    "CHG",
    "WHG",
    "EHG"
)

models <- c(
    "HAPI_samples"
    # "West_Eurasian_samples"
)

modes <- c(
    "mod_freq",
    "no_freq"
)

# load the trajectories
traj <- bind_rows(
    lapply(ancestries, function(ancestry) {
        bind_rows(
            lapply(models, function(model) {
                lapply(modes, function(mode) {
                    prefix <- paste0("clues/ccr5_tags-", model, "/ccr5_tags-", model, "-", ancestry, "-", mode)
                    result <- fromJSON(paste0(prefix, ".json"))
                    clues_trajectory(prefix) %>%
                        mutate(model = model, ancestry = ancestry, mode = mode, pval = result$pval)
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
traj$mode <- factor(traj$mode, levels = c("mod_freq", "no_freq"), labels = c("With modern DAF", "Without modern DAF"))

traj$model <- factor(traj$model, levels = models, labels = str_replace_all(models, "_", " "))

# format the data for display
traj <- traj %>%
    mutate(
        significant = as.numeric(pval < (0.05 / length(ancestries))),
        # format the p-value
        pval = ifelse(
            pval < 0.05,
            sprintf("%.2e", signif(pval, 2)),
            sprintf("%.2f", signif(pval, 2))
        ),
        snp_label = paste0(ancestry, " (p=", pval, ")")
    )

ancestry_colors <- c(
    "ALL" = "#66c2a5",
    "WHG" = "#fc8d62",
    "EHG" = "#8da0cb",
    "CHG" = "#e78ac3",
    "ANA" = "#a6d854"
)

plt <- traj %>%
    # plot the trajectories
    ggplot(aes(x = epoch, y = freq, color = ancestry)) +

    # plot the maximum posterior trajectory
    geom_line(linewidth = 1, na.rm = TRUE) +

    # print the labels
    geom_dl(aes(label = snp_label), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.7), na.rm = TRUE) +

    # display as a grid
    facet_grid(model ~ mode, labeller = labeller(description = label_wrap_gen())) +

    # set the model colours
    scale_color_manual(values = ancestry_colors) +

    # plot non-significant trajectories as transparent
    scale_alpha(range = c(0.3, 1)) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, .32), breaks = seq(0, 1, 0.05), position = "left") +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 270))) +
    labs(
        title = "Ancestry stratified trajectories for CCR5Δ32"
    ) +
    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(0.5, "lines")
    )

ggsave("figure/ancestry_trajectories.png", plt, width = 9, height = 3)

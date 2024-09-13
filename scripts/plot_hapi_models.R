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

source("scrips/clues_utils.R")

# the average generation time in years
gen_time <- 28

# CCR5Δ32 and the paired controls
deletions <- c(
    "rs333", # CCR5Δ32
    "rs61231801",
    "rs66552573",
    "rs67580019",
    "rs143241023",
    "rs150628438",
    "rs369842709",
    "rs556322139"
)

models <- c(
    "Strict_filter",
    "Permissive_filter"
    # "No_filter",
    # "p_Data"
)

modes <- c(
    "mod_freq",
    "no_freq"
)

# load the estimated allele ages
ages <- read_tsv("clues/allele_ages.tsv", show_col_types = FALSE) %>%
    filter(rsid %in% deletions & model %in% models & mode %in% modes)

# load the trajectories
traj <- bind_rows(
    lapply(deletions, function(rsid) {
        bind_rows(
            lapply(models, function(model) {
                lapply(modes, function(mode) {
                    prefix <- paste0("clues/", rsid, "/", rsid, "-", model, "-", mode)
                    result <- fromJSON(paste0(prefix, ".json"))
                    clues_trajectory(prefix) %>%
                        mutate(rsid = rsid, model = model, mode = mode, pval = result$pval)
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
ages$mode <- factor(ages$mode, levels = c("mod_freq", "no_freq"), labels = c("With modern DAF", "Without modern DAF"))

traj$model <- factor(traj$model, levels = models, labels = str_replace_all(models, "_", " "))
ages$model <- factor(ages$model, levels = models, labels = str_replace_all(models, "_", " "))

# format the data for display
traj <- traj %>%
    mutate(
        significant = as.numeric(pval < (0.05 / length(deletions))),
        # format the p-value
        pval = ifelse(
            pval < 0.05,
            sprintf("%.2e", signif(pval, 2)),
            sprintf("%.2f", signif(pval, 2))
        ),
        snp_label = paste0(ifelse(rsid == "rs333", "CCR5Δ32", rsid), " (p=", pval, ")")
    )

# add the significant flag to the allele age df
ages <- ages %>%
    inner_join(
        traj %>% select(rsid, model, mode, significant) %>% unique(),
        by = join_by(rsid, model, mode)
    )

snp_colors <- c(
    "rs333" = "#33a02c",
    "rs61231801" = "#a6cee3",
    "rs66552573" = "#1f78b4",
    "rs67580019" = "#fb9a99",
    "rs143241023" = "#e31a1c",
    "rs150628438" = "#fdbf6f",
    "rs369842709" = "#ff7f00",
    "rs556322139" = "#cab2d6"
)

plt <- traj %>%
    # plot the trajectories
    ggplot(aes(x = epoch, y = freq, color = rsid, alpha = significant)) +

    # plot the maximum posterior trajectory
    geom_line(linewidth = 1, na.rm = TRUE) +

    # print the labels
    geom_dl(aes(label = snp_label), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.7), na.rm = TRUE) +

    # show the allele ages
    geom_point(aes(y = 0), ages) +

    # display as a grid
    facet_grid(model ~ mode, labeller = labeller(description = label_wrap_gen())) +

    # set the model colours
    scale_color_manual(values = snp_colors) +

    # plot non-significant trajectories as transparent
    scale_alpha(range = c(0.3, 1)) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, .20), breaks = seq(0, 1, 0.05), position = "left") +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 470))) +
    labs(
        title = "Ancestry stratified trajectories for the CCR5Δ32 deletion"
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

ggsave("figure/hapi_trajectories.png", plt, width = 9, height = 6)

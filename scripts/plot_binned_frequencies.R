#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.2
quiet(library(Hmisc))

# load the model calls for each deletion
data <- read_delim(fs::dir_ls(path = "data/results_8deletions_22_07_24/", glob = "*.tsv"), delim = "\t", id = "path", show_col_types = FALSE) %>%
    mutate(rsid = str_replace(str_replace(path, "data/results_8deletions_22_07_24/", ""), "_E.tsv", "")) %>%
    select(-path)

# load the sample metadata
samples <- read_tsv("data/neo.impute.1000g.sampleInfo.orig.tsv", guess_max = 100000, col_types = cols()) %>%
    select(sampleId, popId, groupAge, ageAverage)


geno <- data %>%
    pivot_longer(cols = c("No_filter", "Permissive_filter", "Strict_filter"), names_to = "model", values_to = "geno") %>%
    separate_longer_position(cols="geno", width=1) %>%
    # join the sample metadata
    inner_join(samples, by=c("Sample"="sampleId"))

# determine the sample size of the ancients
sample_size <- geno %>% pull(Sample) %>% unique() %>% length()

BIN_SIZE <- 1000

# group by age
binned <- geno %>%
    # focus on CCR5 in the recent past
    filter(rsid == "rs333" & ageAverage < 12000 & model != "No_filter") %>%
    # add a new age bin
    mutate(bin=-round(ageAverage/BIN_SIZE)*BIN_SIZE) %>%
    group_by(model, rsid, bin) %>%
    summarise(sum=sum(geno=="D"), count=n())

# calculate the binomial proportion confidence interval
binned <- bind_cols(binned, Hmisc::binconf(x=binned$sum, n=binned$count, alpha=0.05, return.df=TRUE)) %>%
    rename(daf=PointEst, ci_lower=Lower, ci_upper=Upper)

# constrain the extent of the plotting
xmin <- min(binned$bin)
xmax <- max(binned$bin)
xbreaks <- seq(xmin, xmax, BIN_SIZE)
xlabels <- round(xbreaks / BIN_SIZE)

size_min <- floor(min(binned$count) / 100) * 100
size_max <- ceiling(max(binned$count) / 100) * 100
size_breaks <- 2^(1:ceiling(log2(size_max)))

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

plt <- ggplot(binned, aes(x=bin, y=daf, color = rsid, weight=count)) +

    # geom_smooth(method = "loess", formula = "y ~ x", se = FALSE) +
    # geom_smooth(method="glm", method.args = list(family="quasipoisson"), formula = y ~ splines::ns(x, 3), se = FALSE) +
    geom_smooth(method="glm", method.args = list(family="binomial"), formula = y ~ x, se = FALSE) +
    geom_point(aes(size=count), alpha=.8, position=position_dodge(width = 200)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width=200, position=position_dodge(width = 200)) +

    facet_grid(rsid~model) +

    guides(size=guide_legend(title="Haploid count")) +

    # set the model colours
    scale_color_manual(values = snp_colors) +

    scale_y_continuous(breaks = seq(0, 1, 0.05), position = "left") +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    scale_size_continuous(limits  = c(size_min, size_max), breaks = size_breaks) +

    labs(
        title = "Binned frequency for CCR5Δ32 with binomial regression"
    ) +
    ylab("DAF") +
    xlab("kyr BP") +

    coord_cartesian(ylim=c(0, 0.35)) +

    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

ggsave("figure/binned_frequencies.png", plt, width = 9, height = 3)

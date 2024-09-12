#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse))

tags <- read_tsv(
    "data/ccr5delta32-tags.tsv",
    col_names = c("rsid", "chr", "pos", "ref", "alt", "sample", "geno", "ancestry"),
    show_col_types = FALSE
)

tags %>%
    separate_longer_delim(cols = c(geno, ancestry), delim = "|") %>%
    mutate(geno = as.numeric(geno), ancestry = as.numeric(ancestry)) %>%
    mutate(
        phase = row_number() %% 2,
        geno = ifelse(geno, alt, ref),
        ccr5 = case_when(
            rsid == "rs113341849" & geno == "A" ~ 1,
            rsid == "rs113010081" & geno == "C" ~ 1,
            rsid == "rs79815064" & geno == "G" ~ 1,
            rsid == "rs11574435" & geno == "T" ~ 1,
            .default = 0
        )
    ) %>%
    arrange(sample, phase) %>%
    group_by(sample, phase) %>%
    summarise(
        path = names(which.max(table(ancestry))),
        call = sum(ccr5) / n()
    )

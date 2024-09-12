#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.2

# load the models calls for each control deletion
data <- read_delim(fs::dir_ls(path = "data/results_8deletions_22_07_24/", glob = "*.tsv"), delim = "\t", id = "path", show_col_types = FALSE) %>%
    mutate(rsid = str_replace(str_replace(path, "data/results_8deletions_22_07_24/", ""), "_E.tsv", "")) %>%
    select(-path)

# load the sample metadata
meta.neo <- read_tsv("data/neo.impute.1000g.sampleInfo.orig.tsv", guess_max = 100000, col_types = cols()) %>%
    mutate(age = ifelse(groupAge == "Modern", 0, ageAverage))

# join the sample metadata
data <- data %>%
    inner_join(
        select(meta.neo, sampleId, age),
        by = c("Sample" = "sampleId")
    ) %>%
    mutate(gens = age / 28) %>%
    arrange(rsid, gens)

# --------------------------------------------------------------
# convert all the models into CLUES input format
# --------------------------------------------------------------

for (deletion in unique(data$rsid)) {

    # restrict to the current deletion
    locus <- data %>% filter(rsid == deletion)

    for (model in c("No_filter", "Permissive_filter", "Strict_filter")) {

        # convert the hard-called models into pseudo-likelihoods
        locus %>%
            mutate(call = case_when(
                get({{ model }}) == "RR" ~ "0.000000 -inf -inf",
                get({{ model }}) == "RD" ~ "-inf 0.000000 -inf",
                get({{ model }}) == "DD" ~ "-inf -inf 0.000000"
            )) %>%
            mutate(clues = sprintf("%.6f %s", gens, call)) %>%
            select(clues) %>%
            write_delim(paste0("clues/", deletion, "/", deletion, "-", model, ".ancient"), col_names = F, escape = "none", delim = "")
    }

    # convert the probabilities into log10-scaled likelihoods
    locus %>%
        mutate(call = sprintf("%.6f %.6f %.6f", log10(pRR_Data_n), log10(pRD_Data_n), log10(pDD_Data_n))) %>%
        mutate(clues = sprintf("%.6f %s", gens, call)) %>%
        select(clues) %>%
        write_delim(paste0("clues/", deletion, "/", deletion, "-p_Data.ancient"), col_names = F, escape = "none", delim = "")
}

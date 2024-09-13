#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse))

for (model in c("HAPI_samples", "West_Eurasian_samples")) {
    tags <- read_tsv(
        paste0("data/ccr5_tags-", model, ".tsv"),
        col_names = c("rsid", "chr", "pos", "ref", "alt", "sample", "geno", "ancestry"),
        show_col_types = FALSE
    )

    # call the CCR5 deletion using the tag SNPs, and
    calls <- tags %>%
        separate_longer_delim(cols = c(geno, ancestry), delim = "|") %>%
        mutate(geno = as.numeric(geno), ancestry = as.numeric(ancestry)) %>%
        mutate(
            phase = row_number() %% 2,
            geno = ifelse(geno, alt, ref),
            ancestry = case_when(
                ancestry == 1 ~ "ANA",
                ancestry == 2 ~ "CHG",
                ancestry == 3 ~ "WHG",
                ancestry == 4 ~ "EHG",
            ),
            # convert the SNP calls into support for the CCR5delta32 deletion
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
            # call consensus CCR5delta32 deletion and ancestry
            ancestry = names(which.max(table(ancestry))),
            call = round(mean(ccr5), 0),
            .groups = "drop"
        )

    # load the sample metadata
    meta.neo <- read_tsv("data/neo.impute.1000g.sampleInfo.orig.tsv", guess_max = 100000, col_types = cols()) %>%
        mutate(age = ifelse(groupAge == "Modern", 0, ageAverage))

    # join the sample metadata
    data <- calls %>%
        inner_join(
            select(meta.neo, sampleId, age),
            by = c("sample" = "sampleId")
        ) %>%
        mutate(gens = age / 28)

    # --------------------------------------------------------------
    # convert all the ancestry paths into CLUES input format
    # --------------------------------------------------------------

    # make the pan-ancestry input file
    data %>%
        group_by(sample, gens) %>%
        summarise(call = sum(call), .groups = "drop") %>%
        mutate(call = case_when(
            call == 0 ~ "0.000000 -inf -inf",
            call == 1 ~ "-inf 0.000000 -inf",
            call == 2 ~ "-inf -inf 0.000000"
        )) %>%
        mutate(clues = sprintf("%.6f %s", gens, call)) %>%
        arrange(gens) %>%
        select(clues) %>%
        write_delim(paste0("clues/ccr5_tags-", model, "/ccr5_tags-", model, "-ALL.ancient"), col_names = F, escape = "none", delim = "")

    # now do each of the ancestry paths
    for (anc in unique(data$ancestry)) {
        path <- data %>% filter(ancestry == anc)

        # convert the hard-called models into pseudo-likelihoods
        path %>%
            mutate(call = case_when(
                call == 0 ~ "0.000000 -inf",
                call == 1 ~ "-inf 0.000000",
            )) %>%
            mutate(clues = sprintf("%.6f %s", gens, call)) %>%
            arrange(gens) %>%
            select(clues) %>%
            write_delim(paste0("clues/ccr5_tags-", model, "/ccr5_tags-", model, "-", anc, ".ancient"), col_names = F, escape = "none", delim = "")
    }
}

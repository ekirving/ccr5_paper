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
quiet(library(stringr))

# load the CCR5 deletion calls
ccr5_calls <- read_tsv("data/ancient_samples_09_11_updated.tsv", col_types = cols()) %>%
  replace_na(list(Notes="")) %>%
  select(Sample, Predicted_by_model, Artifacts_filter, Minus_haplo_filter, pRR_Data, pRD_Data, pDD_Data, Notes)

# load the sample metadata
meta.neo <- read_tsv("data/neo.impute.1000g.sampleInfo.orig.tsv", guess_max = 100000, col_types = cols()) %>%
  mutate(age = ifelse(groupAge == "Modern", 0, ageAverage))

# join the sample metadata
ccr5_neo <- ccr5_calls %>%
  inner_join(
    select(meta.neo, sampleId, age),
    by = c("Sample" = "sampleId")
  ) %>%
  mutate(gens = age / 28) %>%
  arrange(gens)

# --------------------------------------------------------------
# convert the 5 different callsets into CLUES input format
# --------------------------------------------------------------

# the pRR_Data / pRD_Data / pDD_Data
ccr5_neo %>%
  mutate(call = sprintf("%.6f %.6f %.6f", log10(pRR_Data), log10(pRD_Data), log10(pDD_Data))) %>%
  mutate(clues = sprintf("%.6f %s", gens, call)) %>%
  select(clues) %>%
  write_delim("clues/p_Data/ccr5-p_Data.ancient", col_names = F, quote_escape = F, delim = "")

# Predicted_by_model
ccr5_neo %>%
  mutate(call = case_when(
    Predicted_by_model == "RR" ~ "0.000000 -inf -inf",
    Predicted_by_model == "RD" ~ "-inf 0.000000 -inf",
    Predicted_by_model == "DD" ~ "-inf -inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", gens, call)) %>%
  select(clues) %>%
  write_delim("clues/Predicted_by_model/ccr5-Predicted_by_model.ancient", col_names = F, quote_escape = F, delim = "")

# Artifacts_filter
ccr5_neo %>%
  mutate(call = case_when(
    Artifacts_filter == "RR" ~ "0.000000 -inf -inf",
    Artifacts_filter == "RD" ~ "-inf 0.000000 -inf",
    Artifacts_filter == "DD" ~ "-inf -inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", gens, call)) %>%
  select(clues) %>%
  write_delim("clues/Artifacts_filter/ccr5-Artifacts_filter.ancient", col_names = F, quote_escape = F, delim = "")

# Minus_haplo_filter
ccr5_neo %>%
  mutate(call = case_when(
    Minus_haplo_filter == "RR" ~ "0.000000 -inf -inf",
    Minus_haplo_filter == "RD" ~ "-inf 0.000000 -inf",
    Minus_haplo_filter == "DD" ~ "-inf -inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", gens, call)) %>%
  select(clues) %>%
  write_delim("clues/Minus_haplo_filter/ccr5-Minus_haplo_filter.ancient", col_names = F, quote_escape = F, delim = "")

# minus_haplotype
ccr5_neo %>%
    mutate(call = ifelse(Notes == "Probably minus haplotype", "-inf 0.000000 -inf", "0.000000 -inf -inf")) %>%
    mutate(clues = sprintf("%.6f %s", gens, call)) %>%
    select(clues) %>%
    write_delim("clues/minus_haplotype/ccr5-minus_haplotype.ancient", col_names = F, quote_escape = F, delim = "")

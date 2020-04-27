#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")

options(echo=TRUE)
options(verbose=TRUE)
args = commandArgs(trailingOnly=TRUE)

library("stringr")
library("dplyr")
library("readr")
library("VariantAnnotation")
library("data.table")
library("parallel")


cosmic_reference_fread_col_specs <- c(
  `Gene name` = "character",
  `Accession Number` = "character",
  `Gene CDS length` = "double",
  `HGNC ID` = "double",
  `Sample name` = "character",
  ID_sample = "double",
  ID_tumour = "double",
  `Primary site` = "character",
  `Site subtype 1` = "character",
  `Site subtype 2` = "character",
  `Site subtype 3` = "character",
  `Primary histology` = "character",
  `Histology subtype 1` = "character",
  `Histology subtype 2` = "character",
  `Histology subtype 3` = "character",
  `Genome-wide screen` = "character",
  GENOMIC_MUTATION_ID = "character",
  LEGACY_MUTATION_ID = "character",
  MUTATION_ID = "double",
  `Mutation CDS` = "character",
  `Mutation AA` = "character",
  `Mutation Description` = "character",
  `Mutation zygosity` = "character",
  LOH = "logical",
  GRCh = "double",
  `Mutation genome position` = "character",
  `Mutation strand` = "character",
  SNP = "character",
  `Resistance Mutation` = "character",
  `FATHMM prediction` = "character",
  `FATHMM score` = "double",
  `Mutation somatic status` = "character",
  Pubmed_PMID = "double",
  ID_STUDY = "double",
  `Sample Type` = "character",
  `Tumour origin` = "character",
  Age = "double",
  HGVSP = "character",
  HGVSC = "character",
  HGVSG = "character"
)
#arg_comsic_ref <- "/scratch/shahlab_tmp/sbeatty/ind231/reference_data/subset_cosmic_coordinates_added.csv"
#args <- c("--cosmic_reference_data=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/subset_cosmic_coordinates_added.csv", "--cpu_count=40", "--output_file=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/subset_cosmic_coordinates_added_mutation_types_parsed.csv")


cosmic_reference_data <- args[str_detect(args,"--cosmic_reference_data")] %>% str_replace("--cosmic_reference_data=","")
output_file <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
cpu_count <- args[str_detect(args,"--cpu_count")] %>% str_replace("--cpu_count=","")

ref_df <- fread(paste(cosmic_reference_data), colClasses=cosmic_reference_fread_col_specs, data.table=FALSE)
patterns_to_remove <-'([c.0-9-ATCG>+-_/*])'
cosmic_mutation_types <- mcmapply(ref_df[,"Mutation.CDS"],USE.NAMES=FALSE, FUN=function(x){str_replace_all(x, patterns_to_remove, "")}, mc.cores=cpu_count)
#cosmic_mutation_types <- ref_df[,"Mutation.CDS"] %>% str_replace_all(patterns_to_remove, "")
cosmic_mutation_types[which(cosmic_mutation_types == "")] <- "snv"
ref_df[,"parsed_mutation_type"] <- cosmic_mutation_types
fwrite(data.frame(ref_df), file=paste(output_file))

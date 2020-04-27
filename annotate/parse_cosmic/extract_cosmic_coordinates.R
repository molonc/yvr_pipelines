#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")

options(echo=TRUE)
options(verbose=TRUE)
args = commandArgs(trailingOnly=TRUE)

library("stringr")
#library("dplyr")
#library("readr")
#library("VariantAnnotation")
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
  MUTATION_ID = "double"(),
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

#args <- c("--random", "--ncpus=20", "--ref_file=/scratch/shahlab_tmp/sbeatty/reference/CosmicMutantExport.tsv", "--outputfile=/scratch/shahlab_tmp/sbeatty/reference/coordinates.csv")

target_file <- args[str_detect(args,"--ref_file")] %>% str_replace("--ref_file=","")
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
cpu_count <- args[str_detect(args,"--cpu_count")] %>% str_replace("--cpu_count=","")

setDTthreads(as.numeric(cpu_count))
ref_df <- fread(paste(target_file), colClasses=cosmic_reference_fread_col_specs, data.table=FALSE)
output <- ref_df[,"Mutation genome position"]
fwrite(data.frame(output), file=paste(output_file))

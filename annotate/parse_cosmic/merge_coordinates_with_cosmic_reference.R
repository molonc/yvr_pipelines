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

#args <- c("--random", "--chr_coordinates=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates_chr.csv", "--start_coordinates=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates_start.csv", "--end_coordinates=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates_end.csv", "--cosmic_reference_file=/scratch/shahlab_tmp/sbeatty/reference/cosmic_sample.tsv", "--output_file=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_parsed.csv")
chr_coordinates <- args[str_detect(args,"--chr_coordinates")] %>% str_replace("--chr_coordinates=","")
start_coordinates <- args[str_detect(args,"--start_coordinates")] %>% str_replace("--start_coordinates=","")
end_coordinates <- args[str_detect(args,"--end_coordinates")] %>% str_replace("--end_coordinates=","")
cosmic_reference_file <- args[str_detect(args,"--cosmic_reference_file")] %>% str_replace("--cosmic_reference_file=","")
output_file <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")

chr <- fread(paste(chr_coordinates), colClasses=c("character"), data.table=FALSE, col.names="chr")
start <- fread(paste(start_coordinates), colClasses=c("double"), data.table=FALSE, col.names="start")
end <- fread(paste(end_coordinates), colClasses=c("double"), data.table=FALSE, col.names="end")

coordinate_columns_df <- cbind(chr,start,end)
rm(chr,start,end)

ref_df <- fread(paste(cosmic_reference_file), colClasses=cosmic_reference_fread_col_specs, data.table=FALSE)


output <- cbind(coordinate_columns_df,ref_df) %>% data.frame
fwrite(output, file=paste(output_file))





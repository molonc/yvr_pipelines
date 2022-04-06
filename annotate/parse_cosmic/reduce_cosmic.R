#!/usr/bin/env Rscript


options(echo=TRUE)
options(verbose=TRUE)
args = commandArgs(trailingOnly=TRUE)

library("stringr")
library("dplyr")
library("readr")
library("VariantAnnotation")
library("data.table")
#library("tidyr")

cosmic_reference_fread_col_specs <- c(
	`chr` = "character",
	`start` = "double",
	`stop`= "double",
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
target_file <- args[str_detect(args,"--cosmic_reference_file")] %>% str_replace("--cosmic_reference_file=","")
output_file <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
cpu_count <- args[str_detect(args,"--cpu_count")] %>% str_replace("--cpu_count=","")

if(length(args) == 0){
  target_file <- "//projects/molonc/aparicio_lab/sbeatty/reference/cosmic_coordinates_added_mutation_types_parsed.csv"
  output_file <- "//projects/molonc/aparicio_lab/sbeatty/reference/cosmic_coordinates_added_mutation_types_parsed_reduced.csv"
  cpu_count <- 35
}
setDTthreads(threads=cpu_count)
cosmic_data <- fread(target_file,colClasses=cosmic_reference_fread_col_specs, nThread=35)
cosmic_data <-  filter(cosmic_data, chr != "")
patterns_to_remove <-'([0-9a-z:.]|[^><AGCT])'

mutations <- mcmapply(cosmic_data[,"HGVSG"],USE.NAMES=FALSE, FUN=function(x){str_replace_all(x, patterns_to_remove, "")}, mc.cores=cpu_count)

split_variant <- function(x, ref_or_alt){
  output <- NA
  if(str_detect(x, ">") == TRUE){
     split_variant_annotation <- str_split(x,">") %>% unlist()
     if(length(split_variant_annotation) == 2){
      if(ref_or_alt == "ref"){
        output <- split_variant_annotation[1]
      } else if (ref_or_alt == "alt"){
        output <- split_variant_annotation[2]
      } else {
        output <- NA
      }
     } else {
        output <- NA
      }
}
output
}
refs <- mcmapply(mutations,USE.NAMES=FALSE, FUN=function(x){split_variant(x, "ref")}, mc.cores=cpu_count)
alts <- mcmapply(mutations,USE.NAMES=FALSE, FUN=function(x){split_variant(x, "alt")}, mc.cores=cpu_count)

cosmic_data[,"ref"] <- refs
cosmic_data[,"alt"] <- alts
cosmic_data_reduced <- cosmic_data %>% group_by(chr, start, end, Mutation.strand, GENOMIC_MUTATION_ID, Genome.wide.screen, parsed_mutation_type, ref, alt,  FATHMM.prediction, FATHMM.score)
#rm(cosmic_data)
cosmic_data_reduced <- cosmic_data_reduced %>%  summarize(n=n())
cosmic_data_reduced <- cosmic_data_reduced %>%  pivot_wider(names_from=Genome.wide.screen, values_from = n) 
cosmic_data_reduced <- cosmic_data_reduced %>% data.frame
cosmic_data_reduced <- cosmic_data_reduced[!is.na(cosmic_data_reduced$start),]
cosmic_data_reduced[c(cosmic_data_reduced$chr == 25),"chr"] <- "MT"
cosmic_data_reduced[c(cosmic_data_reduced$chr == 24),"chr"] <- "Y"
cosmic_data_reduced[c(cosmic_data_reduced$chr == 23),"chr"] <- "X"
cosmic_data_reduced[is.na(cosmic_data_reduced$n), "n"] <- 0
cosmic_data_reduced[is.na(cosmic_data_reduced$Y), "y"] <- 0
cosmic_data_reduced[,"cosmic_patients_total"] <- cosmic_data_reduced[,"n"] + cosmic_data_reduced[,"y"]
cosmic_data_reduced[,"cosmic_patients_gwas"] <- cosmic_data_reduced["y"]
cosmic_data_reduced <- cosmic_data_reduced[,c("chr","start","end", "FATHMM.prediction", "FATHMM.score", "Mutation.strand", "GENOMIC_MUTATION_ID", "parsed_mutation_type", "ref", "alt","cosmic_patients_total", "cosmic_patients_gwas")]
fwrite(cosmic_data_reduced, file=paste(output_file))




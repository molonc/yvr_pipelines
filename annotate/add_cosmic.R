#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")

options(echo=TRUE)
options(verbose=TRUE)
args = commandArgs(trailingOnly=TRUE)

require("stringr")
require("dplyr")
require("readr")
require("VariantAnnotation")
require("data.table")


cols_tsv <- 
cols(
  .default = col_character(),
  case_id = col_character(),
  normal_id = col_character(),
  tumour_id = col_character(),
  chromosome = col_character(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  gene_id = col_character(),
  type = col_character(),
  filter = col_character(),
  ref = col_character(),
  alt = col_character(),
  gt = col_character(),
  pl = col_character(),
  mut_pr = col_double(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  dbsnp = col_character(),
  thousand_genomes = col_character(),
  cosmic = col_character(),
  caller = col_character(),
  trinucleotide_ref = col_character(),
  trinucleotide_alt = col_character(),
  amino_acid_change = col_character(),
  functional_class = col_character(),
  gene_coding = col_character(),
  project = col_character(),
  substitution = col_character(),
  impact = col_character(),
  dna_change = col_character(),
  biotype = col_character()
)

csv_column_types_for_fread <- 
c(
  case_id = "character",
  normal_id = "character",
  tumour_id = "character",
  chromosome = "character",
  start = "double",
  stop = "double",
  gene = "character",
  gene_id = "character",
  type = "character",
  filter = "character",
  ref = "character",
  alt = "character",
  gt = "character",
  pl = "character",
  mut_pr = "double",
  tr = "double",
  ta = "double",
  nr = "double",
  na = "double",
  dbsnp = "character",
  thousand_genomes = "character",
  cosmic = "character",
  caller = "character",
  trinucleotide_ref = "character",
  trinucleotide_alt = "character",
  amino_acid_change = "character",
  functional_class = "character",
  gene_coding = "character",
  project = "character",
  substitution = "character",
  impact = "character",
  dna_change = "character",
  biotype = "character"
)


source("/scratch/shahlab_tmp/sbeatty/yvr_pipelines/annotate/core_functions.R")

target_file <- args[str_detect(args,"--targetfile")] %>% str_replace("--targetfile=","")
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
cosmic_reference_file <- args[str_detect(args,"--cosmic_reference_file")] %>% str_replace("--cosmic_reference_file=","")
if(length(args) == 0){
  target_file <- "clinvar_added/GERM_SNV_INTERSECT_SA1228N.clinvar.csv"
  output_file <- "test.csv"
  cosmic_reference_file <- "/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates_added_mutation_types_parsed_reduced.csv"#
  ncpus <- 20#
} 
type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()
gencode_reference_col_types = cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_character(),
  X4 = col_double(),
  X5 = col_double(),
  X6 = col_character(),
  X7 = col_character(),
  X8 = col_character(),
  X9 = col_character()
)

gene_from_add_info <- function(x){
    output <- NA
    add_info_type <- unlist(strsplit(x, ";")) %>% trimws()
    matches <- str_detect(add_info_type, "gene_name")
    if(length(matches == 1)){
      output <- add_info_type[str_detect(add_info_type, "gene_name")] %>% str_replace("gene_name ","")
    }
    output 
  } 


joined_df <- read_csv(paste(target_file), col_types=cols_tsv) %>% data.frame
joined_ranges <- makeGRangesFromDataFrame(data.frame(joined_df), seqnames.field ="chr" , start.field='start', end.field='end', keep.extra.columns=TRUE)
cosmic_reference_data <- fread(cosmic_reference_file, data.table=FALSE, colClasses=csv_column_types_for_fread)
cosmic_reference_data <- cosmic_reference_data[!is.na(cosmic_reference_data$chr),]
cosmic_ranges <- makeGRangesFromDataFrame(data.frame(cosmic_reference_data),seqnames.field ="chr" , start.field='start', end.field='end')

overlaps <- findOverlaps(joined_ranges, cosmic_ranges, type="equal")

if(length(overlaps) > 0){
overlaps <- overlaps %>% data.frame
overlaps <- extract_exact_matches(matches_df =overlaps, query_df=joined_df, query_ref_allele_column="ref", query_alt_allele_column="alt", subject_df=cosmic_reference_data,subject_ref_allele_column="ref",  subject_alt_allele_column="alt", subject_strand_column="Mutation.strand", query_strand_column="genecode_strand")


  cosmic_reference_data[,"unique_id"] <- NA
  #joined_df <- joined_ranges %>% data.frame
  joined_df[,"unique_id"] <- paste0(joined_df$chr,"@",joined_df$end)
  cosmic_reference_data$unique_id[overlaps$subjectHits] <- joined_df$unique_id[overlaps$queryHits]
  output_df <- left_join(joined_df, cosmic_reference_data, by="unique_id", suffix=c(".variantcaller","cosmic"), KEEP=TRUE)
  names(output_df) <- substitute_column_name("seqnames.variantcaller","chr" ,output_df)
  names(output_df) <- substitute_column_name("start.variantcaller","start" ,output_df)
  names(output_df) <- substitute_column_name("end.variantcaller","end" ,output_df)
  names(output_df) <- substitute_column_name("seqnames","chr" ,output_df)
 } else {
  output_df <- joined_df
 }

fwrite(output_df, file=output_file)



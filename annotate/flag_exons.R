#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")

options(echo=TRUE)
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)

require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)
require("parallel", quietly=TRUE)
require("VariantAnnotation", quietly=TRUE)
require("data.table", quietly=TRUE)
options(echo=TRUE)
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
 
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

substitute_column_name <- function(in_name,out_name, input_data){
  target_data_column_names <- names(input_data)
  matches <- which(target_data_column_names == in_name)
  if(length(matches) == 1){
    target_data_column_names[matches] <- out_name
  }
  target_data_column_names
}


gencode_column_types_for_fread <- 
c("character", "character", "character", "double", "double", "character" ,"character", "character", "character"
)

#args <- c("--random", "--targetfile=/scratch/shahlab_tmp/sbeatty/ind231/genecode_clinvar_annotated/GERM_SNV_INTERSECT_SA1228N.clinvar.gencode.csv", "--outputfile=genecode_comsic_snpeff_annotated/add_exon_test.csv", "--chr_name=chr", "--start_name=start", "--end_name=end")

if(length(args) == 0){
args <- "Rscript //scratch/shahlab_tmp/sbeatty/yvr_pipelines/annotate/flag_exons.R --target_file=cosmic_added/GERM_STRE_INDEL_SA1228T.clinvar.cosmic.csv --outputfile=exons_flagged/GERM_STRE_INDEL_SA1228T.clinvar.cosmic.exons.csv"
args <- str_split(args," ") %>% unlist
}
target_file <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
#target_file <- "cosmic_added/cosmic_added/GERM_SNV_INTERSECT_SA1228N.clinvar.cosmic.csv"
target_file_type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()

if(length(which(str_detect(args, "--chr_name") == TRUE)) == 1){
  chr_name <- args[str_detect(args,"--chr_name")] %>% str_replace("--chr_name=","")
} else {
  chr_name <- "chr"
}

if(length(which(str_detect(args, "--start_name") == TRUE)) == 1){
  start_name <- args[str_detect(args,"--start_name")] %>% str_replace("--start_name=","")
} else {
  start_name <- "start"
}

if(length(which(str_detect(args, "--end_name") == TRUE)) == 1){
  end_name <- args[str_detect(args,"--end_name")] %>% str_replace("--end_name=","")
} else {
  end_name <- "end"
}



# load reference data

gencode_reference <- fread("//projects/molonc/aparicio_lab/sbeatty/reference/gencode.v19.annotation.gtf_withproteinids", skip=5, data.table=FALSE,
  col.names=c("chromosome", "source", "feature_type","start","stop", "score", "strand", "phase", "add_info"), colClasses=gencode_column_types_for_fread)
#gencode_reference_exons <- filter(gencode_reference, feature_type=="exon")


gencode_reference$chromosome <- gsub("chr","", gencode_reference$chromosome)
gencode_reference$add_info <- gsub(regex('\\"'),"",gencode_reference$add_info)
gencode_ranges <-  makeGRangesFromDataFrame(data.frame(gencode_reference), 
         seqnames.field ="chromosome" , start.field="start", end.field="stop", 
      keep.extra.columns=TRUE)


# load sample data 
sample_data_df <- fread(paste(target_file), data.table=FALSE, skip=0)
#names(sample_data_df) <- gsub("chromosome", "chr", names(sample_data_df))
#names(sample_data_df) <- gsub("stop", "end", names(sample_data_df))
sample_data_df[,"unique_id"] <- paste0(sample_data_df$chr,"@",sample_data_df$end)

name_matches <- which(names(sample_data_df) == "chr")
if(length(name_matches) > 1){
match_vector <- names(sample_data_df) == "chr"
match_vector[name_matches[1]] <- FALSE
sample_data_df <- sample_data_df[,names(sample_data_df)[!c(match_vector)]]
} else {
  names(sample_data_df) <- gsub("chr.variantcaller","chr", names(sample_data_df))
}

name_matches <- which(names(sample_data_df) == "start")
if(length(name_matches) > 1){
match_vector <- names(sample_data_df) == "start"
match_vector[name_matches[1]] <- FALSE
sample_data_df <- sample_data_df[,names(sample_data_df)[!c(match_vector)]]
}

name_matches <- which(names(sample_data_df) == "end")
if(length(name_matches) > 1){
match_vector <- names(sample_data_df) == "end"
match_vector[name_matches[1]] <- FALSE
sample_data_df <- sample_data_df[,names(sample_data_df)[!c(match_vector)]]
}

sample_data_ranges <- makeGRangesFromDataFrame(sample_data_df, 
       seqnames.field =chr_name , start.field=start_name, end.field=end_name, keep.extra.columns=TRUE)

overlaps <- findOverlaps(sample_data_ranges, gencode_ranges, type="any")
sample_data_df[,"gencode_feature_type"] <- NA
sample_data_df[,"gencode_add_info"] <- NA

if(length(overlaps) > 0){
  #sample_data_df[c(overlaps@from %>% unique),"is_exon"] <- TRUE
  sample_data_df[,"gencode_feature_type"][overlaps@from] <- gencode_reference$feature_type[overlaps@to]
  sample_data_df[,"gencode_add_info"][overlaps@from] <- gencode_reference$add_info[overlaps@to]

 } 
fwrite(sample_data_df, file=output_file)

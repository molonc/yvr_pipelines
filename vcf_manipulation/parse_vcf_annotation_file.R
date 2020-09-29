#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)

require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)
require()
#require("parallel", quietly=TRUE)
#require("VariantAnnotation", quietly=TRUE)
#require("data.table", quietly=TRUE)
#options(echo=TRUE)

args = commandArgs(trailingOnly=TRUE)

#target_vcf_path <- "/scratch/shahlab_tmp/sbeatty/ind231/data_export/vcf_subset.vcf"
#output_file_name <- "vcf_output.vcf"
#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
##line_to_add <- "##INFO=<ID=FATHMM,Number=.,Type=String,Description=\"FATHMM score\">"
#new_header_line_ID <- "FATHMM,Number"
#new_header_line_description <- "FATHMM score"
#new_annotation_column <- "FATHMM"


if(length(args) == 0){
args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/parse_vcf_annotation_file.R --target_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/ind231_summary_V14.SA1228.annotations_subset.csv --output_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/ind231_summary_V14.parsed_annotations.tsv.gz"
args <- str_split(args," ") %>% unlist
}

target_file <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")

annotation_input_file.column_specs <- cols(
  .default = col_character(),
  TB_ID = col_double(),
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  task = col_double(),
  mut_pr = col_double(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  Clinvar_AlleleID = col_double(),
  GeneID = col_double(),
  NumberSubmitters = col_double(),
  Clinvar_VariationID = col_double(),
  cosmic_patients_total = col_double(),
  cosmic_patients_gwas = col_double(),
  stop_codon_gain_potential = col_logical(),
  stop_codon_loss_potential = col_logical(),
  CCDS = col_logical(),
  FATHMM.score = col_double()
)


annotation_input_file <- read_csv(target_file, col_types=annotation_input_file.column_specs)
annotation_input_file_columns <- names(annotation_input_file)
#names(annotation_input_file)[1] <- paste("#", names(annotation_input_file)[1], sep="")

temporary_annotations_file_name <-   target_file %>% str_split("/") %>% unlist %>% last %>% str_replace(".csv",".temp.tsv")
write_tsv(annotation_input_file, path=temporary_annotations_file_name, col_names=FALSE)
system(paste0("bgzip -c ", temporary_annotations_file_name, " > ", output_file))
system(paste0("tabix -s1 -b2 -e3 ", output_file))

/scratch/shahlab_tmp/sbeatty/ind231/ind231_summary_V14.csv


target_file %>% str_split("/") %>% unlist %>% last
target_columsn <- gencode_gene_name


##INFO=<ID=1000Gen,Number=.,Type=String,Description="1000Gen flag">

type = 
number = 
ID = 
Description

museq_col_specs <- cols(
  case_id = col_character(),
  normal_id = col_character(),
  tumour_id = col_character(),
  chromosome = col_double(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  gene_id = col_character(),
  type = col_character(),
  filter = col_character(),
  ref = col_character(),
  alt = col_character(),
  gt = col_character(),
  pl = col_number(),
  mut_pr = col_double(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  dbsnp = col_character(),
  thousand_genomes = col_character(),
  cosmic = col_logical(),
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
a <- read_tsv("/scratch/shahlab_tmp/danlai/APARICIO-590/SING_MUTA/mutationseq_pipeline/OUTPUT/RUN/SA1228N_museq/outputs/results/TASK_8_PARSE_museq_parsed.tsv", col_types=cols(chromosome = col_character(), cosmic=col_character()))

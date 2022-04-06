#!/usr/bin/env Rscript

options(echo=TRUE)
options(verbose=TRUE)
args = commandArgs(trailingOnly=TRUE)

library("stringr")
library("dplyr")
library("readr")
#library("VariantAnnotation")
library("data.table")
library("parallel")


.cosmic_pos_split <- function(x, return_value){
  output <- x %>% str_split(":") %>% unlist
  if(return_value == "chr"){
    output <- output %>% first
  } else if (return_value == "start"){
    output <- output[2] %>% str_split("-") %>% unlist %>% first
  } else if (return_value == "end"){
    output <- output[2] %>% str_split("-") %>% unlist  %>% last
  } else {
    output <- "error"
  }
  output
}

#args <- c("--input_file=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates.csv", "--outputfile=/scratch/shahlab_tmp/sbeatty/ind231/reference_data/cosmic_coordinates_chr.csv", "--output_coordinate=start")
input_file <- args[str_detect(args,"--input_file")] %>% str_replace("--input_file=","")
#!str_detect(paste(target_file), paste(getwd()))
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
output_coordinate <- args[str_detect(args,"--output_coordinate")] %>% str_replace("--output_coordinate=","")
cpu_count <- args[str_detect(args,"--ncpus")] %>% str_replace("--ncpus=","")


input <- fread(paste(input_file), data.table=FALSE, col.names="position")
output <- input$position 
output[1:length(input$position)] <- mcmapply(input$position, FUN=function(x){.cosmic_pos_split(x, output_coordinate)}, mc.cores=cpu_count)
output <- data.frame(output)
names(output) <- paste(output_coordinate)
fwrite(output,file=paste(output_file))
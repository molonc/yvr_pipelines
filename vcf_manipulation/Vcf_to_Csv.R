options(echo=TRUE)
options(verbose=FALSE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)
require("tidyr", quietly=TRUE)



req_libs <- c("stringr", "readr", "dplyr","tidyr") 
for (lib in req_libs) {
  if (!require(lib, character.only = TRUE)) {
    cat("Installing missing library, may take some time, please be patient just this once")
      install.packages(lib)
  }
  if (!require(lib, character.only = TRUE)) {
    stop("Error installing required packages, please contact Daniel Lai at dalai@bccrc.ca to figure this out...")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
	args <-  "Rscript scripts/Vcf_to_Csv.R --target_file=277845T2.PAIR_MUTA.vcf"
	args <- args %>% str_split(" ") %>% unlist

}

target_file <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")

target_file_per_line <- read_lines(target_file)

tsv <- read_tsv(target_file,skip=which(str_detect(target_file_per_line, "#CHROM") == TRUE)[1]-1, col_types=cols('#CHROM' = col_character()))

spread_row <- function(row_i){
	info_tag_value_and_names <-  row_i[,"INFO"] %>% str_split(";") %>% unlist
	info_tag_names <- mapply(info_tag_value_and_names, USE.NAMES=FALSE, FUN=function(x){str_split(x,"=")[[1]][1]})
	info_tag_values <- mapply(info_tag_value_and_names, USE.NAMES=FALSE, FUN=function(x){str_split(x,"=")[[1]][2]})

	if(length(which(names(tsv) == "FORMAT")) > 0){
		spread_row <- tidyr::separate(row_i,col=last(names(row_i)),sep=":", into=unlist(str_split(row_i[1,"FORMAT"],":")), remove=TRUE) %>% dplyr::select(-FORMAT)  %>%  tidyr::separate(.,col="INFO",sep=";", into=info_tag_names, remove=TRUE) 
	} else {
		spread_row <- tidyr::separate(row_i,col="INFO",sep=";", into=info_tag_names, remove=TRUE)
	}
	for(tag in info_tag_names){
		spread_row[,tag] <- spread_row[,tag]  %>% str_replace_all(.,paste0(tag,"="),"")
	}
	spread_row
}


if(nrow(tsv) > 0){
	spread_rows <- lapply(1:nrow(tsv), FUN=function(x){spread_row(tsv[x,])})
	output_dataframe <- bind_rows(spread_rows)
} else {
	output_dataframe <- data.frame(message="No variants located within the target file")
}


write_csv(output_dataframe, file=str_replace(last(unlist(str_split(target_file,"/"))),".vcf",".csv"))

write_tsv(data.frame(readLines(target_file,n=which(str_detect(target_file_per_line, "#CHROM") == TRUE)[1]-1)), file=str_replace(last(unlist(str_split(target_file,"/"))),".vcf","_header.txt"), col_names=FALSE)

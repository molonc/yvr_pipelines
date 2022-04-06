
require("readr")
require("stringr")
require("dplyr")


test_file <- c("//projects/molonc/aparicio_lab/sbeatty/BXE/BXE-264/files_task6/SA575T.PAIR_STRE_INDEL.TASK_10_STRELKA_FLAG_COSMIC_INDL_strelka.passed.somatic.indels.annotSnpEff.annotMA.flagDBsnp.flag1000gen.flagCosmic.vcf")

temp_header_path <- test_file %>% str_split("/") %>% unlist %>% last %>% paste0("./",.,".temp.hdr")

command <- paste("bcftools view -h",test_file, ">",temp_header_path)

print(command)

system(command)

header <- readLines(temp_header_path)
info_lines <- header[str_detect(header,"##INFO=<")]

parse_info_line <- function(info_line_string){
	output <- info_line_string %>% str_split(",") %>% unlist %>% first %>% str_split("=") %>% unlist %>% last
	output
}


info_tags <- mapply(info_lines, USE.NAMES=FALSE, FUN=parse_info_line)



standard_tags <- header %>% last %>% str_split("\t") %>% unlist %>% str_replace("#","")


tags <- c(info_tags, standard_tags)

c <-  tags %>% paste("-F",.) %>% paste(collapse=" ")

a <- "-F PR -F TC -F TR -F TA -F NR -F NA -F ND -F NI -F ANN -F LOF -F NMD -F MA -F DBSNP -F 1000Gen -F Cosmic -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO"


# sing muta 

"-F PR -F TC -F TR -F TA -F NR -F NA -F ND -F NI -F GT -F PL -F ANN -F LOF -F NMD -F MA -F DBSNP -F 1000Gen -F Cosmic -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO"

# sing stre

"-F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -F MQ -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO"
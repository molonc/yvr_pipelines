### Realign: A snakemake workflow to convert bam files from BWA-ALN to BWA-MEM 

Workflow:

1. List of bam files specified as a list
2. Sort by query name using samtools
3. Convert bam to fastq using bedtools
4. Align using bwa-mem
5. Coordinate sort using samtools
6. Mark duplicates using picard tools
7. Index using samtools. 

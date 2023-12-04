##Chromatin Workshop

#Location
/grid/genomicscourse/home/shared/chromatin_workshop

##github
https://github.com/AndressaOL/CSHL_course

##
#FIRST STEP - MultiQC (1) Quality Control of Sequencing using FastQC/MultiQC)
data from folder fastq

#SECOND STEP - 2) Processing data & Genome Mapping

data from folders: genome & fastq
data fastq: SRR5063143_naive_H3K27ac_PE.fastq
data genome:  index_chrm22  & Homo_sapiens.GRCh38.dna.chromosome.22.fa

##
2.3) Post-mapping data
genome folder : chrom22.sizes

#THIRD STEP: 3) Peak Calling

bam files: folders bam_files
SRR5063149_naive_H3K4me3_treat.bam
SRR5063154_naive_input_treat.bam
SRR5063143_naive_H3K27ac_treat.bam
SRR5063153_naive_input_treat.bam
#

#peak_data
#peak files (outputs)
H3K27ac_peaks.broadPeak #peak
H3K27ac_peaks.broadPeak_result.txt #frips
H3K4me3_peaks.narrowPeak
H3K4me3_peaks.narrowPeak_result.txt
##

#folder bw_data
Plots: 1) deeptools:Correlation and Heatmap plots
SRR5063143_naive_H3K27ac_norm.bw
SRR5063144_naive_H3K27ac_norm.bw
SRR5063149_naive_H3K4me3_norm.bw
SRR5063150_naive_H3K4me3_norm.bw
SRR5063153_naive_input_norm.bw
SRR5063154_naive_input_norm.bw
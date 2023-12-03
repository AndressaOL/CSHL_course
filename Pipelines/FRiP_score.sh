#!/bin/bash
#***************************************************
##calculating the FRiP Score
##Modified from https://yiweiniu.github.io/blog/2019/03/Calculate-FRiP-score/
##use chmod 777 or chmod +x to give acess to execute the pipeline
#$1 bam 
#2 Peak file (Bed)
##Author: Andressa Oliveira de Lima aolima@uw.edu
##*************************************************


##First step 
#Count total reads in BAM
samtools view -c  $1 > $1_samtools.txt

## Second Step
##convert BAM (BAM used to call peaks) to BED
bedtools bamtobed -i $1 | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'  > $1.tagAlign.bed
bedtools sort -i  $1.tagAlign.bed >   $1.sorted.tagAlign.bed 

#Third Step
##Count the overlap regions peaks and bam file 
bedtools sort -i $2 | bedtools merge -i stdin |\
bedtools intersect -u -a $1.sorted.tagAlign.bed -b stdin | wc -l > $2_total.txt

##Calculate the FRiPs
v1=$(<"$2_total.txt")
v2=$(<"$1_samtools.txt")
 #result=$(($v1 / $v2))
result=$(awk "BEGIN {print $v1/$v2}")
echo "$result" > "$2_result.txt"

	
##delete files
rm *samtools.txt 
rm *tagAlign.bed 
rm *total.txt

#!/bin/bash
####################

#modify the bed file format to use as input on GREAT
for i in *Peak; do
awk  -v OFS='\t' '{print "chr"$1,$2,$3}' $i > $i.bed
done

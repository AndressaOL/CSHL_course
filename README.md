# CSHL_course

## CSHL Fall 2023 Computational Genomics Course
Welcome! This is the material and tutorials for the Chromatin Workshop.
The database adopted in this course is under the reference: "Enhancer Chromatin and 3D Genome Architecture Changes from Naive to Primed Human Embryonic Stem Cell States".
https://www.cell.com/stem-cell-reports/pdfExtended/S2213-6711(19)30126

## Data info:
- **Cell Type**: Naïve cells
- **Histone modification**: H3K4me3: Promoters & H3K27ac: Promoters and Enhancers
- **Library info**: **1)** SE-fastq files; **2)** 3 M rads
- **Data File**: Total 6 files (2 for each histone modification & 2 input files)

### 1) Quality Control of Sequencing using FastQC/MultiQC
- **Documentation**: *FastQC*: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ & *MultiQC*: https://multiqc.info/
- FastQC: ~3 -4 min <br /> 
`fastqc *.fastq`
- MultiQC: ~ 1 min <br />
`multiqc *.zip`

### 2) Processing data & Genome Mapping
**Chromap** for aligning and preprocessing high throughput chromatin profiles (*ATAC-seq & ChIP-seq*); **1)** Trimm the low-quality reads and adaptors; **2)** Remove duplicated reads; **3)** Perform the mapping. https://github.com/haowenz/chromap

- Build the indexed genome (available !!!) ~ 1 min
` chromap -i -r genome.fa -o index`
- Flags:
**-i**: indexing genome flag 
**-r**: reference genome
**-o**: output name

### 2.1) Reference Genome
- Genome assembly: GRCh38.p14 <br />
&#x1F538; **Only the chromosome 22 human genome** (*only because it is one of the smallest!!*)
- Chromosome 22 info : https://useast.ensembl.org/Homo_sapiens/Location/Chromosome?r=22 <br />

### 2.2) Processing & mapping data 
- Chromap: *~35 sec* <br />
`chromap --preset chip -x index_chrm22 -r Homo_sapiens.GRCh38.dna.chromosome.22.fa -q 20 --min-read-length 10   -1 SRR5063143_naive_H3K27ac_PE.fastq  -o  SRR5063143_naive_H3K27ac_chromap.bed` <br />
- Flags:
**--preset chip**:Mapping chip reads
**-x**:index 
**-r**:reference genome 
**-q**:Min MAPQ (*Quality*)
**--min-read-length**:min length read
**-1**:Single-end (*if use paired-end also included -2*) 
**-o**: output file
*--trim-adapters(not used)*

**check files**: Output file (*BED format*) <br /> 
`head SRR5063143_naive_H3K27ac_chromap.bed`

 &#x1F539; chrm; start; end; N; q; strand. <br />  
  22 &nbsp; 10510250 &nbsp; 10510300 &nbsp; N &nbsp; 59 &nbsp; + <br /> 
  22 &nbsp; 10510252 &nbsp; 10510302 &nbsp; N &nbsp; 46 &nbsp; - <br /> 
  22 &nbsp; 10511600 &nbsp; 10511650 &nbsp; N &nbsp; 60 &nbsp; + <br /> 

  ### 2.3) Pos-mapping data 
*2.3.1)* Convert bed to bam *~2sec* <br /> 
`bedtools bedtobam  -i  SRR5063143_naive_H3K27ac_chromap.bed  -g net/hawkins/vol1/home/aolima/CSHL_Course/genome/chrom22.sizes > SRR5063143_naive_H3K27ac_chromap.bam` <br /> 
&#x1F538; *-g flag*: it is the sizes of each chromossome

&#x1F539; **Extra** save the size for each chromossome  <br />
- **MUST!!** use the same version of reference genome use on the analysis <br />
`samtools faidx genome.fa <br /> <br /> `
`cut -f1,2 genome.fa.fai > sizes.genome <br />` 

**check files**: Output file (*BAM format*) <br /> 
`samtools view SRR5063143_naive_H3K27ac_chromap.bam | head -n 5` 

*2.3.2)* Sort .**bam** & index generation **.bai** & convert to ***.bw** (*BigWig*) *xsec*  <br />
*a.)* `samtools sort SRR5063143_naive_H3K27ac_chromap.bam  -o SRR5063143_naive_H3K27ac_treat.bam` <br />
*b.)* `samtools index SRR5063143_naive_H3K27ac_treat.bam` <br />
*c.)* `bamCoverage -p max -b SRR5063143_naive_H3K27ac_treat.bam  --normalizeUsing RPKM  -v  -o SRR5063143_naive_H3K27ac_norm.bw` <br />
&#x1F538; *a.)* sort the bam files; *b.)* create a index; *c.)* convert the bam to bw & normalize data RPKM (deeptools) <br />

&#x1F539; **Extra** Remove the Chrm MT <br />
- Chromosome MT (Mitocondrial) can cause noise in the *calling peaks* should be remove from the *.bam files  <br />
`samtools index ${sorted.bam.file}`  <br />
`samtools idxstats ${sorted.bam.file} | cut -f1 | grep -v Mt | xargs samtools view -b ${sorted.bam.file}  > ${sorted-noMT.bam.file}  <br />
 &#x1F538; Mt depend the reference genome *(check the reference and annotation genome)*; idxstats index create. <br />
 
#### 3) Peak Calling 
**MACS2** the Model-based Analysis of ChIP-Seq (MACS) for chormatin data analysis https://pypi.org/project/MACS2/ <br />
  










  






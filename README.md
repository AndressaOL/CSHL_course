# Chromatin Workshop  12/4/2023

## CSHL Fall 2023 Computational Genomics Course
Welcome! This is the material and tutorials for the Chromatin Workshop.
The database adopted in this course is under the reference: "Enhancer Chromatin and 3D Genome Architecture Changes from Naive to Primed Human Embryonic Stem Cell States". https://www.sciencedirect.com/science/article/pii/S2213671119301262?via%3Dihub
**Data info:**
- **Cell Type**: Naïve cells
- **Histone modification**: H3K4me3: *Promoters* & H3K27ac: *Promoters and Enhancers*
- **Library info**: **1)** SE-fastq files; **2)** 3 M reads
- **Data File**: Total 6 files (2 for each histone modification & 2 input files)
##

## Tools and Packages Required: <br />
- deepTools:https://deeptools.readthedocs.io/en/develop  <br />
- Samtools: https://www.htslib.org/  <br />
- Chromap: https://github.com/haowenz/chromap <br />
- bedtools: https://bedtools.readthedocs.io/en/latest/  <br />
- MACS2: https://pypi.org/project/MACS2/  <br />
- GREAT: http://great.stanford.edu/public/html/  <br />
- SRplot: https://www.bioinformatics.com.cn/en  <br />

## Install deepTools <br />
**1)** Install miniconda in your home directory <br />
1.1) `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh` *remember to change the permissions* <br />
&#x1F538; check here for others versions: https://docs.conda.io/projects/miniconda/en/latest/ <br />
 1.2) `bash Miniconda3-latest-Linux-x86_64.sh` *install miniconda* <br />
 1.3) activate base: `source /grid/genomicscourse/home/oliveira/miniconda3/bin/activate` *check if is installed conda list* <br />
 **2)** Install deepTools <br />
 2.1) `conda create --name deeptools` *create a new environment* <br />
 2.2) `conda install -c bioconda deeptools` <br /> 
 
## 1) Quality Control of Sequencing using FastQC/MultiQC
- **Documentation**: *FastQC*: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ & *MultiQC*: https://multiqc.info/
- FastQC: ~3 -4 min <br /> 
`fastqc *.fastq`
- MultiQC: ~ 1 min <br />
`multiqc *.zip`

### 2) Processing data & Genome Mapping
**Chromap** for aligning and preprocessing high throughput chromatin profiles (*ATAC-seq & ChIP-seq*); **1)** Trimm the low-quality reads and adaptors; **2)** Remove duplicated reads; **3)** Perform the mapping. https://github.com/haowenz/chromap <br />

- Build the indexed genome *(available!!)* ~ 1 min
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

`head SRR5063143_naive_H3K27ac_chromap.bed` <br />
- **chrm; start; end; N; q; strand**. <br />
22 &nbsp; 10510250 &nbsp; 10510300 &nbsp; N &nbsp; 59 &nbsp; + <br /> 
22 &nbsp; 10510252 &nbsp; 10510302 &nbsp; N &nbsp; 46 &nbsp; - <br /> 
22 &nbsp; 10511600 &nbsp; 10511650 &nbsp; N &nbsp; 60 &nbsp; + <br /> 

### 2.3) Pos-mapping data 
2.3.1) Convert bed to bam *~2sec* <br /> 
`bedtools bedtobam  -i  SRR5063143_naive_H3K27ac_chromap.bed  -g net/hawkins/vol1/home/aolima/CSHL_Course/genome/chrom22.sizes > SRR5063143_naive_H3K27ac_chromap.bam` **-g flag**: sizes for each chromossome; *Extra: save the size for each chromossome*  <br />

&#x1F538;**MUST!!** use the same version of reference genome use on the analysis <br />
```
samtools faidx genome.fa 
cut -f1,2 genome.fa.fai > sizes.genome
``` 
2.3.2) Sort .**bam** & index generation **.bai** & convert to ***.bw** (*Big*Wig*) *~3min*  <br />
*a.)* `samtools index  SRR5063143_naive_H3K27ac_chromap.bam`  <br />
*b.)* `samtools sort SRR5063143_naive_H3K27ac_chromap.bam  -o SRR5063143_naive_H3K27ac_treat.bam` <br />
*c.)* `samtools index SRR5063143_naive_H3K27ac_treat.bam` <br />
*d.)* `bamCoverage -p max -b SRR5063143_naive_H3K27ac_treat.bam  --normalizeUsing RPKM  -v  -o SRR5063143_naive_H3K27ac_norm.bw` **a.)** sort the bam files; **b.)** create a index; **c.)** convert bam to bw & normalize data RPKM (deeptools) <br />

- **Extra** Remove the Chrm MT; &#x1F538; Mt depends the reference genome *(check the reference and annotation genome)*; idxstats: create the index. <br />
Chromosome MT (Mitochondrial)) can cause noise in the *calling peaks* should be removed from the *.bam files  <br />
```
samtools index ${sorted.bam.file} 
samtools idxstats ${sorted.bam.file} | cut -f1 | grep -v Mt | xargs samtools view -b ${sorted.bam.file}  > ${sorted-noMT.bam.file}
```

 **check files**: Output file (*BAM format*); 
 - *Check the biwig files in the genome browser*  <br />
`samtools view SRR5063143_naive_H3K27ac_chromap.bam | head -n 5` <br />
- Use the IGV app: https://igv.org/app/  <br />
- Select the hg38 genome and select the chromosome 22 *(chr22)* <br />
- upload the *.bw files from the HPC to your personal PC *can use SCP or STFP*  <br />
- upload the *.bw in the *track* function  <br />
- Let's have fun!! check the **FBXO7** gene  <br />
  *https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000100225;r=22:32474676-32498829
 
## 3) Peak Calling 
**MACS2** the Model-based Analysis of ChIP-Seq (MACS) for chormatin data analysis https://pypi.org/project/MACS2/ <br />
**Analysis for ChIP-seq; ATAC-seq; Cut&Tag**. *The parameters depend on the data type.*  <br />

**3.1)** Samples info:  <br />
- SRR5063143_naive_H3K27ac_treat.bam *(Performed here!!)*  <br />
- Extra files: SRR5063144_naive_H3K27ac_treat.bam; SRR5063149_naive_H3K4me3_treat.bam; SRR5063150_naive_H3K4me3_treat.bam; SRR5063153_naive_input_treat.bam; SRR5063154_naive_input_treat.bam *(available !!)*  <br />
- Histone modification" *H3K27ac:* broad peaks; *H3K4me3* narrow peaks. <br />   

**3.2)** MACS2 *~2 min* <br />  

3.2.1) H3K4me3 *(Narrow Peaks)*  <br />  
`macs2 callpeak  -t  SRR5063149_naive_H3K4me3_treat.bam -c SRR5063154_naive_input_treat.bam -f BAM  -g hs  --nomodel  -n H3K4me3 --outdir ${your_path_directory}  2> H3K4me3_macs2.log` **flag**: -g:*effective genome size (hs)* <br />  
&#x1F538; MACS2 has effective human genome size, non-model genome uses the effective genome size <br />  
- *hs*:2.7e9
- *mm*:1.87e9
- *ce*:9e7
- *dm*:1.2e8
**check files**: MACS2 generates several outputs; only check *.log and *.narrowPeaks <br />  
To check the output narrowPeaks file uses: `wc-l` to count the number of peaks and `head` & `ls -ll` to check the output file <br /> 

3.2.2) H3K27ac *(Broad Peaks)* <br /> 
`macs2 callpeak  -t  SRR5063143_naive_H3K27ac_treat.bam -c SRR5063153_naive_input_treat.bam -f BAM  -g hs -n H3K27ac  --nomodel  --broad --outdir ${your_path_directory} 2> H3K27ac_broad_macs2.log` <br /> 

##QC Analysis *~ 3 min* 
**Fraction of reads in peaks (FRiP):** FRiP Score essential to evaluate the Peaks Quality. *more details:* https://yiweiniu.github.io/blog/2019/03/Calculate-FRiP-score/ <br />
- Request data: *.bam files & *Peaks files (Narrow or broad)
- To calcualte the FRiPs run the script *(Pipelines folder)*
```
bash FRiP_score.sh  SRR5063143_naive_H3K27ac_treat.bam  H3K27ac_peaks.broadPeak
bash FRiP_score.sh SRR5063149_naive_H3K4me3_treat.bam H3K4me3_peaks.narrowPeak
``` 

## Let's visualize our results with Plots
**1)** Correlation plots and Heatmap. Correlation matrix bewteen the replicates (*QC analysis*) and Heatmap (*visualize the signal intensity:Input; HK3me4;HK27ac*) 






















  






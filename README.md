# Cozzolino2023
Scripts used for data processing and analysis in manuscript about the effect of Mediator kinase inhibition on interferon signaling in Down syndrome. With the exception of the splicing analysis, all RNA-seq based analysis was performed by Kira Cozzolino with guidance from Lynn Sanford, Taylor Jones, Mary Allen, and Sam Hunter in the Dowell and Allen labs. Splicing analysis was performed by Benjamin Erickson in the Bentley lab. Separately, ANOVA analysis of cytokine screen results was performed by Taylor Jones using the deposited scripts. 

RNA-SEQ DATA
Data Processing
- Indexes: Indexes were constructed using hisat2 version 2.1.0 using the hisat2-build command on hg38. 
- Quality control: Data was received as fastq files. Quality of all files was checked using fastqc version 0.11.5
- Trimming: Reads were trimmed using bbduk, specifically bbmap version 38.05. The command and parameters for bbduk were as follows:

for fastq in /path/to/fastq/files/*R1.fastq.gz; do
sample=$(basename ${fastq} _R1.fastq.gz)

bbduk.sh -Xmx40g \
    t=32 \
    in=${fastq} \
    in2=/path/to/fastq/files/${sample}_R2.fastq.gz \
    out=/path/to/output/files/${sample}_R1.trim.fastq.gz \
    out2=/path/to/output/files/${sample}_R2.trim.fastq.gz \
    ref=${bbmap_adapters} \
    ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
    maq=10 minlen=25 \
    tpe tbo \
    literal=AAAAAAAAAAAAAAAAAAAAAAA \
    stats=${sample}.trimstats.txt \
    refstats=${sample}.refstats.txt

- Mapping: Reads were mapped using hisat2 version 2.1.0. The hisat2 command was run first with the flags -p 32, using hg38 as a reference genome. Samtools version 1.10 was then used to convert .sam files to .bam and .flagstat files. This was done by first running the samtools view command with the flags  -@ 32 -bS -o. The samtools sort command was next run on the files using the flags -@ 32. The samtools index command was used with the flags -@ 32, and then the flagstat command was used.
- Counting: Featurecounts was used to count mapped reads for each sample.  Featurecounts R script with parameters used in this study has been uploaded. Necessary libraries: Rsubread.
- Quantification and statistics was performed using deseq2. Scripts used for both within-genotype and cross-genotype (containing normalizaton for ploidy of genes on chromsome 21) have been uploaded. Visualization of data was also performed in R. Necessary libraries: DeSeq2, dplyr, textshape, tidyverse, ggplot2, data.table, plyr, RCp;prBrewer, circlize, and ComplexHeatmap. 

- Splicing analysis: Raw fastq files were processed separately by Benjamin Erickson and analyzed to detect splicing changes. The workflow was as follows: duplicates were removed and adaptors trimmed using bbTools (v39.01), trimmed reads were mapped against hg38 using hisat2 (2.1.0), and mapped reads were processed with rMATS (4.0.1) using the deposited script. rMATS results were filtered based on inclusion cutoffs of FDR<0.05, absolute(IncLevelDifference)<0.2, and >=2 reads/replicate. Sashimi plots were generated from rMATS results using a modified script based on ggashimi.py

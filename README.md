# Cozzolino2023
Scripts used for data processing and analysis in manuscript about the effect of Mediator kinase inhibition on interferon signaling in Down syndrome

Data Processing
- Quality control: Data was received as fastq files. Quality of all files was checked using fastqc version 0.11.5
- Trimming: Reads were trimmed using bbduk, specifically bbmap version 38.05. The command and parameters for bbduk were as follows:

for fastq in /path/to/fastq/files/*R1.fastq.gz; do
sample=$(basename ${fastq} _R1.fastq.gz)

bbduk.sh -Xmx40g \
    t=32 \
    in=${fastq} \
    in2=/path/to/fascs{sample}_R2.fastq.gz \
    out=/scratch/Users/kico4293/RNAseq_April2022/rep3/fastq/${sample}_R1.trim.fastq.gz \
    out2=/scratch/Users/kico4293/RNAseq_April2022/rep3/fastq/${sample}_R2.trim.fastq.gz \
    ref=${bbmap_adapters} \
    ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
    maq=10 minlen=25 \
    tpe tbo \
    literal=AAAAAAAAAAAAAAAAAAAAAAA \
    stats=${sample}.trimstats.txt \
    refstats=${sample}.refstats.txt

- 

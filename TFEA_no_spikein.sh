##################### LOAD MODULES ######################################
#source ~/venv_TFEA_tj/bin/activate

export R_LIBS_USER=/Users/kico4293/R/x86_64-pc-linux-gnu-library/4.2

module purge
module load python/3.6.3
module load python/3.6.3/matplotlib/1.5.1
module load python/3.6.3/scipy/0.17.1
module load python/3.6.3/numpy/1.14.1
module load python/3.6.3/htseq/0.9.1
module load python/3.6.3/pybedtools/0.7.10

module load samtools/1.8
module load bedtools/2.25.0
module load meme/5.0.3

###################### SET PATHS #######################################

GENOME='/scratch/Shares/dowell/genomes/hg38/hg38.fa'
MOTIF='/scratch/Shares/dowell/motifs/HOCOMOCODatabaseFIMO/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme'
#MOTIF='/scratch/Users/tajo5912/dbapplication/scripts/best_curated_Human.meme'

BAMS='/scratch/Users/kico4293/PROseq_May2023/NascentFlow_hg38/mapped/bams'
BED='/scratch/Users/kico4293/PROseq_May2023/PROseq_MUMERGE.bed'
OUTDIR='/scratch/Users/kico4293/PROseq_May2023/TFEA_Output'

ROOTNAME1=$1
ROOTNAME2=$2

########################################################################

mkdir -p ${OUTDIR}/${ROOTNAME1}_vs_${ROOTNAME2}

TFEA --sbatch kico4293@colorado.edu --cpus 16 \
--output ${OUTDIR}/${ROOTNAME1}_vs_${ROOTNAME2} \
--combined_file ${BED} \
--bam1 ${BAMS}/${ROOTNAME1}_1.sorted.bam ${BAMS}/${ROOTNAME1}_2.sorted.bam \
--bam2 ${BAMS}/${ROOTNAME2}_1.sorted.bam ${BAMS}/${ROOTNAME2}_2.sorted.bam \
--label1 ${ROOTNAME1} --label2 ${ROOTNAME2} \
--genomefasta ${GENOME} \
--fimo_motifs ${MOTIF} \
--debug


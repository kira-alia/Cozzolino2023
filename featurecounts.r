library("Rsubread")

# featurecounts_files.txt is a tab-delimited file that contains 
# all of the filenames of the mapped BAM files to count
# and a second column with the sample names
filelist_input <- "/directions/to/featurecounts_files.txt"
gtf_input <- "/directions/to/hg38_refseq_genenames_included.gtf"
outfile <- "/directions/to/outfile/for/output.txt"

inputData <- read.table(filelist_input, sep="\t")

# Put file and sample lists into variable
fileList <- as.character(t(inputData["V1"]))
sampleList <- as.character(t(inputData["V2"]))

# Load in GTF file for metadata and get unique gene/transcript pairs
gtf_table <- read.table(gtf_input)
gtf_rows <- unique(gtf_table[c("V1","V10","V13")])
for (rownum in 1:length(gtf_rows$V10)){gtf_rows$V1[rownum] = gtf_table$V1[match(gtf_rows$V13[rownum],gtf_table$V13)]}

# Run featureCounts on these BAM files - look up featureCounts manual to understand parameters
coverage <- featureCounts(files=fileList,
                    annot.ext=gtf_input,
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="transcript_id",
                    useMetaFeatures=TRUE,
                    allowMultiOverlap=TRUE,
                    largestOverlap=TRUE,
                    countMultiMappingReads=FALSE,
                    isPairedEnd=TRUE,
                    strandSpecific=1,
                    nthreads=8)

# Write out annotation and count data as tab delimited txt files
# Count data is just counts
# Annotation data has the gene length information for calculating TPM

colnames(coverage$counts) <- sampleList
coverage$annotation["TranscriptID"] <- coverage$annotation["GeneID"]
coverage$annotation["GeneID"] <- gtf_rows["V10"]

write.table(x=data.frame(coverage$annotation[,c("GeneID","TranscriptID","Length")],
                         gtf_rows$V1,coverage$counts,stringsAsFactors=False),
            file=outfile,row.names=FALSE,sep='\t',quote=FALSE)

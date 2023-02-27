library(DESeq2)
library(dplyr)
library(textshape)
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(plyr)

#directions to the metadata
cutoff=0.1
metadatafile="/path/to/metadate/file"
coveragetablefile="/path/to/featurecounts/output"
outdir="/path/to/directory/"

#directions to the annotation file, pull out T21 genes
masterannotationdf=read.table("/path/to/directory/containing/hg38_refseq_genenames_included.gtf",
                              sep="\t", header = TRUE)
masterannotationdf
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
masterannotationdf_only21 <- separate(masterannotationdf_only21, col = gene_id,into = c("gene_string","name","transcript","transcript_id"),sep=" ")
masterannotationdf_only21$name <- gsub(x = masterannotationdf_only21$name,";","")
unique(masterannotationdf_only21$name)
aneuploidygenes <- unique(masterannotationdf_only21[["name"]]) #Only chr21 genes

#read the metadata, set controls
metadata <- read.csv(metadatafile, header=TRUE, sep=",", stringsAsFactors = TRUE)
metadata$ploidy <- relevel(metadata$ploidy , "ploidy_typical")
head(metadata)
metadata$treatmentIFN <- relevel(metadata$treatmentIFN , "CTRL")
metadata$treatmentCA <- relevel(metadata$treatmentCA , "CTRL")
head(metadata)

#read and check the coverage table file, dimensions should be sample number + 2
coveragetable <- read.csv(coveragetablefile,sep="\t",row.names = 2)
head(coveragetable)
countdat <- coveragetable
dim(countdat)
head(countdat)

baseploidy <- 2
alt_ploidy <-3

#filter for most highly expressed isoform
countdat$sums <- rowSums(countdat[,c(3:ncol(countdat))])/countdat$Length
countdat <- (countdat[order(-countdat$sums,countdat$GeneID),])
countdat <- countdat[!(duplicated(countdat$GeneID)),]
rownames(countdat) <- countdat$GeneID
countdat <- countdat[,-c(1,2,ncol(countdat))]
head(countdat)

#T21 normalization
controlgenes <-ifelse(rownames(ddsboth) %in% aneuploidygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
ploidy_trisomy21 = ifelse(rownames(ddsboth) %in% aneuploidygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(ddsboth))

normFactors <- matrix(do.call(cbind, mget(paste0(ddsboth$ploidy))),
                      ncol=ncol(ddsboth),nrow=nrow(ddsboth),
                      dimnames=list(1:nrow(ddsboth),
                                    1:ncol(ddsboth)))
head(countdat)
View(normFactors)
normFactors <- normFactors/baseploidy
normFactors
unique(normFactors)

#re-running with T21 normalization factors
model.matrix(~treatmentIFN*treatmentCA*ploidy,data=metadata)
ddsboth <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~treatmentIFN*treatmentCA*ploidy)
ddsboth <- DESeq(ddsboth)
DEddsboth <- ddsboth
resultsNames(ddsboth)

ddsCollapsed_normfactor <-estimateSizeFactors(ddsboth,normMatrix=normFactors, 
                                              controlGenes=controlgenes) 
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) 
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) 
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor)  
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) 

#view results - below is results output for T21 vs D21 comparison. Other results names include the ones for comparing differential IFN- or CA-responsiveness in T21 vs D21
resboth=results(ddsCollapsed_normfactor, name="ploidy_ploidy_trisomy21_vs_ploidy_typical")
summary(resboth)
head(resboth)
resdata <- as.data.frame(resboth)
resdata$chr <- ifelse(rownames(resdata) %in% aneuploidygenes,"chr21","not_chr21")
View(resdata)
results_sig <- subset(resboth, padj < 0.05)
results_sig$GeneID <- rownames(results_sig)
results_extrasig <- subset(results_sig, abs(log2FoldChange)>1)
results_extrasig$GeneID <- rownames(results_extrasig)
#view(results_sig)

# Volcano Plot
ggplot(resdata, aes(log2FoldChange, -log10(padj))) + 
  geom_point() +
  theme_classic() + ##use the classic theme template 
  xlab("log2 Fold Change") + ##label the x-axis
  ylab("-log10 Adjusted P-value") + ##label the y-axis
  ggtitle("Sample Title") + ##to add title
  theme(plot.title = element_text(hjust = 0.5)) + ##to center title
  geom_point(data=as.data.frame(results_extrasig),aes(log2FoldChange, -log10(padj)), color="red",size=1.25, alpha=1) +
  xlim(-12,10)

#basic MA plot of data and export as data table
DESeq2::plotMA(resboth,  main="resboth")
head(resboth[ order( resboth$padj ), ])
outfilename=paste(outdir, "filename.bed", sep="")
write.table(resboth, file = outfilename, sep="\t",  quote = FALSE, row.names = TRUE, col.names = TRUE)

#Generate ranked dataframe with gene names in first column, second column is log-p, value multiplied by sign of fold change
rnkdf <- tibble(gene = rownames(resdata),
                rnk = -log10(resdata$padj) * sign(resdata$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()

#Write out table without additional info for GSEA
outfilename=paste(outdir, "samplename.rnk", sep="")
write.table(rnkdf, file = outfilename, 
            append = FALSE, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")

#Prettier MA plot
ggplot(resdata, aes(log(baseMean,2), log2FoldChange)) + 
  geom_point(color="gray20", alpha = 0.5, size = 1) + ##change opacity of points
  theme_bw() + ##set the theme of plot 
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle("Differentially CA-responsive Genes - T21 vs D21") + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18,face = "bold"),
        axis.text = element_text(size = 16)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=as.data.frame(results_sig),aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75) +
ggsave("MA_Plot_T21vsD21_CA.png", plot = last_plot(), path = outdir,
       scale = 1, width = NA, height = NA, units = c("mm"),
       dpi = 300, limitsize = TRUE)

#violin plot of Log2FC on chr21
lesschrs <- c("chr21")
medresdata <- resdata[!(is.na(resdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]

ggplot() +
  geom_violin(data=resdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  theme_classic(base_size = 30) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1,
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1,
              position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=resdata,aes(x=chr, y=log2FoldChange),fun.y=median,
               geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

#View histogram of p-values (for FDR visualization)
hist(as.numeric(resboth$pvalue), breaks=50, col="grey")
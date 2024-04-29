library(DESeq2)
library(dplyr)
library(textshape)
library(tidyverse)
library(ggplot2)
library(data.table)

#directions to the metadata
cutoff=0.1
metadatafile="directions/to/metadata/file"
coveragetablefile="directions/to/featurecounts/output"
outdir="directions/to/out/directory"

#read the metadata
metadata <- read.csv(metadatafile, header=TRUE, sep=",", stringsAsFactors = TRUE)
head(metadata)
metadata$genotype <- relevel(metadata$ploidy , "ploidy_typical")
metadata$treatmentIFN <- relevel(metadata$treatmentIFN , "CTRL")
metadata$treatmentCA <- relevel(metadata$treatmentCA , "CTRL")
metadata$group <- relevel(metadata$group , "groupA") #this will change depending on samples being compared
head(metadata)

#read and check the coverage table file, dimensions should be sample number + 2
coveragetable <- read.csv(coveragetablefile,sep="\t",row.names = 2)
head(coveragetable)
countdat <- coveragetable
dim(countdat)
head(countdat)

#assign variables for comparison - sample2 is control sample in each comparison
sample2="groupA"
sample1="groupB"
contrastcol = "group"

#filter for most highly expressed isoform of each gene
countdat$sums <- rowSums(countdat[,c(3:ncol(countdat))])/countdat$Length
countdat <- (countdat[order(-countdat$sums,countdat$GeneID),])
countdat <- countdat[!(duplicated(countdat$GeneID)),]
head(countdat)
rownames(countdat) <- countdat$GeneID
countdat <- countdat[,-c(0,1,2,ncol(countdat))]
head(countdat)

# Batch correct - this won't be needed for all comparisons
cov1 <- factor(metadata$group)
cov1

sva_corrected <- ComBat_seq(counts = as.matrix(countdat),
                            batch=metadata$batch,
                            group=cov1)
sva_corrected

metadata$batch
write.csv(sva_corrected, file=paste(outdir,"filename_corrected.csv",sep=""))
countdata_corrected <- as.data.frame(sva_corrected)
countdata_corrected <- round(countdata_corrected,digits=0)
write.csv(countdata_corrected, file=paste(outdir,"sva_corrected_rounded.csv",sep=""))

#run Deseq on RNA-seq
ddsboth <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~group) #use this for all samples; if you're sva correcting use countdata_corrected instead
ddsboth <- DESeq(ddsboth)
DEddsboth <- ddsboth

#view result statistics, designate subset of genes for labeling on volcano plot
resboth=results(DEddsboth, contrast=c(contrastcol,sample1,sample2))
summary(resboth)
head(resboth)
resdata <- as.data.frame(resboth)
results_sig <- subset(resboth, padj < 0.05)
results_sig$GeneID <- rownames(results_sig)
results_extrasig <- subset(results_sig, abs(log2FoldChange)>1)
results_extrasig$GeneID <- rownames(results_extrasig)
sizeFactors(ddsboth)

# Volcano Plot
ggplot(resdata, aes(log2FoldChange, -log10(padj))) + 
  geom_point() +
  theme_classic() + ##set the theme 
  xlab("log2 Fold Change") + ##label the x-axis
  ylab("-log10 Adjusted P-value") + ##label the y-axis
  ggtitle("Sample title") + ##to add title
  theme(plot.title = element_text(hjust = 0.5)) + ##to center title
  geom_point(data=as.data.frame(results_extrasig),aes(log2FoldChange, -log10(padj)), color="red",size=1.5, alpha=1) +
  xlim(-5,5) 

#MA plot
ggplot(resdata, aes(log(baseMean,2), log2FoldChange)) + 
  geom_point(color="gray20", alpha = 0.5, size = 1) + ##change opacity of points
  theme_bw() + ##set the theme of plot 
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle("Sample Title") + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18,face = "bold"),
        axis.text = element_text(size = 16)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=as.data.frame(results_sig),aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75) +
  ggsave("Sample Title", plot = last_plot(), path = outdir,
         scale = 1, width = NA, height = NA, units = c("mm"),
         dpi = 300, limitsize = TRUE)

#basic MA plot of data and export as data table
DESeq2::plotMA(resboth,  main="resboth")
head(resboth[ order( resboth$padj ), ])
outfilename=paste(outdir, "filename.bed", sep="")
write.table(resboth, file = outfilename, sep="\t",  quote = FALSE, row.names = TRUE, col.names = TRUE)

#Generate ranked dataframe with gene names in first column, second column is log-p, value multiplied by sign of fold change
rnkdf <- tibble(gene = rownames(resdata),
                rnk = -log10(resdata$padj) * sign(resdata$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% na.omit()

#Write out table without additional info for GSEA
outfilename=paste(outdir, "samplename.rnk", sep="")
write.table(rnkdf, file = outfilename, 
            append = FALSE, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t")

library(tidyverse)
library(ggplot2)
library(data.table)

#directions to the metadata
#coveragetablefile="/directions/to/fc_counts_human.txt"
outdir="/directions/to/outdir"
ifnfile="/directions/to/inf_up.csv"
normalcounts="/directions/to/normalizedcounts.csv"

normcounts <- read.csv(normalcounts)
head(normcounts)

ifngenes <- read.csv(ifnfile)
head(ifngenes)
# #read and check the coverage table file, dimensions should be sample number + 2
# coveragetable <- read.csv(coveragetablefile,sep="\t",row.names = 2)
# 
# #### pull in norm counts and lens
# #normcounts <- read.csv(paste0(outdir,'normalized_nmID.csv'), header=TRUE, row.names = 'X')
# 
# countdat <- coveragetable
# countdat$sums <- rowSums(countdat[,c(3:ncol(countdat))])/countdat$Length
# countdat <- (countdat[order(-countdat$sums,countdat$GeneID),])
# countdat <- countdat[!(duplicated(countdat$GeneID)),]
# head(countdat)
# 
# countdatlen <- countdat %>% select(c('GeneID', 'Length'))
# head(countdatlen)
# dim(countdatlen)
# dim(normcounts)
# 
# write.csv(as.data.frame(countdatlen), file=paste(outdir,"normalized_nmID_len.csv",sep=""))
# write.csv(as.data.frame(normcounts), file=paste(outdir,"normalized_nmID.csv",sep=""))

#head(normcounts)
#countsandlen <- (merge(normcounts, countdatlen, by.x=0, by.y=0))
#head(countsandlen)
#dim(countsandlen)

#countsandlen$D21_VV_1_rpkm <- countsandlen$D21_VV_1/countsandlen$Length
# countsandlen$D21_VV_2_rpkm <- countsandlen$D21_VV_2/countsandlen$Length
# countsandlen$D21_VV_3_rpkm <- countsandlen$D21_VV_3/countsandlen$Length
# countsandlen$D21_VIFN_1_rpkm <- countsandlen$D21_VIFN_1/countsandlen$Length
# countsandlen$D21_VIFN_2_rpkm <- countsandlen$D21_VIFN_2/countsandlen$Length
# countsandlen$D21_VIFN_3_rpkm <- countsandlen$D21_VIFN_3/countsandlen$Length
# 
# countsandlen$T21_VV_1_rpkm <- countsandlen$T21_VV_1/countsandlen$Length
# countsandlen$T21_VV_2_rpkm <- countsandlen$T21_VV_2/countsandlen$Length
# countsandlen$T21_VV_3_rpkm <- countsandlen$T21_VV_3/countsandlen$Length
# countsandlen$T21_VIFN_1_rpkm <- countsandlen$T21_VIFN_1/countsandlen$Length
# countsandlen$T21_VIFN_2_rpkm <- countsandlen$T21_VIFN_2/countsandlen$Length
# countsandlen$T21_VIFN_3_rpkm <- countsandlen$T21_VIFN_3/countsandlen$Length

# countsandlen$D21_VV_1_logRPKM <- log(countsandlen$D21_VV_1_rpkm, 2)
# countsandlen$D21_VIFN_1_logRPKM <- log(countsandlen$D21_VIFN_1_rpkm, 2)
# countsandlen$T21_VV_1_logRPKM <- log(countsandlen$T21_VV_1_rpkm, 2)
# countsandlen$T21_VIFN_1_logRPKM <- log(countsandlen$T21_VIFN_1_rpkm, 2)

plotdf <- filter(normcounts, T21_VV_avg >0 & D21_VV_avg >0 & D21_VIFN_avg >0 & T21_VIFN_avg>0)
head(plotdf)

ifnvec <- as.vector(ifngenes$GeneID)
ifnvec

names(plotdf)
min(plotdf$D21_VV_avg)
min(plotdf$T21_VV_avg)

# create multiple linear model
plotdfIFN <- filter(plotdf, GeneID %in% ifnvec)
head(plotdfIFN)
IFNlm_fit <- lm(T21_VV_avg ~ 0  + D21_VV_avg, data=plotdfIFN)
summary(IFNlm_fit)
IFNcf <- coef(IFNlm_fit)
InterceptIFN <- 0
SlopeIFN <- IFNcf[1]

IFNlm_fit_TX <- lm(T21_VIFN_avg ~ 0  + D21_VIFN_avg, data=plotdfIFN)
summary(IFNlm_fit_TX)
IFNcf_TX <- coef(IFNlm_fit_TX)
InterceptIFN_TX <- 0
SlopeIFN_TX <- IFNcf_TX[1]

lm_fit <- lm(T21_VV_avg ~ 0 + D21_VV_avg, data=plotdf)
summary(lm_fit)
cf <- coef(lm_fit)
Intercept <- 0
Slope <- cf[1]
Slope

lm_fittx <- lm(T21_VIFN_avg ~ 0 + D21_VIFN_avg, data=plotdf)
summary(lm_fittx)
cftx <- coef(lm_fittx)
Intercepttx <- 0
Slopetx <- cftx[1]
Slopetx

#plot everything
ggplot(plotdf, aes(D21_VV_avg, T21_VV_avg)) +
  geom_point(color='grey', alpha = 0.5) + ##change opacity of points
  geom_point(data=plotdf,aes(D21_VIFN_avg, T21_VIFN_avg), color='grey',size=1, alpha=0.75) +  
  theme_classic() + ##use the classic theme template
  xlab('log2 Average Normalized Counts') + ##label the x-axis
  ylab('log2 Fold Change') + ##label the y-axis
  ggtitle(paste0('MA plot')) + ##to add title
  scale_x_log10() + scale_y_log10() +
  geom_point(data=plotdfIFN,aes(D21_VV_avg, T21_VV_avg), color='red',size=1, alpha=0.75) +
  geom_point(data=plotdfIFN,aes(D21_VIFN_avg, T21_VIFN_avg), color='green',size=1, alpha=0.9) +
  geom_abline(
    slope=Slopetx,
    intercept=Intercepttx) +
  geom_abline(
    slope=Slope,
    intercept=Intercept) +
  geom_abline(
    slope=SlopeIFN,
    intercept=InterceptIFN, color='blue') +
  geom_abline(
    slope=SlopeIFN_TX,
    intercept=InterceptIFN_TX, color='green')

#just vehicle
ggplot(plotdf, aes(D21_VV_avg, T21_VV_avg)) +
  geom_point(color='grey', alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template
  xlab('D21 Vehicle (RPKM)') + ##label the x-axis
  ylab('T21 Vehicle (RPKM)') + ##label the y-axis
  ggtitle(paste0('MA plot')) + ##to add title
  scale_x_log10() + scale_y_log10() +
  geom_point(data=plotdfIFN,aes(D21_VV_avg, T21_VV_avg), color='red',size=1, alpha=0.75) +
  geom_abline(
    slope=Slope,
    intercept=Intercept, color="black") +
  geom_abline(
    slope=SlopeIFN,
    intercept=InterceptIFN, color='red')

#just treatment
ggplot(plotdf, aes(D21_VIFN_avg, T21_VIFN_avg)) +
  geom_point(color='grey', alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template
  xlab('D21 IFN (RPKM)') + ##label the x-axis
  ylab('T21 IFN (RPKM)') + ##label the y-axis
  ggtitle(paste0('MA plot')) + ##to add title
  scale_x_log10() + scale_y_log10() +
  geom_point(data=plotdfIFN,aes(D21_VIFN_avg, T21_VIFN_avg), color='red',size=1, alpha=0.9) +
  geom_abline(
    slope=Slopetx,
    intercept=Intercepttx, color="black") +
  geom_abline(
    slope=SlopeIFN_TX,
    intercept=InterceptIFN_TX, color='red')

genevector <- c('IRF1', 'CXCL9')
plotdfIFN$diff <- plotdfIFN$D21_VIFN_avg - plotdfIFN$D21_VV_avg
qset<- quantile(plotdfIFN$diff, 0.9)
dim(plotdfIFN)
plotdfIFNsmall <- filter(plotdfIFN, diff > qset)
dim(plotdfIFNsmall)
#plotdfIFNsmall <- filter(plotdfIFN, D21_VV_avg > 10, T21_VV_avg > 10)
#plotdfIFNsmall <- filter(plotdfIFNsmall, D21_VV_avg < 100, T21_VV_avg < 100)
plotdfIFNsmall
plotdfIFNsmall <- filter(plotdfIFN, 
                    GeneID %in% genevector)

plotdfIFNsmall$col1 <- 1
plotdfIFNsmall$col2 <- 2

ggplot() +   
  geom_segment(aes(x = plotdfIFNsmall$col1,y = plotdfIFNsmall$D21_VV_avg, 
      xend = plotdfIFNsmall$col2,yend = plotdfIFNsmall$D21_VIFN_avg),
  size = 0.3, alpha=0.2, color='blue') +
  geom_segment(aes(x = plotdfIFNsmall$col1,y = plotdfIFNsmall$T21_VV_avg, 
                   xend = plotdfIFNsmall$col2,yend = plotdfIFNsmall$T21_VIFN_avg),
               size = 0.3, alpha=0.2, color='red') +
  geom_point(data=plotdfIFNsmall,aes(col1, T21_VV_avg), color='red',size=2, alpha=0.5) +
  geom_point(data=plotdfIFNsmall, aes(col1, D21_VV_avg), color='blue', size=2, alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template
  xlab('') + ##label the x-axis
  ylab('') + ##label the y-axis
  ggtitle(paste0('Gene Expression')) + ##to add title
  geom_point(data=plotdfIFNsmall,aes(col2, D21_VIFN_avg), color='blue',size=2, alpha=0.9) +
  geom_point(data=plotdfIFNsmall,aes(col2, T21_VIFN_avg), color='red',size=2, alpha=0.9)


test %>%
  ggplot() +
  geom_point(aes(x=x,y=y),color='red') +
  geom_point(aes(x=vec.x,y=vec.y),color='blue') +
  geom_segment(
    aes(x = x,y = y, xend = vec.x,yend = vec.y),
    arrow = arrow(length = unit(0.03,units = "npc")),
    size = 1
  )


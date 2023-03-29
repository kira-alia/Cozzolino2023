library(argparse)
library(tidyverse)

#parser <- ArgumentParser()
#parser$add_argument("-c", "--cyloop", action="store", dest="cyloop",
                   # help="LOOP")
#parser$add_argument("-r", "--relevel", action="store", dest="relevel",
                    #help="relevel")

#args <- parser$parse_args()
#cyloop <- args$cyloop
#relevel <- args$relevel

#outdir<-paste0('/Users/tajo5912/anova_down/',relevel,'/')
#setwd(outdir)

df <- read.csv(file='C:/Users/kirac/Desktop/TaatjesLab/CytokineScreen/anova_down/av_values_for_anova.txt',
               sep='\t')
head(df)
cydf <- filter(df, cytokine == cyloop)
cydf$sample <- factor(cydf$sample)
cydf$sample <- relevel(cydf$sample , relevel)

tapply(cydf$value, cydf$sample, mean) #for each factor use the count and mean function for each group
# #could also find variance
# tapply(cydf$value, cydf$sample, var) #this is the within group variance
# 
boxplot(cydf$value ~ cydf$sample)
#one way anova, does not assume equal variances and does a welsh correction
oneway.test(cydf$value~cydf$sample)
#another way that gets you way more info
aov.out = aov(value~sample, data=cydf) #aov runs the lm command
sink(file=paste0(outdir,cyloop,'.txt'))
summary.lm(aov.out) #this more resembles what we have seen before
sink()







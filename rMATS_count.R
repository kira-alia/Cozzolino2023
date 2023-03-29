# benjamin.erickson@cuanschutz.edu

# Counting total numbers from rMATS output

library(tidyverse)

sigcounts <- function(comp = "rMATS_folder", type = "SE", min_count = 2, max_FDR=0.05, max_diffrence = 0.2){
  psis <- read_tsv(paste0(comp,'/',type,'.MATS.JC.txt'),show_col_types = FALSE) 
  
  out_length <- length(str_split(psis[1,"IJC_SAMPLE_1"],",")[[1]])
  
  sumofcol <- function(df, col1,col2, colnum) {
    mutate(df, !!{{colnum}} := !! {{col1}} + !! {{col2}})
  }
  
  psis <- psis %>%
    #Split the replicate read counts that are separated by commas into different columns
    separate(., col = IJC_SAMPLE_1, into = paste0('IJC_S1R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_1, into = paste0('SJC_S1R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = IJC_SAMPLE_2, into = paste0('IJC_S2R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_2, into = paste0('SJC_S2R', 1:out_length), sep = ',', remove = T, convert = T)
  
  # Now sum to get reads in each condition for each event and filter.
  psis.filtered <- psis
  for(i in 1:out_length){
    psis.filtered <- sumofcol(psis.filtered,
                              as.name(paste0("IJC_S1R",i)),as.name(paste0("SJC_S1R",i)),
                              as.name(paste0("S1R",i,"counts"))) %>% 
      filter(.,!!as.name(paste0("S1R",i,"counts")) >= min_count) %>% 
      sumofcol(.,
               as.name(paste0("IJC_S2R",i)),as.name(paste0("SJC_S2R",i)),
               as.name(paste0("S2R",i,"counts"))) %>% 
      filter(.,!!as.name(paste0("S2R",i,"counts")) >= min_count)
    
  }
  
  # Defining sensitive exons #only those whose PSI decreases < max_diffrence
  psis.sensitive1 <- filter(psis.filtered, FDR < max_FDR,  IncLevelDifference < -max_diffrence) 
  # Defining sensitive exons #only those whose PSI increase > max_diffrence
  psis.sensitive2 <- filter(psis.filtered, FDR < max_FDR,  IncLevelDifference > max_diffrence) 
  
  # IncLevelDifference > 0, explanation for the different PSI types :
  # SE less exon skiping in condition 1
  # MXE:
  # + include exon1 skip exon 2
  # - include exon2 skip exon 1
  # A5SS more downstream
  # A3SS more upstream
  # RI retains intron
  
  psis.sensitive1 %>%  write_tsv(.,paste0(comp,'/',type,'.MATS.JC.filter_less.txt'),col_names = T)
  psis.sensitive2 %>%  write_tsv(.,paste0(comp,'/',type,'.MATS.JC.filter_more.txt'),col_names = T)
  
  v1 = paste("less_diff",type,comp,sep="_")
  v2 = paste("greater_diff",type,comp,sep="_")
  tibble(!!v1 := nrow(psis.sensitive1), !!v2 :=nrow(psis.sensitive2))
  
}

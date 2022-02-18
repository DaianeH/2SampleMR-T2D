library(metafor)
library(plyr) 
library(meta) 
library(rmeta) 
library(TwoSampleMR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library("scales")
library(tidyr)
library(stringr)


#Command line usage: Rscript 2mr.r [trait-T2D] [disease-secondary outcome] [disease_dataset-summary statistics]


args <- commandArgs(trailingOnly = TRUE)

df_all <- data.frame(stringsAsFactors = F)
df1 <- data.frame(stringsAsFactors = F)

trait <- as.character(args[1])
disease <- as.character(args[2])

#summary stats of T2D
exp_sT2D <- read_exposure_data(
   filename = "/path/to/MVP.T2D.EUR.dbGaP.txt",
   sep = "\t",
   snp_col = "SNP",
   beta_col = "BETA",
   se_col = "SE",
   effect_allele_col = "EA",
   other_allele_col = "NEA",
   eaf_col = "EAF",
   pval_col = "P"
)

#Read all variants that are part of tissue-grouped SNP sets
T2D_variants <- read.table("path/to/T2D_snps_all.txt", header = T)
T2D_variants <- unique(T2D_variants)

#Select only variants that are part of tissue-grouped SNP sets from summary stats
exp_T2D_sub <- merge(x=exp_T2D,y=T2D_variants,by="SNP",all.y=TRUE)

#summary stats of secondary outcome
outcome<- read_outcome_data(
  snps = exp_T2D_sub$SNP,
  filename = args[3], 
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pval"
)

outcome_sub <- merge(x=outcome,y=T2D_variants,by="SNP",all.y=TRUE)

#Harmonise
full_harm_outcome <- harmonise_data(
  exposure_dat = exp_T2D_sub, 
  outcome_dat = outcome_sub
)

#Read path where all tissue-group SNP sets are
files_SNPs <- list.files(path = paste("path/to/files_tissues/"), pattern = ".txt")

#Run 2 sample MR in each tissue-group SNP set
for (i in 1:(length(files_SNPs))){
   SNPs_tissue <- read.table(file=paste("path/to/files_tissues/",files_SNPs[i], sep=""), stringsAsFactors = F, header = T)
  T2DSNPs_outcome <- merge(x=full_harm_outcome, y=SNPs_tissue, by="SNP", all.y=TRUE)
  T2DSNPs_outcome <- T2DSNPs_outcome[!is.na(T2DSNPs_outcome$effect_allele.exposure), ]
  result1 <- mr(T2DSNPs_outcome, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_weighted_mode"))
  result1 <- result1[result1$method == 'Inverse variance weighted',]
  result1$outcome <- files_SNPs[i]
  result1$b <- result1$b
  result1$se <- result1$se
  result1$odds_ratio <- exp(as.numeric(result1$b))
  result1$ci_lower <- exp(as.numeric(result1$b)-(1.96*as.numeric(result1$se)))
  result1$ci_upper <- exp(as.numeric(result1$b)+(1.96*as.numeric(result1$se)))
  result1$disease <- disease
  result1$N <- nrow(T2DSNPs_outcome)

  df1 <- rbind(df1,result1)
}


#Plot results
p <-  ggplot(df1, mapping =aes(x = odds_ratio, y = reorder(outcome, desc(outcome)))) +
      geom_line(stat = "identity")+
      geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
      geom_errorbarh(aes(xmax = ci_upper, xmin = ci_lower), size = .5, height = .2, color = "gray50") +
      geom_point(size = 3.5, color = "blue") +
      theme_bw()+
      theme(panel.grid.minor = element_blank(), axis.text=element_text(size=12)) +
      ylab("") +
      xlab("Odds ratio") +
      scale_x_continuous(trans='log10') +
      ggtitle(paste(disease, trait, sep="-"))
print(p)
ggsave(filename = paste(disease,"_", trait, ".pdf", sep=""))
dev.off()

#Append results 
df_all <- rbind(df_all,df1)

#Write results
write.table(df_all, file=paste("T2D_MR_", disease, "_", trait, ".txt", sep=""), sep ='\t', col.names=T, row.names=F, quote=F)

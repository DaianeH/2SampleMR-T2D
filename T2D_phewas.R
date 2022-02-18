library(PheWAS)
library(data.table)

geno<-fread("/path/to/phewas_T2D_tissue_grouped_GRSs.txt",header=T)
phecodes<-fread("/path/to/PheWAS_BioMe_.txt")
geno$MASKED_MRN <- as.integer(geno$MASKED_MRN)
pheno<-merge(geno,phecodes,by.x="MASKED_MRN",by.y=c("id"))
phenotypes=names(pheno)[-(1:40)]

df <- data.frame(stringsAsFactors = F)
tissue <- unique(pheno$trait_tissue)
for (i in 1:(length(tissue))){
   disease_by_disease <- pheno[pheno$trait_tissue==tissue[1],]
   trait <- tissue[1]
   #Run phewas
   results=phewas(phenotypes=phenotypes,genotypes=c("grs"), covariates=c("Age","Sex","median_BMI", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), min.records=10, data=disease_by_disease,cores=3,significance.threshold=c("bonferroni","fdr"))
   results_d=addPhecodeInfo(results)
   sig_results <- results_d[results_d$bonferroni&!is.na(results_d$p),]

   #Write results, only significant and also complete table
   write.table(sig_results,file=paste(trait, ".txt", sep=""), sep = '\t', col.names = T, row.names = F, quote = F)
   write.table(results_d, file=paste(trait, "_all.txt", sep=""), sep='\t', col.names=T, row.names=F, quote=F)

   #Plot results, with legends on top 10
   results_top <- (results[order(results$p),]) 
   phewas_plot <- phewasManhattan(results, OR.direction = T, title=paste("Tissue: ", trait, sep=""), annotate.size=3, annotate.level=results_top[order(results_top$p),"p"][10])
   phewas_plot
   dev.off()
}

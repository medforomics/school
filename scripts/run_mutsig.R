#/cm/shared/apps/R/gcc/3.6.1/bin/Rscript

library(MutationalPatterns)
library(tidyverse)

divisionRel<-function(df){
  sum_df<-sapply(df,sum)
  for (i in 1:ncol(df)){
     df[,i]<-100*round((df[,i]/sum_df[i]),3)
  }
  return(df)
}

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

pathdir <- '.'

vcf_files <- list.files(path=pathdir,pattern = "snps.vcf")
labs <- gsub(".snps.vcf","", vcf_files, perl=TRUE)
vcfs <- read_vcfs_as_granges(vcf_files, labs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

cancer_signatures <- read.table("/project/shared/bicf_workflow_ref/seqprg/musica/aux_files/signatures_probabilities.txt", sep = "\t", header = TRUE)
proposed_etiology <- read.table("/project/shared/bicf_workflow_ref/seqprg/musica/aux_files/proposed_etiology_COSMIC_signatures.txt",sep="\t",header=FALSE)[,2]
known_cancer_signatures<-read.table("/project/shared/bicf_workflow_ref/seqprg/musica/aux_files/cancermatrix.tsv",header=TRUE,sep="\t",row.names=1)

cancer_signatures_aux <- cancer_signatures[order(cancer_signatures[,1]),]
cancer_signatures <- as.matrix(cancer_signatures_aux[,4:33])
cancer_signatures_mut_types <- as.matrix(cancer_signatures_aux[,1:3])

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
relcontribution <- divisionRel(as.data.frame(fit_res$contribution))
tbl <- data.frame(Signature = 1:length(proposed_etiology), Proposed_Etiology = proposed_etiology, relcontribution)
tbl <- tbl[tbl[,3] > 0,] 
write.table(tbl,file = "mutational_signature.txt", sep = "\t", quote=F, row.names=F)
BiocManager::install("MutationalPatterns")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)

vcf_files <- c("../200508_A00295_0339_AHNTYJDMXX/PMABM000GYN_PMABM000GYL_PMCRZ442MBG_WXS_PASS.vcf","../200508_A00295_0339_AHNTYJDMXX/PMABM000GYN_PMABM000GYL_PMCRZ442MBG_WXS_PASS.vcf")

sample_names <- c("PMABM000GYN","PMABM000GYN2")

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(vcfs)

type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

plot_96_profile(mut_mat, condensed = TRUE)


#sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table("sigProfiler_exome_SBS_signatures.csv", sep = ",", header = TRUE)
somaticType <- sapply(c(1:nrow(cancer_signatures)), function(x) paste0(substr(cancer_signatures[x,2],1,1),"[",cancer_signatures[x,1],"]",substr(cancer_signatures[x,2],3,3)))
new_order = match(row.names(mut_mat), somaticType)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = somaticType
cancer_signatures = as.matrix(cancer_signatures[,3:67])

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 10)

plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")

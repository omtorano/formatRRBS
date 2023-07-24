## code to prepare `DATASET` dataset goes here
gtf <- read.delim("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/GCF_016745375.1_EPA_FHM_2.0_genomic.gtf.gz",
                  comment.char="#", sep="\t", header=FALSE)
gtf_gene <- gtf[grep("gene", gtf$V3), ]
usethis::use_data(gtf_gene, internal = TRUE)

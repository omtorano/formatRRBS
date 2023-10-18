## code to prepare `DATASET` dataset goes here
gtf <- read.delim("C:/Users/otorano/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/EE2/GCF_016745375.1_EPA_FHM_2.0_genomic.gtf.gz",
                  comment.char="#", sep="\t", header=FALSE)
gtf_gene <- gtf[grep("gene", gtf$V3), ]
gtf_gene <- suppressWarnings(tidyr::separate(gtf_gene, V9, " ", into = c('1', 'V9'))[, c(1, 2, 3, 4, 5, 6, 7, 8, 10)])
gtf_gene$V9 <- gsub(";", "", gtf_gene$V9)

usethis::use_data(gtf_gene, internal = TRUE, overwrite = TRUE)

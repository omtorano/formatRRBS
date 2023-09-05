#' Format percent methylation from RRBS files for multi omic level integration moli
#'
#' @param x A list of RRBS data frames
#' @param meta A metadata file with rownames matching names of x and a TRT column
#' @param coverage Percent coverage (0.8 = 80%)
#' @param window Window of interest in bases upstream and downstream of genes
#'
#' @return A formatted .RDS file for use in moli, and csv of stats about features removed
#' @export
format_moliRRBS <- function(x, meta, coverage = 0.8, window = 1000) {

  # Collect stats on the features removed with each filter
  stat_rrbs <- t(as.data.frame(lapply(x, nrow)))

  # Remove features with <10 counts
  for (i in names(x)) {
    x[[i]]<-x[[i]][(x[[i]]$V5 + x[[i]]$V6) >= 10, ]
  }
  stat_rrbs <- cbind(stat_rrbs, t(as.data.frame(lapply(x, nrow))))

  # Create unique feature id col
  for(i in 1:length(x)) {
    x[[i]][, 7] <- paste(x[[i]][, 1], x[[i]][, 2])
    colnames(x[[i]])[4] <- names(x)[i]
  }

  # Limit based on percent coverage
  meta <- meta[rownames(meta) %in% names(x), ]
  match_features <- x[[1]][, 7]
  for (i in names(x)[-1]) {
    match_features<-c(match_features, x[[i]][, 7])
  }
  obj <- table(match_features)

  # Need to have at least %coverage*total samples, rounded up to be conservative
  obj_sub <- names(obj[obj >= ceiling(coverage * nrow(meta))])
  # Now need to check the distribution of features within groups
  for(j in unique(meta$TRT)) {
    test_list<-vector()
    for ( i in names(x[names(x) %in% meta[meta$TRT == j, 1]])){
      test_list<-append(test_list,x[[i]]$"V7"[x[[i]]$"V7" %in% obj_sub])
    }
    t_tl <- table(test_list)
    obj_sub <- names(t_tl[t_tl >= ceiling(length(meta[meta$TRT == j, 1]) * coverage)])
  }

  # Limit features to those with %coverage
  sub_x <- list()
  for(i in names(x)) {
    sub_x[[i]]<-x[[i]][x[[i]][, 7] %in% obj_sub,]
  }
  stat_rrbs <- cbind(stat_rrbs, t(as.data.frame(lapply(sub_x, nrow))))

  # Find chromosomes that are 1k upstream and downstream
  obj_sub1 <- tidyr::separate(data = as.data.frame(obj_sub),
                              as.data.frame(obj_sub)[, 1], " ",
                              into = c("chrom", "loc")
                              )

  # Need to search within specific chrom, numbers are not unique
  # The below loop might not label the chromosome correctly, chrom will be
  # relabeled as downstream chrom if chromosomes are less than window bp apart.
  obj_sub2 <- obj_sub1[gtf_gene$V1 %in% obj_sub1$chrom, ]
  obj_sub2$keep <- rep(0, length(obj_sub2$chrom))
  for (i in unique(obj_sub2$chrom)) {
    a <- gtf_gene[gtf_gene$V1 %in% i, ]
    b <- obj_sub2[obj_sub2$chrom == i, ]
    for (j in 1:length(a[, 1])) {
      k <- b$loc[as.numeric(b$loc) >= a[j, 4] - window & as.numeric(b$loc) <= a[j, 5] + window]
      if(length(c) > 0) {
        obj_sub2[obj_sub2$chrom == i & obj_sub2$loc %in% k, 3] <- gtf_gene$V9[j]
      }
    }
  }
  obj_sub3 <- obj_sub2[!obj_sub2$keep==0, ]
  obj_sub3 <- paste(obj_sub3$chrom, obj_sub3$loc)

  # Limit to features within window downstream to window upstream of gene
  sub_x <- list()
  for (i in names(x)) {
    sub_x[[i]] <- x[[i]][x[[i]][, 7] %in% obj_sub3, ]
  }
  stat_rrbs <- cbind(stat_rrbs,t(as.data.frame(lapply(sub_x,nrow))))
  colnames(stat_rrbs) <- c("total", "min10ct", "percent_cov", "gene_region")
  rrbs_merge <- merge(sub_x[[1]][, c(4, 7)], sub_x[[2]][, c(4, 7)], by="V7", all = TRUE)
  for (i in names(sub_x)[-c(1:2)]){
    rrbs_merge<-merge(rrbs_merge, sub_x[[i]][, c(4, 7)], by = "V7", all = TRUE)
  }
  rownames(rrbs_merge) <- rrbs_merge$V7
  rrbs_merge <- rrbs_merge[, -1]

  # Double check
  for (i in trts) {
    print(i)
    print(1 - max(rowSums(is.na(rrbs_merge[, grep(i,colnames(rrbs_merge))]))) / ncol(rrbs_merge[, grep(i,colnames(rrbs_merge))]))
  }
  print(paste((ncol(rrbs_merge) - max(rowSums(is.na(rrbs_merge)))) / nrow(meta),
              "missing values across all samples,", coverage, "specified"), sep = " "
        )

  # Save output
  saveRDS(rrbs_merge, paste0(ctrl, "_", window, "_", coverage, ".RDS"))
  utils::write.csv(stat_rrbs, paste0(ctrl, window, "_", coverage, "stats.csv"))
  }

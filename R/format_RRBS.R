#' Format percent methylation from RRBS files for multi-omic integration
#'
#' @param x A list of RRBS data frames
#' @param meta A metadata file with rownames matching names of x and a TRT column
#' @param total_count Minimum counts (methylated + unmethylated) at each feature, default = 10
#' @param coverage Coverage (enter as decimal, 0.8 = 80%)
#' @param window Window of interest in bases upstream and downstream of genes, requires upstream and downstream value, default 1000, 1000
#' @param values Return beta (100*(M/(M+U)), gives beta distribution), M-values (log2((M+1)/(U+1)), Gaussian), or all
#'
#' @return A formatted .RDS file for use in moli, and csv of stats about features removed
#' @export
format_moliRRBS <- function(x, meta, total_count = 10, coverage = 0.8, window = c(1000, 1000), values = c("beta", "M", "all")) {

  #stop if metadata not provided
  if(is.null(meta)) stop("No metadata provided")

  #stop if metadata not provided
  if(is.null(meta$TRT)) stop("No TRT column, fix metadata")

  # Collect stats on the features removed with each filter
  stat_rrbs <- t(as.data.frame(lapply(x, nrow)))

  # Remove features with < specified total counts
  for (i in names(x)) {
    x[[i]]<-x[[i]][(x[[i]]$V5 + x[[i]]$V6) >= total_count, ]
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
    for ( i in names(x[names(x) %in% meta[meta$TRT == j, 1]])) {
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

  # Remove features that are all un/methylated across all samples

  obj_sub <-

  # Find chromosomes that are in specified upstream and downstream window
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
      k <- b$loc[as.numeric(b$loc) >= a[j, 4] - window[2] & as.numeric(b$loc) <= a[j, 5] + window[1]]
      if(length(c) > 0) {
        obj_sub2[obj_sub2$chrom == i & obj_sub2$loc %in% k, 3] <- gtf_gene$V9[j]
      }
    }
  }
  obj_sub3 <- obj_sub2[!obj_sub2$keep == 0, ]
  obj_sub3 <- paste(obj_sub3$chrom, obj_sub3$loc)

  # Limit to features within window downstream to window upstream of gene
  sub_x <- list()
  for (i in names(x)) {
    sub_x[[i]] <- x[[i]][x[[i]][, 7] %in% obj_sub3, ]
  }
  stat_rrbs <- cbind(stat_rrbs,t(as.data.frame(lapply(sub_x,nrow))))
  colnames(stat_rrbs) <- c("total", "min10ct", "percent_cov", "gene_region")
  # Output options
  if (values == "beta") {
    rrbs_merge <- merge(sub_x[[1]][, c(4, 7)], sub_x[[2]][, c(4, 7)], by="V7", all = TRUE)
    for (i in names(sub_x)[-c(1:2)]) {
      rrbs_merge <- merge(rrbs_merge, sub_x[[i]][, c(4, 7)], by = "V7", all = TRUE)
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
    saveRDS(rrbs_merge, paste0(ctrl, "_", window, "_", coverage, "_beta.RDS"))
    utils::write.csv(stat_rrbs, paste0(ctrl, "_", window, "_", coverage, "stats.csv"))
  }

  if (values == "M") {
    for (i in names(sub_x)) {
      sub_x[[i]][, 8] <- log2((sub_x[[i]][, 5] + 1) / (sub_x[[i]][, 6] + 1))
      colnames(sub_x[[i]])[8] <- paste0(colnames(sub_x[[i]])[4], "_Mval")
    }
    rrbs_merge <- merge(sub_x[[1]][, c(8, 7)], sub_x[[2]][, c(8, 7)], by="V7", all = TRUE)
    for (i in names(sub_x)[-c(1:2)]) {
      rrbs_merge <- merge(rrbs_merge, sub_x[[i]][, c(8, 7)], by = "V7", all = TRUE)
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
    saveRDS(rrbs_merge, paste0(ctrl, "_", window, "_", coverage, "_Mval.RDS"))
    utils::write.csv(stat_rrbs, paste0(ctrl, "_", window, "_", coverage, "stats.csv"))
  }
  if (values == "all") {
    for (i in names(sub_x)) {
      sub_x[[i]][, 8] <- log2((sub_x[[i]][, 5] + 1) / (sub_x[[i]][, 6] + 1))
      colnames(sub_x[[i]])[8] <- paste0(colnames(sub_x[[i]])[4], "_Mval")
    }
    rrbs_merge <- merge(sub_x[[1]][, c(4, 8, 7)], sub_x[[2]][, c(4, 8, 7)], by="V7", all = TRUE)
    for (i in names(sub_x)[-c(1:2)]) {
      rrbs_merge <- merge(rrbs_merge, sub_x[[i]][, c(4, 8, 7)], by = "V7", all = TRUE)
    }
    rownames(rrbs_merge) <- rrbs_merge$V7
    rrbs_merge <- rrbs_merge[, -1]

    # Double check
    for (i in trts) {
      paste0(i," "
            (1 - max(rowSums(is.na(rrbs_merge[, grep(i,colnames(rrbs_merge))]))) / ncol(rrbs_merge[, grep(i,colnames(rrbs_merge))]))*100,
             "%"," missing")
    }
    print(paste0((((ncol(rrbs_merge) - max(rowSums(is.na(rrbs_merge)))) / nrow(meta)) / 2)*100,
                "%"," missing values across all samples,", coverage*100, "%", " specified"))
    # Save output
    saveRDS(rrbs_merge[, grep("Mval", colnames(rrbs_merge))], paste0(ctrl, "_", window, "_", coverage, "_Mval.RDS"))
    saveRDS(rrbs_merge[, grep("Mval", invert = TRUE, colnames(rrbs_merge))], paste0(ctrl, "_", window, "_", coverage, "_beta.RDS"))
    utils::write.csv(stat_rrbs, paste0(ctrl, "_", window, "_", coverage, "stats.csv"))
  }



  }
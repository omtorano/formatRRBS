#' Format percent methylation from RRBS files for multi-omic integration
#'
#' @param x A list of RRBS data frames
#' @param meta A metadata file with rownames matching names of x and a TRT column
#' @param total_count Minimum counts (methylated + unmethylated) at each feature, default = 10
#' @param coverage Coverage (enter as decimal, 0.8 = 80%)
#' @param window Window of interest in bases upstream and downstream of genes, requires upstream and downstream value, default = 1000, 1000
#'
#' @return Two formatted .RDS files of beta (100*(M/(M+U)) gives beta distribution) and M-values (log2((M+1)/(U+1)), Gaussian) for use in downstream analysis, and csv of stats about features removed.
#' @export
format_RRBS <- function(x, meta, total_count = 10, coverage = 0.8, window = c(1000, 1000)) {

  #stop if data are not provided
  if(is.null(x)) stop("No data provided")

  #stop if metadata not provided
  if(is.null(meta)) stop("No metadata provided")

  #stop if metadata not provided
  if(is.null(meta$TRT)) stop("No TRT column, fix metadata")

  #stop if metadata rownames =/= data names
  #if(rownames(meta) != names(x)) stop("No TRT column, fix metadata")

  # Collect stats on the features removed with each filter
  stat_rrbs <- t(as.data.frame(lapply(x, nrow)))

  # Remove features with < specified total counts
  for (i in names(x)) {
    x[[i]]<-x[[i]][(x[[i]]$V5 + x[[i]]$V6) >= total_count, ]
  }
  stat_rrbs <- cbind(stat_rrbs, t(as.data.frame(lapply(x, nrow))))

  # Create unique feature id col
  for(i in 1:length(x)) {
    x[[i]][, 7] <- paste(x[[i]][, 1], x[[i]][, 2], sep = "-")
    colnames(x[[i]])[4] <- names(x)[i]
  }

  # Limit based on percent coverage
  meta <- meta[rownames(meta) %in% names(x), ]
  obj <- table(unlist(lapply(x, function(i) i[, 'V7', drop = FALSE])))

  # Need to have at least %coverage*total samples, rounded up to be conservative
  obj_sub <- names(obj[obj >= ceiling(coverage * nrow(meta))])
  # Now need to check the distribution of features within groups
  for(j in unique(meta$TRT)) {
    test_list <- vector()
    for ( i in names(x[names(x) %in% rownames(meta[meta$TRT == j, ])])) {
      test_list <- append(test_list,x[[i]]$"V7"[x[[i]]$"V7" %in% obj_sub])
    }
    t_tl <- table(test_list)
    obj_sub <- names(t_tl[t_tl >= ceiling(length(rownames(meta[meta$TRT == j, ])) * coverage)])
  }

  # Limit features to those with % coverage
  sub_x <- list()
  for(i in names(x)) {
    sub_x[[i]] <- x[[i]][x[[i]][, 7] %in% obj_sub, ]
  }

  stat_rrbs <- cbind(stat_rrbs, t(as.data.frame(lapply(sub_x, nrow))))

  # Find chromosomes that are in specified upstream and downstream window
  obj_sub1 <- tidyr::separate(data = as.data.frame(obj_sub),
                              as.data.frame(obj_sub)[, 1], "-",
                              into = c("chrom", "loc")
                              )

  # Need to search within specific chrom, numbers are not unique
  obj_sub2 <- obj_sub1[gtf_gene$V1 %in% obj_sub1$chrom, ]
  obj_sub2$gene <- rep("NA", length(obj_sub2$chrom))
  obj_sub2$chrom_start <- rep("NA", length(obj_sub2$chrom))
  obj_sub2$chrom_end <- rep("NA", length(obj_sub2$chrom))
  for (i in unique(obj_sub2$chrom)) {
    v1 <- as.numeric(obj_sub2[obj_sub2$chrom %in% i, 2]) # numeric feature location from RRBS unique features
    dlist <- lapply(v1, function(x) which(data.table::between(x, gtf_gene[gtf_gene$V1 %in% i, 4] - window[1], gtf_gene[gtf_gene$V1 %in% i, 5] + window[2])))
    obj_sub2[obj_sub2$chrom %in% i , ][rep(seq(length(dlist)), lengths(dlist)), 3] <- gtf_gene[gtf_gene$V1 %in% i, ][unlist(dlist), 9]
    obj_sub2[obj_sub2$chrom %in% i , ][rep(seq(length(dlist)), lengths(dlist)), 4] <- gtf_gene[gtf_gene$V1 %in% i, ][unlist(dlist), 4]
    obj_sub2[obj_sub2$chrom %in% i , ][rep(seq(length(dlist)), lengths(dlist)), 5] <- gtf_gene[gtf_gene$V1 %in% i, ][unlist(dlist), 5]
  }

  utils::write.csv(obj_sub2, paste0(ctrl, "_", window[1], "-", window[2], "_", coverage, "_", "Unique_Features_to_Genes.csv"))
  obj_sub3 <- obj_sub2[!obj_sub2$gene == "NA", ]
  rownames(obj_sub3) <- paste(obj_sub3$chrom, obj_sub3$loc, sep = "-")

  # Limit to features within window downstream to window upstream of gene
  sub_x <- list()
  for (i in names(x)) {
    sub_x[[i]] <- x[[i]][x[[i]][, 7] %in% rownames(obj_sub3), ]
  }
  stat_rrbs <- cbind(stat_rrbs, t(as.data.frame(lapply(sub_x, nrow))))

  # Merge and limit to variable features
  rrbs_merge <- dplyr::full_join(sub_x[[1]][, c(4, 7)], sub_x[[2]][, c(4, 7)], by = "V7", keep = FALSE)
  for (i in names(sub_x)[-c(1:2)]) {
    rrbs_merge <- dplyr::full_join(rrbs_merge, sub_x[[i]][, c(4, 7)], by = "V7", keep = FALSE)
  }
  rownames(rrbs_merge) <- rrbs_merge$V7
  rrbs_merge <- rrbs_merge[, -2]

  rrbs_merge <- rrbs_merge[!rowSums(rrbs_merge, na.rm = TRUE) == 0, ]
  rrbs_merge <- rrbs_merge[!matrixStats::rowRanges(as.matrix(rrbs_merge), na.rm = TRUE)[, 1] == 100, ]

  stat_rrbs <- cbind(stat_rrbs, t(nrow(rrbs_merge) - as.data.frame(lapply(rrbs_merge, function(y) sum(length(which(is.na(y))))))))
  colnames(stat_rrbs) <- c("total", "min10ct", "percent_cov", "gene_region", "variable_counts")

  # Double check
  for (i in trts) {
    print(i)
    print(1 - max(rowSums(is.na(rrbs_merge[, grep(i,colnames(rrbs_merge))]))) / ncol(rrbs_merge[, grep(i,colnames(rrbs_merge))]))
  }
  print(paste((ncol(rrbs_merge) - max(rowSums(is.na(rrbs_merge)))) / nrow(meta),
              "missing values across all samples,", coverage, "specified"), sep = " ")

  # Save beta output
  saveRDS(rrbs_merge, paste0(ctrl, "_", window[1], "-", window[2], "_", coverage, "_beta.RDS"))
  utils::write.csv(stat_rrbs, paste0(ctrl, "_", window[1], "-", window[2], "_", coverage, "stats.csv"))

  # Calculate and save M-values
  for (i in names(sub_x)) {
    sub_x[[i]][, 8] <- log2((sub_x[[i]][, 5] + 1) / (sub_x[[i]][, 6] + 1))
    colnames(sub_x[[i]])[8] <- paste0(colnames(sub_x[[i]])[4], "_Mval")
  }
  rrbs_merge_M <- dplyr::full_join(sub_x[[1]][, c(8, 7)], sub_x[[2]][, c(8, 7)], by = "V7", keep = FALSE)
  for (i in names(sub_x)[-c(1:2)]) {
    rrbs_merge_M <- dplyr::full_join(rrbs_merge_M, sub_x[[i]][, c(8, 7)], by = "V7", keep = FALSE)
  }
  rownames(rrbs_merge_M) <- rrbs_merge_M$V7
  rrbs_merge_M <- rrbs_merge_M[, -2]
  rrbs_merge_M <- rrbs_merge_M[rownames(rrbs_merge), ]

  # Save output
  saveRDS(rrbs_merge_M, paste0(ctrl, "_", window[1], "-", window[2], "_", coverage, "_Mval.RDS"))


  }

#' Import RRBS files for formatting
#'
#' @param ctrl Character(s) denoting control treatment in file names
#' @param trts List of character(s) denoting treatment groups in file names
#' @param dir Directory with bismark.cov.gz files
#'
#' @return A list of imported files called 'x', ctrl and trts variables
#' @export
import_moliRRBS <- function(ctrl, trts, dir = getwd()) {
  temp <- list.files(dir, full.names = TRUE)
  x <- lapply(temp[grepl(paste(c(ctrl, trts), collapse = '|'), temp)], data.table::fread, data.table = FALSE)
  names(x) <- gsub("_bismark_bt2.bismark.cov.gz", "",
                   list.files(dir, full.names = FALSE)[grepl(paste(c(ctrl, trts), collapse = '|'),
                   list.files(dir, full.names = TRUE))]
                   )
  assign("x", x, envir = .GlobalEnv)
  assign("ctrl", ctrl, envir = .GlobalEnv)
  assign("trts", trts, envir = .GlobalEnv)
}



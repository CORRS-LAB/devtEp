#' Select Site for the given file
#' p-Value, t-Stats, abds_diffs
#' @param file_name The file name of the data
#' @examples
#' file_name = tempfile()
#' data(demo_devEpt_1k)
#' write.table(demo_devEpt_1k[1:10, ], file_name, row.names=F,  sep = '\t', quote = FALSE)
#' uu <- selectSite(file_name)
selectSite <- function(file_name) {
    ref.dat <- read.table(file_name, header = TRUE, nrow = 2)
    df <- ref.dat[, -c(1, 2)]
    celltypes <- read.table(text = colnames(df), sep = ".")[, 2]
    u_celltypes <- sort(unique(celltypes))
    uu <- getresult(file_name)
    colnames(uu[[1]]) <- u_celltypes
    colnames(uu[[2]]) <- u_celltypes
    return(list(
        pvalue = (1 - (uu[[1]] |> pt(df = 2 * length(u_celltypes) - 1))) * 2,
        ttest = uu[[1]],
        abs_diffs = uu[[2]]
    ))
}

#' Demo Dataset for selectSite
#' 
#' There are 1023 rows and 15 columns in the dataset.
#' @name demo_devEpt_1k
#' @docType data
NULL
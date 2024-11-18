#' Select top N rows from a data.table/ data.frame
#' @param input The input data.table/ data.frame
#' @param N The number of rows to select
#' @param decreasing Whether to select the top N rows or the bottom N rows
#' @return The the selected indices of the top N rows
#' @examples
#' dt.input <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
#' select_top_N(dt.input, 2)
#' select_top_N(dt.input, 2, decreasing = FALSE)
select_top_N <- function(input, N, decreasing = TRUE) {
    return(
        data.matrix(
            unique(
                unlist(
                    apply(
                        input,
                        2,
                        function(x) order(x, decreasing = decreasing)[1:N])))))
}

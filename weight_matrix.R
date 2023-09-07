#' Calculate the weight matrix of cells
#' @description Calculate the cells weight matrix by leveraging highly
#' variable genes
#' @usage weight_matrix(V, DE_index, nu=10000)
#' @param V A gene by cell gene expression matrices
#' @param DE_index The index/names of highly variable genes
#' @param nu A parameter in the heat kernel formula
#' @return The cell weight matrix
#' @examples
#' data(target_data)
#' names_variable <- select_genes(list(t(target_data)), ngenes=200, verbose=F)
#' Q_cell <- weight_matrix(V=target_data, names_variable)
#' @export
#' @importFrom stats dist

weight_matrix <- function(V, DE_index, nu=10000){
    M <- V[DE_index,]
    sim_matrix_euc <- dist(t(M), method = "euclidean", diag = TRUE, upper=TRUE)
    Q <- exp(-sim_matrix_euc/nu)
    Q <- as.matrix(Q)
    return(Q)
}


































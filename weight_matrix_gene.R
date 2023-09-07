#' Calculate the weight matrix of genes
#' @description Calculate the genes weight matrix
#' @usage weight_matrix_gene(X, nv=100)
#' @param X A gene by cell gene expression matrices
#' @param nv A parameter in the heat kernel formula
#' @return The gene weight matrix
#' @examples
#' data(target_data)
#' Q_gene <- weight_matrix_gene(X=target_data)
#' @export
#' @import scater
#' @importFrom stats dist

weight_matrix_gene <- function(X, nv=100){
    X.pca <- scater::calculatePCA(t(log2(X+1)))
    sim_matrix <- dist(X.pca, method = "euclidean", diag = TRUE, upper=TRUE)
    sim_matrix <- as.matrix(sim_matrix)
    Q <- exp(-sim_matrix/nv)
    return(Q)
}






























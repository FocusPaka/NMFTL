#' Select hightly variable genes for cell graphic regularization term
#' @description Select a part of variable genes used to obtain
#' cell weight matrix
#' @usage select_genes(X.list, ngenes = 200, verbose = FALSE)
#' @param X.list A list of cells-by-genes gene expression matrices
#' @param ngenes Integer, number of highly variable genes to select
#' @param verbose Boolean scalar, whether to show Seurat gene selection
#' messages (default FALSE)
#' @return A character that includes the names of ngenes selected variable
#' genes
#' @examples
#' data(target_data)
#' names_variable <- select_genes(list(t(target_data)), ngenes=200,
#' verbose=FALSE)
#' @export
#' @import Seurat


select_genes <- function(X.list, ngenes = 200, verbose = FALSE) {
    m = length(X.list)
    obj.list = lapply(1:m, function(j) {
        obj = Seurat::CreateSeuratObject(counts = t(X.list[[j]]))
        obj = Seurat::NormalizeData(obj, verbose = verbose)
        obj = Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                           nfeatures = ngenes,
                                           verbose = verbose)
        obj
    })
    selected.genes = Seurat::SelectIntegrationFeatures(object.list = obj.list,
                                                       nfeatures = ngenes,
                                                       verbose = TRUE)
    return(selected.genes)
}

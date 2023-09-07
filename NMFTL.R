#' A method that impute scRNA-seq data with non-negative matrix factorization
#' and transfer learning
#' @description A iteratively algorithm that used to realize NMFTL method
#' @usage NMFTL(X_src, X_tar, L_gene, L_cell, S, rank, options)
#' @param X_src A genes-by-cells source gene expression matrices
#' @param X_tar A genes-by-cells target gene expression matrices
#' @param L_gene The Laplacian matrix of genes
#' @param L_cell The Laplacian matrix of cells
#' @param S The transfer matrix from source data to target data
#' @param rank The rank of factorized matrix
#' @param options A sets of parameters should be set beforhand, including
#' initialization of Nonnegative Matrix Factorization Algorithms
#' (`U`,`V_src`,`V_tar`), `verbose`, the maximum iteration times `max_epoch`,
#' parameters in model (`lambda_1`,`lambda_2`,`lambda_3`)
#' @return list(.) return the shared gene representation matrix U, the cell
#' representation matrix V_src from source data, the cell representation matrix
#' V_tar from target data ...
#' @examples
#' data(source_data)
#' data(target_data)
#' ncells1 <- dim(source_data)[2]
#' ncells2 <- dim(target_data)[2]
#' genes2 <- select_genes(list(t(target_data)), ngenes=200, verbose=FALSE)
#' dat <- cbind(source_data,target_data)
#' Q <- weight_matrix(V=target_data, genes2)
#' Deg <- diag(colSums(Q))
#' L_cell <- Deg - Q
#' Q_gene <- weight_matrix_gene(X=target_data)
#' Deg_gene <- diag(colSums(Q_gene))
#' L_gene <- Deg_gene - Q_gene
#' options <- list()
#' rank <- ncol(dat)/10
#' init <- NNDSVD(abs(dat), rank, seed_method='nndsvd',densify='none')
#' options$U <- init$W
#' options$V_src <- (init$H)[,1:ncells1]
#' options$V_tar <- (init$H)[,(ncells1+1):(ncells1+ncells2)]
#' options$verbose <- 2
#' options$max_epoch <- 100
#' options$lambda_1 <- 0
#' options$lambda_2 <- 0.1      # CELL
#' options$lambda_3 <- 0.1       # GENE
#' S_pearson <- cor(source_data, target_data, method = 'pearson')
#' result <- NMFTL(X_src=source_data, X_tar=target_data, L_gene, L_cell,
#' S=abs(S_pearson), rank, options)
#' @export

NMFTL <- function(X_src, X_tar, L_gene, L_cell, S, rank, options){
    cat('# NMFTL: started ...\n')
    # initialize factors
    U <- options$U
    V_src <- options$V_src
    V_tar <- options$V_tar
    A <- diag(rowSums(S))
    B <- diag(colSums(S))

    # initialize
    epoch <- 1

    # cost
    f_val <- differ <- numeric(options$max_epoch)
    f_val[epoch] <- norm(X_tar- U %*% V_tar, type = 'F')^2 / 2
    differ[epoch] <- f_val[epoch]

    if(options$verbose > 1){
        cat(paste0('NMFTL: Epoch = 1, cost = ', f_val[epoch], ', optgap =',
                   abs(differ[epoch]), '\n'))
    }

    # main loop
    while((eval(abs(differ[epoch])) > 100) & (epoch <= options$max_epoch)){
        # update U
        temp1 <- options$lambda_3 * L_gene %*% U
        temp1_m <- pmin(temp1, 0)
        temp1_p <- pmax(temp1, 0)
        U <- U * (X_src %*% t(V_src)+X_tar %*% t(V_tar)-temp1_m) /
            (U %*% V_src %*% t(V_src) + U %*% V_tar %*% t(V_tar) + temp1_p)
        U <- U + (U < 1e-16) * 1e-16

        # update V_src
        V_src <- V_src * (t(U) %*% X_src +  options$lambda_1 * V_tar %*% t(S))/
            (t(U) %*% U %*% V_src + options$lambda_1 * V_src %*% A)
        V_src <- V_src + (V_src < 1e-16) * 1e-16

        # update V_tar
        temp2 <- options$lambda_2 * V_tar %*% L_cell
        temp2_m <- pmin(temp2, 0)
        temp2_p <- pmin(temp2, 0)
        V_tar <- V_tar * (t(U) %*% X_tar + options$lambda_1 *
                              V_src %*% S - temp2_m) / (t(U) %*% U %*% V_tar +
                                                            temp2_p +
                                                            options$lambda_1 *
                                                            V_tar %*% B)
        V_tar <- V_tar + (V_tar < 1e-16) * 1e-16

        epoch <- epoch + 1
        f_val[epoch] <- norm(X_tar- U %*% V_tar, type = 'F')^2 / 2
        differ[epoch] <- f_val[epoch] - f_val[epoch-1]

        if(options$verbose > 1){
            cat(paste0('NMFTL: Epoch = ', epoch, ', cost = ', f_val[epoch],
                       ', optgap = ',abs(differ[epoch]), '\n'))
        }
    }
    return(list(U = U, V_src = V_src, V_tar = V_tar, cost = f_val,
                difference = differ))
}

#' The initialization of the non-negative matrix factorization
#' @description This function implements the NNDSVD algorithm for
#  initialization of Nonnegative Matrix Factorization Algorithms.
#' @usage NNDSVD(A, k, seed_method='nndsvd', densify=c('none', 'average', 'random'))
#' @param A the input nonnegative m x n matrix A
#' @param k the rank of the computed factors W,H
#' @param seed_method A character, defalt is nndsvd
#' @param densify A character that select which algorithme
#' @return A nonnegative m x k matrix W, a nonnegative k x n matrix H
#' @examples
#' data(source_data)
#' data(target_data)
#' dat <- cbind(source_data, target_data)
#' init <- NNDSVD(abs(dat), rank=ncol(dat)/10, seed_method='nndsvd',
#' densify='none')
#' @export
#' @importFrom stats runif


NNDSVD <- function(A, k, seed_method='nndsvd', densify=c('none', 'average', 'random')){
    # % This function implements the NNDSVD algorithm described in [1] for
    # % initialization of Nonnegative Matrix Factorization Algorithms.
    # %
    # % [W,H] = nndsvd(A,k,flag);
    # %
    # % INPUT
    # % ------------
    #     %
    # % A    : the input nonnegative m x n matrix A
    # % k    : the rank of the computed factors W,H
    # % flag : indicates the variant of the NNDSVD Algorithm
    # %        flag = 0 --> NNDSVD
    # %        flag = 1 --> NNDSVDa
    # %        flag = 2 --> NNDSVDar
    # %
    # % OUTPUT
    # % -------------
    #     %
    # % W   : nonnegative m x k matrix
    # % H   : nonnegative k x n matrix
    # %
    # %
    # % References:
    #     %
    # % [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    # %     start for nonnegative matrix factorization, Pattern Recognition,
    # %     Elsevier
    # %
    # % This code is kindly provided by the authors for research porpuses.
    # % - Efstratios Gallopoulos (stratis@ceid.upatras.gr)
    # % - Christos Boutsidis (boutsc@cs.rpi.edu)
    # %
    # % For any problems or questions please send an email to boutsc@cs.rpi.edu
    # %--------------------------------------------------------------------------

    if(seed_method=="nndsvd"){
        densify <- match.arg(densify)
        flag <- which(densify == c('none', 'average', 'random')) - 1
        # %----------------------check the input matrix------------------------------
        # % n = numel(A) 返回数组 A 中的元素数目 n 等同于 prod(size(A))。
        if(any(A<0)){
            cat('The input matrix contains negative elements !')
        }

        # %size of the input matrix
        m <- nrow(A)
        n <- ncol(A)

        # %the matrices of the factorization
        W <- matrix(0, m, k)
        H <- matrix(0, k, n)

        # % 1st SVD --> partial SVD rank-k to the input matrix A.
        s <- RSpectra::svds(A, k);
        U <- s$u; S <- s$d; V <- s$v
        # dim(U) (150,5) length(S) 5 dim(V) (150, 5)

        # %choose the first singular triplet to be nonnegative
        W[,1] <- sqrt(S[1]) * abs(U[,1]);
        H[1,] <- sqrt(S[1]) * abs(t(V[,1]));

        # % 2nd SVD for the other factors (see table 1 in our paper)
        for (i in 2:k) {
            uu <- U[,i]; vv <- V[,i];
            uup <- .pos(uu); uun <- .neg(uu);
            vvp <- .pos(vv); vvn <- .neg(vv);
            n_uup <- .norm(uup);
            n_vvp <- .norm(vvp);
            n_uun <- .norm(uun) ;
            n_vvn <- .norm(vvn) ;
            termp <- n_uup * n_vvp; termn <- n_uun * n_vvn;
            if(termp >= termn){
                W[,i] <- sqrt(S[i] * termp) * uup / n_uup;
                H[i,] <- sqrt(S[i] * termp) * vvp / n_vvp;
            }else{
                W[,i] <- sqrt(S[i] * termn) * uun / n_uun;
                H[i,] <- sqrt(S[i] * termn) * vvn / n_vvn;
            }
        }

        # %------------------------------------------------------------
        # %actually these numbers are zeros
        W[W<0.0000000001] <- 0.1;
        H[H<0.0000000001] <- 0.1;

        # % NNDSVDa: fill in the zero elements with the average
        if(flag == 1){
            ind1 <- W==0
            ind2 <- H==0
            average <- mean(A)
            W[ind1] <- average
            H[ind2] <- average
        }else if(flag == 2){
            # % NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
            ind1 <- W==0
            ind2 <- H==0
            n1 <- sum(ind1)
            n2 <- sum(ind2)

            average <- mean(A)
            W[ind1] <- average * runif(n1, min = 0, max = 1) / 100
            H[ind2] <- average * runif(n2, min = 0, max = 1) / 100
        }

        # return matrices W and H
        out<-list(W=W, H=H)
    }else{
        #W<-matrix(rbinom(k*nrow(A),round(max(A)/10),0.5),ncol = k)
        W <- matrix(runif(k*nrow(A),0,1),nrow = nrow(A),ncol = k)
        #H<-matrix(runif(k*ncol(A)),nrow = k)
        H <- matrix(runif(k*ncol(A),0,1),nrow = k,ncol=ncol(A))
        out <- list(W=W, H=H)
    }
    return(out)
}



###% Auxliary functions
.pos <- function(x){ as.numeric(x>=0) * x }
.neg <- function(x){ - as.numeric(x<0) * x }
.norm <- function(x){ sqrt(drop(crossprod(x))) }







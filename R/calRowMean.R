#'@title Compute the means of replicates
#'@description Duplicated franctions A and B are
#'summarized by taking their mean for each protein.
#'After taking the mean, the data log2 transformed.
#'Further, the 5 main fractions are used to check correlation
#'between input datas. It is a helper function.
#'@param d.df data.frame; A data frame of 10 fraction profiles
#' consisting of replicate A and B.
#'@export
#'@examples {
#'
#'r.df <- calRowMean(SubCellBarCode::hcc827Ctrl)
#'
#'}
#'@return r.df

calRowMean <- function(d.df){
    r.means <- lapply(seq_len(5), function(x){
        k <- 2 * x -1
        t.df <- rowMeans(d.df[, c(k:(k+1))])
    })

    r.df <- data.frame(do.call("cbind", r.means))
    colnames(r.df) <- c("Cyto", "Nsol", "NucI", "Horg", "Lorg")
    r.df <- log2(r.df)
    return(r.df)
}


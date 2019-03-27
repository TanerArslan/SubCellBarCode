#'@title Load the fractionated proteomics data
#'@description Sampled median normalized TMT ratios are checked if
#'there is any "NA" valeus. If any, the corresponding
#'row is filtered out. Later, the data is normalized by taking log2.
#'@param protein.data data.frame; fractionated proteomics data where
#'data contains 10 columns of duplicated 5 fractionations and
#'rownames must be gene-centric protein names
#'@importFrom stats na.omit
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl[1:10,])
#'}
#'@return protein.data.df

loadData <- function(protein.data){

    if(! is.data.frame(protein.data))
        stop('Input must be a data frame format! Type ?loadData')

    if(! ncol(protein.data) == 10)
        stop('Input data must have 10 columns! Type ?loadData')

    if (! is.character(rownames(protein.data)))
        stop('Rownames must be character!')


    if(sum(is.na(protein.data)) > 0){

        protein.data <- na.omit(protein.data)
        if(sum(protein.data < 0) > 0){
            protein.data.df <- protein.data
            warning("Input data.frame has negative values.
Presumably, data was normalized. Therefore, log normalization
was not performed.")
        }else{
            protein.data.df <- log2(protein.data)
            cat("Rows with NA values were filtered.")
        }


    }else{
        if(sum(protein.data < 0) > 0){
            protein.data.df <- protein.data
            warning("Input data.frame has negative values.
Presumably, data was normalized. Therefore, log
normalization was not performed.")
        }else{
            protein.data.df <- log2(protein.data)
        }

    }
    return(protein.data.df)
}


#'@title Convert identifier to gene symbol
#'@description Identifier for each feature should be converted
#'into gene symbols unless they are not gene symbols
#'@param df data.frame; fractionated proteomics data where
#'data contains 10 columns of duplicated 5 fractionations and
#'rownames must be identifier e.g. UNIPROT, Entrez ID
#'@param id caharacter; identifier id for each protein
#'@import org.Hs.eg.db
#'@import AnnotationDbi
#'@export
#'@examples {
#'
#'df <- data.frame(Uniprot = c("A4D0S4","A8TX70","O00305","O00337"),
#'Organism = rep("Homo Sap.", 4))
#'
#'rownames(df) <- df$Uniprot
#'}
#'
#'@return df
#'


convert2symbol <- function(df, id = "UNIPROT"){

    if (! id %in% columns(org.Hs.eg.db))
        stop('Enter valid identifier.')

    if (! is.character(rownames(df)))
        stop('Rownames must be character!')

    uniprot.id <- as.character(unique(rownames(df)))

    #converting ids
    annot.df<- AnnotationDbi::select(org.Hs.eg.db,
                            keys = uniprot.id,
                            keytype = id,
                            columns = "SYMBOL")

    ## filter-1 ids
    #filter if there is no corresponding symbol
    annot.df <- annot.df[complete.cases(annot.df),]

    #filter if there multiple id corresponding one one symbol
    #remove one of them
    annot.df <- annot.df[!duplicated(annot.df$SYMBOL),]

    #filter-3 if one id matches to multiple symbol
    dup.ids <- annot.df$UNIPROT[duplicated(annot.df$UNIPROT)]
    annot.df <- subset(annot.df, ! annot.df$UNIPROT %in% dup.ids)

    df <- df[as.character(annot.df$UNIPROT),]
    rownames(df) <- annot.df$SYMBOL

    return(df)
}


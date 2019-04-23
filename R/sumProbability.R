#'@title Sum compartment test data probabilities to
#'neighborhood probabilities
#'@description Compartment levels classifications on the  test data are
#'summed up to associated neighborhood levels. It is a helper function.
#'@param df data.frame; test data classifications at the neighborhood
#'level and probablity vectors for each protein.
#'@export
#'@examples {
#'
#'#create mock data
#'df <- data.frame(Protein = "TP53",
#'svm.pred = "N1",
#'S1 = as.numeric(0.02),
#'S2 = as.numeric(0.02),
#'S3 = as.numeric(0.02),
#'S4 = as.numeric(0.02),
#'N1 = as.numeric(0.72),
#'N2 = as.numeric(0.02),
#'N3 = as.numeric(0.02),
#'N4 = as.numeric(0.02),
#'C1 = as.numeric(0.02),
#'C2 = as.numeric(0.02),
#'C3 = as.numeric(0.02),
#'C4 = as.numeric(0.02),
#'C5 = as.numeric(0.02),
#'M1 = as.numeric(0.02),
#'M2 = as.numeric(0.02))
#'
#'rownames(df) <- "TP53"
#'
#'sum.df <- sumProbability(df)
#'
#'}
#'@return summed.df

sumProbability <- function(df){

    t.secretory.df <- data.frame(df[, colnames(df)[3:6]])
    t.secretory.df$Secretory <- rowSums(t.secretory.df)
    t.nuclear.df <- data.frame(df[, colnames(df)[7:10]])
    t.nuclear.df$Nuclear <- rowSums(t.nuclear.df)
    t.cytosol.df <- data.frame(df[, colnames(df)[11:15]])
    t.cytosol.df$Cytosol <- rowSums(t.cytosol.df)
    t.Mitochondria.df <- data.frame(df[, colnames(df)[16:17]])
    t.Mitochondria.df$Mitochondria <- rowSums(t.Mitochondria.df)

    summed.df <- data.frame(Proteins = rownames(df),
                            df[,colnames(df)[seq_len(2)]],
                            Secretory = t.secretory.df$Secretory,
                            Nuclear = t.nuclear.df$Nuclear,
                            Cytosol = t.cytosol.df$Cytosol,
                            Mitochondria = t.Mitochondria.df$Mitochondria)
    #temp neighborhood df
    t.n.df <- summed.df[,4:7]
    summed.df$svm.pred <- colnames(t.n.df)[apply(t.n.df, 1, which.max)]
    return(summed.df)
}

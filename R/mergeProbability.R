#'@title Merge compartment probabilities to neighborhood probabilities
#'@description Compartment levels classifications are summed up to
#'associated neighborhood levels. It is a helper function.
#'@param df data.frame; all predictions at the neighborhood level and
#'probablity vectors for each protein
#'@export
#'@examples {
#'
#'#create mock data
#'df <- data.frame(Protein = "TP53",
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
#'merged.df <- mergeProbability(df)
#'
#'}
#'@return merged.df

mergeProbability <- function(df){

    t.secretory.df <- data.frame(df[, colnames(df)[2:5]])
    t.secretory.df$Secretory <- apply(t.secretory.df, 1, sum)
    t.nuclear.df <- data.frame(df[, colnames(df)[6:9]])
    t.nuclear.df$Nuclear <- apply(t.nuclear.df, 1, sum)
    t.cytosol.df <- data.frame(df[, colnames(df)[10:14]])
    t.cytosol.df$Cytosol <- apply(t.cytosol.df, 1, sum)
    t.Mitochondria.df <- data.frame(df[, colnames(df)[15:16]])
    t.Mitochondria.df$Mitochondria <- apply(t.Mitochondria.df, 1, sum)

    merged.df <- data.frame(Proteins = rownames(df),
                            svm.pred.all = df[,colnames(df)[1]],
                            Secretory = t.secretory.df$Secretory,
                            Nuclear = t.nuclear.df$Nuclear,
                            Cytosol = t.cytosol.df$Cytosol,
                            Mitochondria = t.Mitochondria.df$Mitochondria)
    # temp neihborhood df
    t.n.df <- merged.df[,3:6]
    merged.df$svm.pred.all <- colnames(t.n.df)[apply(t.n.df, 1, which.max)]
    rownames(merged.df) <- merged.df$Proteins
    return(merged.df)
}

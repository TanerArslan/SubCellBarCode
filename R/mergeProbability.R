#'@title Merge compartment probabilities to neighborhood probabilities
#'@description Compartment levels classifications are summed up to
#'associated neighborhood levels. It is a helper function.
#'@param df data.frame; all predictions at the neighborhood level and
#'probablity vectors for each protein
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'set.seed(7)
#'c.prots <- sample(c.prots, 365)
#'cls <- svmClassification(c.prots, df, markerProteins)
#'
#'all.A <- cls[[1]]$all.prot.pred
#'
#'all.n.repA <- replacePrediction(all.A, column = "svm.pred.all")
#'
#'m.all.repA <- mergeProbability(all.n.repA)
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

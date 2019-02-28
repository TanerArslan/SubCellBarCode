#'@title Merge compartment and neighborhood classification
#'@description Compartment and neighborhood classifications are merged
#'for the single output.
#'@param compartmentCls data.frame; all predictions, including
#'unclassified as well, and
#'probablity vectors for each protein in compartment classification
#'@param neighborhoodCls data.frame; all predictions, including
#'unclassified as well, and
#'probablity vectors for each protein in compartment classification
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'set.seed(1)
#'
#'cls <- svmClassification(c.prots[1:700], df, markerProteins)
#'
#'test.A <- cls[[1]]$svm.test.prob.out
#'test.B <- cls[[2]]$svm.test.prob.out
#'
#'t.c.df <- computeThresholdCompartment(test.A, test.B)
#'
#'t.n.df <- computeThresholdNeighborhood(test.A, test.B)
#'
#'all.A <- cls[[1]]$all.prot.pred
#'all.B <- cls[[2]]$all.prot.pred
#'
#'c.cls.df <- applyThresholdCompartment(all.A, all.B, t.c.df)
#'
#'n.cls.df <- applyThresholdNeighborhood(all.A, all.B, t.n.df)
#'
#'cls.df <- mergeCls(c.cls.df, n.cls.df)
#'}
#'@return cls.df


mergeCls <- function(compartmentCls, neighborhoodCls){

    # get the same order of the proteins
    compartmentCls <- compartmentCls[rownames(neighborhoodCls),]

    if( ! identical(rownames(compartmentCls), rownames(neighborhoodCls)) )
        stop('Make sure your inputs files posses same rownames')


    cls.df <- data.frame(Protein = compartmentCls$Proteins,
                        NeighborhoodCls = neighborhoodCls$svm.pred.all,
                        CompartmentCls = compartmentCls$svm.pred,
                        neighborhoodCls[,3:6],
                        compartmentCls[,3:17])

    cls.df <- cls.df[order(cls.df$Protein, decreasing = FALSE),]

}

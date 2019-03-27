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
#'#create mock data
#'com.df <- data.frame(Proteins = "TP53",
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
#'rownames(com.df) <- "TP53"
#'
#'neig.df <- data.frame(Proteins = "TP53",
#'svm.pred.all = "Nuclear",
#'Secretory = as.numeric(0.01),
#'Nuclear = as.numeric(0.95),
#'Cytosol = as.numeric(0.02),
#'Mitochondria = as.numeric(0.02))
#'
#'rownames(neig.df) <- "TP53"
#'
#'cls.df <- mergeCls(com.df, neig.df)
#'
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

#'@title Apply thresholds to neighborhood classification
#'@description Apply thresholds for all predictions at the
#'neighborhood level to increase the
#'true positive rate and remove poor classification.
#'@param all.repA data.frame; all predictions and probablity vectors
#'for each protein in replicate A
#'@param all.repB data.frame; all predictions and probablity vectors
#'for each protein in replicate B
#'@param threshold.df data.frame; collection od precision and recall
#'values for each neighborhood
#'@importFrom stats aggregate
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'set.seed(7)
#'c.prots <- sample(c.prots, 500)
#'cls <- svmClassification(c.prots, df, markerProteins)
#'
#'test.A <- cls[[1]]$svm.test.prob.out
#'test.B <- cls[[2]]$svm.test.prob.out
#'
#'t.n.df <- computeThresholdNeighborhood(test.A, test.B)
#'
#'all.A <- cls[[1]]$all.prot.pred
#'all.B <- cls[[2]]$all.prot.pred
#'
#'n.cls.df <- applyThresholdNeighborhood(all.A, all.B, t.n.df)
#'}
#'@return n.cls.df

applyThresholdNeighborhood <- function(all.repA, all.repB, threshold.df){

    #upgrade compartment labels to neighborhood labels for prediction
    all.n.repA <- SubCellBarCode::replacePrediction(df = all.repA,
                                            column = "svm.pred.all")
    all.n.repB <- SubCellBarCode::replacePrediction(df = all.repB,
                                            column = "svm.pred.all")

    #sum up compartment level predictions to neighborhood predictions
    m.all.repA <- SubCellBarCode::mergeProbability(all.n.repA)
    m.all.repB <- SubCellBarCode::mergeProbability(all.n.repB)

    m.all.repB <- m.all.repB[rownames(m.all.repA), ]

    #merge the replicates by averaging the probabilities
    all.repAB.mean <- rbind(m.all.repA, m.all.repB)[,-2]
    all.repAB.mean <- aggregate(.~Proteins, data = all.repAB.mean, mean)
    all.repAB <- data.frame(Proteins = all.repAB.mean$Proteins,
                    svm.pred.all = rep("Unclassified", nrow(all.repAB.mean)),
                    all.repAB.mean[,2:5])
    rownames(all.repAB) <- all.repAB$Proteins
    all.repAB <- all.repAB[rownames(m.all.repA), ]

    #merge two replicates and average them
    n.repA.m <- m.all.repA[m.all.repA$svm.pred.all == m.all.repB$svm.pred.all, ]
    n.repB.m <- m.all.repB[m.all.repB$svm.pred.all == m.all.repB$svm.pred.all, ]

    combined.reps <- rbind(n.repA.m, n.repB.m)
    combined.df <- data.frame(Proteins = combined.reps$Proteins,
                                combined.reps[, 3:6])
    averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
    rownames(averaged.reps) <- averaged.reps$Proteins
    averaged.reps <- averaged.reps[rownames(n.repA.m), ]

    combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                    svm.pred.all = n.repA.m$svm.pred.all,
                                    averaged.reps[, 2:5])

    #apply the theresholds for each neighborhood
    neighborhoods <- c("Secretory","Nuclear", "Cytosol", "Mitochondria")

    confident.classification <- lapply(neighborhoods, function(m){
        # temp precision
        t.p <- unname(unlist(threshold.df[threshold.df$Neighborhood == m, ][2]))
        #temp recall
        t.r <- unname(unlist(threshold.df[threshold.df$Neighborhood == m, ][3]))
        if (! is.na(t.p)){
            t.value <- max(t.p, t.r)
            temp.df <- combined.rep.A.B[combined.rep.A.B$svm.pred.all == m, ]
            up.threshold.df <- temp.df[temp.df[m] >= t.value, ]
        }
    })

    conf.df <- do.call("rbind", confident.classification)

    ##adding "unclassified proteins"
    no.class <- subset(all.repAB, ! rownames(all.repAB) %in% rownames(conf.df))

    n.cls.df <- rbind(conf.df, no.class)

}



#'@title Apply thresholds to compartments
#'@description Apply thresholds for all predictions to increase the
#'true positive rate and remove poor classification.
#'@param all.repA data.frame; all predictions and probablity vectors
#'for each protein in replicate A
#'@param all.repB data.frame; all predictions and probablity vectors
#'for each protein in replicate B
#'@param threshold.df data.frame; collection od precision and recall
#'values for each compaartment
#'@importFrom stats aggregate
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'set.seed(7)
#'c.prots <- sample(c.prots, 550)
#'cls <- svmClassification(c.prots, df, markerProteins)
#'
#'test.A <- cls[[1]]$svm.test.prob.out
#'test.B <- cls[[2]]$svm.test.prob.out
#'
#'t.c.df <- computeThresholdCompartment(test.A, test.B)
#'
#'all.A <- cls[[1]]$all.prot.pred
#'all.B <- cls[[2]]$all.prot.pred
#'
#'c.cls.df <- applyThresholdCompartment(all.A, all.B, t.c.df)
#'}
#'@return c.cls.df

applyThresholdCompartment <- function(all.repA, all.repB, threshold.df){

    if( ! identical(rownames(all.repA), rownames(all.repB)) )
        stop('Replicates have to posses same indexed rownames')

    all.repA$Proteins <- rownames(all.repA)
    all.repB$Proteins <- rownames(all.repB)

    #merge the replicates by averaging the probabilities
    repAB.mean <- rbind(all.repA, all.repB)[,-1]
    repAB.mean <- aggregate(.~Proteins, data = repAB.mean, mean)
    repAB <- data.frame(Proteins = repAB.mean$Proteins,
                        svm.pred = rep("Unclassified", nrow(repAB.mean)),
                        repAB.mean[,seq_len(16)])
    rownames(repAB) <- repAB$Proteins
    repAB <- repAB[rownames(all.repA), ]

    all.repA.match <- all.repA[all.repA$svm.pred.all == all.repB$svm.pred.all, ]
    all.repB.match <- all.repB[all.repB$svm.pred.all == all.repA$svm.pred.all, ]
    combined.repAB <- rbind(all.repA.match, all.repB.match)
    combined.df <- data.frame(Proteins = combined.repAB$Proteins,
                                combined.repAB[, 2:16])
    averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
    rownames(averaged.reps) <- averaged.reps$Proteins
    averaged.reps <- averaged.reps[rownames(all.repA.match), ]
    averaged.repAB <- data.frame(Proteins = averaged.reps$Proteins,
                                svm.pred = all.repA.match$svm.pred.all,
                                averaged.reps[, 2:16])


    #apply the theresholds for each compartment
    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                        "C1", "C2", "C3", "C4", "C5", "M1", "M2")

    confident.classification <- lapply(compartments, function(m){
        #temp precision
        t.p <- unname(unlist(threshold.df[threshold.df$Compartment == m, ][2]))
        #temp recall
        t.r <- unname(unlist(threshold.df[threshold.df$Compartment == m, ][3]))
        t.value <- max(t.p, t.r)
        if(is.numeric(t.value)){
            temp.df <- averaged.repAB[averaged.repAB$svm.pred == m, ]
            up.threshold.df <- temp.df[temp.df[m] >= t.value, ]
        }
    })

    confident.df <- do.call("rbind", confident.classification)
    confident.df <- confident.df <- confident.df[complete.cases(confident.df),]

    # adding "unclassified proteins"
    no.clss <- subset(repAB, ! rownames(repAB) %in% rownames(confident.df))[-3]

    c.cls.df <- rbind(confident.df, no.clss)

}

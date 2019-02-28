#'@title Probability threshold for neighborhood classification
#'@description Thresholds for each neighborhood are decided to get
#' confident predictions.
#'@param test.repA data.frame; test predictions, observation and
#'probablity vectors for each protein in replicate A
#'@param test.repB data.frame; test predictions, observation and
#'probablity vectors for each protein in replicate B
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
#'t.n.df <- computeThresholdNeighborhood(test.A, test.B)
#'}
#'@return threshold.neighborhood.df

computeThresholdNeighborhood <- function(test.repA, test.repB){

    couple.lsit <- list(c("Secretory", "S1"), c("Secretory", "S2"),
                        c("Secretory", "S3"), c("Secretory", "S4"),
                        c("Nuclear", "N1"), c("Nuclear", "N2"),
                        c("Nuclear", "N3"), c("Nuclear", "N4"),
                        c("Cytosol", "C1"), c("Cytosol", "C2"),
                        c("Cytosol", "C3"), c("Cytosol", "C4"),
                        c("Cytosol", "C5"), c("Mitochondria", "M1"),
                        c("Mitochondria", "M2"))

    #upgrade compartment labels to neighborhood labels
    replaceRows <- function(df, column = c("Observation", "svm.pred")){
        multiple.lst <- lapply(couple.lsit, function(f){
            temp.df <- df[df[column] == unname(unlist(f[2])), ]
            temp.df[[column]] <- as.character(unname(unlist(f[1])))
            temp.df
        })
        replaced.df <- do.call("rbind", multiple.lst)
    }

    replaceObervation <- replaceRows(df = test.repA, column = "Observation")
    neighborhood.repA <- replaceRows(df = replaceObervation,
                                        column = "svm.pred")

    replaceObervation <- replaceRows(df = test.repB, column = "Observation")
    neighborhood.repB <- replaceRows(df = replaceObervation,
                                        column = "svm.pred")

    # concatanate the probabilities for corresponding neighborhood
    sumProbability <- function(df){

        t.secretory.df <- data.frame(df[, colnames(df)[3:6]])
        t.secretory.df$Secretory <- apply(t.secretory.df, 1, sum)
        t.nuclear.df <- data.frame(df[, colnames(df)[7:10]])
        t.nuclear.df$Nuclear <- apply(t.nuclear.df, 1, sum)
        t.cytosol.df <- data.frame(df[, colnames(df)[11:15]])
        t.cytosol.df$Cytosol <- apply(t.cytosol.df, 1, sum)
        t.Mitochondria.df <- data.frame(df[, colnames(df)[16:17]])
        t.Mitochondria.df$Mitochondria <- apply(t.Mitochondria.df, 1, sum)

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

    sum.repA <- sumProbability(neighborhood.repA)
    sum.repB <- sumProbability(neighborhood.repB)

    sum.repB <- sum.repB[rownames(sum.repA), ]

    #merge two replicates and average them
    n.test.repA.match <- sum.repA[sum.repA$svm.pred == sum.repB$svm.pred, ]
    n.test.repA.match$Proteins <- rownames(n.test.repA.match)
    n.test.repB.match <- sum.repB[sum.repB$svm.pred == sum.repA$svm.pred, ]
    n.test.repB.match$Proteins <- rownames(n.test.repB.match)

    combined.reps <- rbind(n.test.repA.match, n.test.repB.match)
    combined.df <- data.frame(Proteins = combined.reps$Proteins,
                                combined.reps[, 4:7])
    averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
    rownames(averaged.reps) <- averaged.reps$Proteins
    averaged.reps <- averaged.reps[rownames(n.test.repA.match), ]

    combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                    Observation = n.test.repA.match$Observation,
                                    svm.pred = n.test.repA.match$svm.pred,
                                    averaged.reps[, 2:5])

    #estimate the neihborhood theresholds by calculating precision and recall
    cls.levels <- c("Secretory", "Nuclear", "Cytosol", "Mitochondria")

    results <- lapply(cls.levels, function(l){
        cls.df <- combined.rep.A.B[combined.rep.A.B$svm.pred == l, ]
        cls.obs.df <- neighborhood.repA[neighborhood.repA$Observation == l, ]
        parameters <- lapply(seq(0, 1, 0.002), function(t){
            u.df <- cls.df[cls.df[l] >= t, ]
            p.cls <- sum(u.df$Observation == u.df$svm.pred)/nrow(u.df)
            cls.down.df <- cls.df[cls.df[l] < t,]
            r.cls <- (sum(u.df$Observation == u.df$svm.pred))/nrow(cls.obs.df)
            f.score <- (2 * p.cls * r.cls) / (p.cls + r.cls)
            values <- list(Precision = p.cls,
                            Recall = r.cls,
                            fscore= f.score,
                            threshold = t,
                            Compartment = l)

        })
        result.df <- data.frame(do.call(rbind.data.frame, parameters))
        up.precision <- result.df[result.df$Precision >= 0.95, ]
        recall.threshold <- result.df[!duplicated(result.df[, c('Recall')]), ]
        threshold.df <- list(Neighborhood = l,
                                Precision = up.precision$threshold[1],
                                Recall = recall.threshold$threshold[2])
    })

    threshold.neighborhood.df <- data.frame(do.call(rbind, results))

    colnames(threshold.neighborhood.df)[2:3] <- c("PrecisionBasedThreshold",
                                                    "RecallBasedThreshold")
    threshold.neighborhood.df$OptedThreshold <- apply(threshold.neighborhood.df,
                                1, function(x) max(unlist(x[2]), unlist(x[3])))
    return(threshold.neighborhood.df)


}

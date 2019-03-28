#'@title Probability threshold for compartment classification
#'@description Thresholds for each compartment are decided to get
#'confident predictions.
#'@param test.repA data.frame; test predictions, observation and
#'probablity vectors for each protein in replicate A
#'@param test.repB data.frame; test predictions, observation and
#'probablity vectors for each protein in replicate B
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
#'}
#'@return threshold.compartment.df

computeThresholdCompartment <- function(test.repA, test.repB){

    if( ! identical(rownames(test.repA), rownames(test.repB)) )
        stop('Replicates have to posses same indexed rownames')

    #merge two replicates and average them
    test.repA.match <- test.repA[test.repA$svm.pred == test.repB$svm.pred, ]
    test.repA.match$Proteins <- rownames(test.repA.match)
    test.repB.match <- test.repB[test.repA$svm.pred == test.repB$svm.pred, ]
    test.repB.match$Proteins <- rownames(test.repB.match)

    combined.reps <- rbind(test.repA.match, test.repB.match)
    combined.df <- data.frame(Proteins = combined.reps$Proteins,
                                combined.reps[, 3:17])
    averaged.reps <- aggregate(.~Proteins, data = combined.df, mean)
    rownames(averaged.reps) <- averaged.reps$Proteins
    averaged.reps <- averaged.reps[rownames(test.repA.match),]

    combined.rep.A.B <- data.frame(Proteins = averaged.reps$Proteins,
                                    Observation = test.repA.match$Observation,
                                    svm.pred = test.repA.match$svm.pred,
                                    averaged.reps[, 2:16])


    #estimate the compartments theresholds  by calculating precision and recall
    cls.levels <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                        "C1", "C2", "C3", "C4", "C5", "M1", "M2")

    results <- lapply(cls.levels,function(l){
        cls.df <- combined.rep.A.B[combined.rep.A.B$svm.pred==l, ]
        cls.obs.df <- test.repA[test.repA$Observation==l, ]
        if(nrow(cls.df) > 0 & nrow(cls.obs.df) > 0){
            parameters <- lapply(seq(0, 1, 0.005), function(t){
            u.df <- cls.df[cls.df[l] >= t, ]
            p.cls <- sum(u.df$Observation == u.df$svm.pred) / nrow(u.df)
            cls.down.df <- cls.df[cls.df[l] < t, ]
            r.cls <- (sum(u.df$Observation == u.df$svm.pred)) / nrow(cls.obs.df)
            f.score <- (2 * p.cls * r.cls) / (p.cls + r.cls)
            values <- list(Precision = p.cls,
                            Recall = r.cls,
                            fscore = f.score,
                            threshold = t,
                            Compartment = l)
            })

            result.df <- data.frame(do.call(rbind.data.frame, parameters))
            up.precision <- result.df[result.df$Precision>= 0.9, ]
            recall.threshold <- result.df[!duplicated(result.df[,c('Recall')]),]
            threshold.df <- list(Compartment = l,
                                Precision = up.precision$threshold[1],
                                Recall = recall.threshold$threshold[2])
        }})

    threshold.compartment.df <- data.frame(do.call(rbind, results))

    missing.comp <- setdiff(cls.levels, threshold.compartment.df$Compartment)

    if(length(missing.comp) > 1){
        miss.df <- data.frame(Compartment = missing.comp,
                            Precision = rep("NA", length(missing.comp)),
                            Recall = rep("NA", length(missing.comp)))
        threshold.compartment.df <- rbind(threshold.compartment.df, miss.df)
        comps <- threshold.compartment.df$Compartment
        rownames(threshold.compartment.df) <- comps
        threshold.compartment.df <- threshold.compartment.df[cls.levels,]

    }else if(length(missing.comp) == 1){
        miss.df <- data.frame(Compartment = missing.comp,
                            Precision = "NA",
                            Recall = "NA")
        threshold.compartment.df <- rbind(threshold.compartment.df, miss.df)
        comps <- threshold.compartment.df$Compartment
        rownames(threshold.compartment.df) <- comps
        threshold.compartment.df <- threshold.compartment.df[cls.levels,]
    }else{
        threshold.compartment.df <- threshold.compartment.df
    }
    colnames(threshold.compartment.df)[2:3] <- c("PrecisionBasedThreshold",
                                                    "RecallBasedThreshold")
    threshold.compartment.df$OptedThreshold <- apply(threshold.compartment.df,
                                1, function(x) max(unlist(x[2]), unlist(x[3])))

    return(threshold.compartment.df)
}

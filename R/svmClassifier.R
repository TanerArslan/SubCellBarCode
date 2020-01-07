#'@title Protein subcellular localization classification
#'@description Support Vector Machine classifier is trained and
#'used for prediction of protein subcellular localization
#'@param markerProteins character; robust marker proteins along
#' with subcellular localization that are present in the given data.
#'@param protein.data data.frame; fractionated proteomics data
#'@param markerprot.df data.frame; collection of marker proteins
#'along with corresponding subcellular localization
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
#'}
#'@import caret
#'@import e1071
#'@importFrom stats predict
#'@return all.classifications

svmClassification <- function(markerProteins, protein.data, markerprot.df){

    prot.df <- protein.data[markerProteins, ]

    marker.df <- markerprot.df[markerProteins, ]

    subcell.loc.data <- cbind(prot.df, marker.df[seq_len(2)])

    #splitting the data
    label.train <- caret::createDataPartition(subcell.loc.data$Compartments,
                        times = 1,
                        p=0.7,
                        list = FALSE)

    df.train <- subcell.loc.data[label.train, ]
    df.test <- subcell.loc.data[-label.train, ]

    #build classifier for both replicate A and B, respectively.
    build.classifier <- lapply(c("A", "B"), function(x) {
        rep.train.df <- df.train[, grepl(sprintf("\\.%s\\.", x),
                                            colnames(df.train))]
                replicate.train.label <- df.train$Compartments

                rep.test.df <- df.test[, grepl(sprintf("\\.%s\\.", x),
                                            colnames(df.test))]
                replicate.test.label <- df.test$Compartments

                # tune the best parameters
                svm.tune <- e1071::tune(svm,
                                    train.x = rep.train.df,
                                    train.y = factor(replicate.train.label),
                                    kernel="radial",
                                    ranges=list(cost = 10^(seq_len(4) - 2),
                                    gamma = c(.5, 1, 1.5, 2)))

                # build the svm model
                model.svm <- e1071::svm(factor(replicate.train.label) ~ .,
                                        data = rep.train.df,
                                        probability=TRUE,
                                        kernel="radial",
                                        cost=svm.tune$best.parameters[1],
                                        gamma=svm.tune$best.parameters[2],
                                        scale = FALSE,
                                        cross=10)

                # predict and test on test data
                svm.pred <- predict(model.svm, rep.test.df, probability = TRUE)
                svm.pred.test.df <- data.frame(svm.pred)
                svm.test.prob <- data.frame(attr(svm.pred, "probabilities"))
                test.observation <- data.frame(df.test$Compartments)
                svm.test.prob.out <- cbind(test.observation,
                                            svm.pred.test.df,
                                            svm.test.prob)
                colnames(svm.test.prob.out)[1] <- 'Observation'
                svm.test.prob.out <- svm.test.prob.out[,
                                    c("Observation", "svm.pred", "S1", "S2",
                                    "S3", "S4", "N1", "N2", "N3", "N4", "C1",
                                    "C2", "C3", "C4", "C5", "M1", "M2")]

                if(sum(svm.test.prob.out$Observation ==
                        svm.test.prob.out$svm.pred) / nrow(df.test) < 0.50){

                    cat("Overall prediction accuracy is < % 50.
                        Downstream analysis will not be accurate enough.
                        We highly recommend you to perform wet-lab analyis
                        again.")
                }

                # predict all proteins
                all.prot.df <- protein.data[,grepl(sprintf("\\.%s\\.", x),
                                                        colnames(df.train))]
                svm.pred.all <- predict(model.svm,
                                        all.prot.df,
                                        probability = TRUE)

                svm.all.label <- data.frame(svm.pred.all)
                svm.all.prob <- data.frame(attr(svm.pred.all, "probabilities"))
                all.prot.pred <- cbind(svm.all.label, svm.all.prob)
                all.prot.pred <- all.prot.pred[,c("svm.pred.all", "S1", "S2",
                                    "S3", "S4", "N1", "N2", "N3", "N4", "C1",
                                    "C2", "C3", "C4", "C5", "M1", "M2")]

                # merge both test predicitions and all protein predictions
                all.classifications <- list(svm.test.prob.out=svm.test.prob.out,
                                        all.prot.pred = all.prot.pred,
                                        model = model.svm)

                return(all.classifications)

            })
}

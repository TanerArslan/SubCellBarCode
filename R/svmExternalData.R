#'@title Peptide/exon/transcript centric or
#'PTM enriched classification
#'@description Peptide/exon/transcript centric or
#'PTM enriched classification is applied to predict
#'localization of them.
#'@param df, data frame fractionated additional data
#'@param modelA, model for the replicate A classification
#'@param modelB, model for the replicate B classification
#'@import caret
#'@import e1071
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
#'modelA <- cls[[1]]$model
#'modelB <- cls[[2]]$model
#'
#'exon.cls <- svmExternalData(SubCellBarCode::hcc827exon,
#'modelA = modelA, modelB = modelB)
#'}
#'@return c.cls.df

svmExternalData <- function(df, modelA, modelB){

    Cls <- lapply(c("A", "B"), function(x) {
        if(x == "A"){
            all.prot.df <- df[,grepl(sprintf("\\.%s\\.", x),
                                        colnames(df))]

            #take log normalization if it was not
            if(sum(all.prot.df < 0) >= 0){
                all.prot.df <- log2(all.prot.df)
            }

            # get the model
            svm.mod1 <- modelA

            # make prediction
            svm.pred.all <- predict(svm.mod1,
                                    all.prot.df,
                                    probability = TRUE)

            svm.all.label <- data.frame(svm.pred.all)
            svm.all.prob <- data.frame(attr(svm.pred.all, "probabilities"))
            all.prot.pred <- cbind(df[,1], svm.all.label, svm.all.prob)
            colnames(all.prot.pred)[1] <- "Gene_Symbol"
            all.prot.pred <- all.prot.pred[,c("Gene_Symbol","svm.pred.all",
                                            "S1", "S2", "S3", "S4", "N1",
                                            "N2", "N3", "N4", "C1",
                                            "C2", "C3", "C4", "C5",
                                            "M1", "M2")]


        }else{
            all.prot.df <- df[,grepl(sprintf("\\.%s\\.", x),
                                     colnames(df))]

            #take log normalization if it was not
            if(sum(all.prot.df < 0) >= 0){
                all.prot.df <- log2(all.prot.df)
            }

            # get the model
            svm.mod2 <- modelB

            # make prediction
            svm.pred.all <- predict(svm.mod2,
                                    all.prot.df,
                                    probability = TRUE)

            svm.all.label <- data.frame(svm.pred.all)
            svm.all.prob <- data.frame(attr(svm.pred.all, "probabilities"))
            all.prot.pred <- cbind(df[,1], svm.all.label, svm.all.prob)
            colnames(all.prot.pred)[1] <- "Gene_Symbol"
            all.prot.pred <- all.prot.pred[,c("Gene_Symbol","svm.pred.all",
                                            "S1", "S2", "S3", "S4", "N1",
                                            "N2", "N3", "N4", "C1",
                                            "C2", "C3", "C4", "C5",
                                            "M1", "M2")]
        }
    })}

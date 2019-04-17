#'@title Replace compartment predictions to neighborhood predictions
#'@description Compartment level classifications are replaced with
#'neighborhood level assignment. It is a helper function.
#'@param df data.frame; all predictions at the compartment level and
#'probablity vectors for each protein
#'@param column character; selected column in the data frame, df
#'@export
#'@examples {
#'
#'#define mock data frame
#'df <- data.frame(svm.pred.all = c("S1","S2","S3","S4",
#'"N1","N2","N3","N4",
#'"C1","C2","C3","C4","C5",
#'"M1","M2"))
#'
#'df$svm.pred.all <- as.character(df$svm.pred.all)
#'
#'df <- replacePrediction(df, column = "svm.pred.all")
#'}
#'@return replaced.df

replacePrediction <- function(df,
                        column = c("svm.pred.all", "Observation", "svm.pred")){

    couple.lsit <- list(c("Secretory", "S1"), c("Secretory", "S2"),
                        c("Secretory", "S3"), c("Secretory", "S4"),
                        c("Nuclear", "N1"), c("Nuclear", "N2"),
                        c("Nuclear", "N3"), c("Nuclear", "N4"),
                        c("Cytosol", "C1"), c("Cytosol", "C2"),
                        c("Cytosol", "C3"), c("Cytosol", "C4"),
                        c("Cytosol", "C5"), c("Mitochondria", "M1"),
                        c("Mitochondria", "M2"))

    multiple.lst <- lapply(couple.lsit, function(f){
        temp.df <- df[df[column] == unname(unlist(f[2])), ]
        if(nrow(temp.df) > 0){
            temp.df[[column]] <- as.character(unname(unlist(f[1])))
            temp.df
        }
    })
    replaced.df <- do.call("rbind", multiple.lst)
}


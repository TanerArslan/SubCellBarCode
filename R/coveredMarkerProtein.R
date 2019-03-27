#'@title Evaluate marker protein coverage
#'@description Given the proteomics data, number of overlapped marker
#'proteins is calculated.
#'Bar plot for each compartment is plotted.
#'@param proteinIDs  character; gene symbol id
#'@param markerproteins character; 3365 proteins gene symbol ids
#'@import ggplot2
#'@importFrom graphics plot
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'}
#'@return covered.proteins

calculateCoveredProtein <- function(proteinIDs, markerproteins){

    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                        "C1", "C2", "C3", "C4", "C5", "M1", "M2")

    color.code <- c("gold", "orange", "salmon", "tomato2", "grey90",
                "grey70", "grey50", "grey30", "lightblue", "aquamarine",
                "cyan", "deepskyblue2", "turquoise3", "burlywood4", "tan4")

    compartment.size <- c(358, 351, 252, 174, 192, 121, 231, 198, 242,
                            132, 220, 215, 341, 69, 269)


    covered.proteins <- intersect(proteinIDs, markerproteins)

    if( length(covered.proteins) < 1)
        stop('There is no overlap between marker proteins and data!')

    c.marker.df <- SubCellBarCode::markerProteins[covered.proteins, ]

    coverageCompWise <- lapply(seq_len(length(compartments)), function(x){
    temp.df <- c.marker.df[c.marker.df$Compartments == compartments[x], ]
    values <- list(Compartments = compartments[x],
            ColorCode = color.code[x],
            ProteinCoverage = 100 * ((dim(temp.df)[1]) / compartment.size[x]))
    })

    coverage.df <- as.data.frame(do.call("rbind", coverageCompWise))

    #check if there is not enough enrichemnt in any compartment
    non.enriched.loc <- coverage.df[coverage.df$ProteinCoverage < 20, ]
    if(nrow(non.enriched.loc) == 1){
        warning("There is not enough enrichment at: ",
                as.character(non.enriched.loc$Compartments),
                "\nWe recommend you to perform the fractionation, again.")
    }else if(nrow(non.enriched.loc) > 1){
        comp <- paste(as.character(non.enriched.loc$Compartments),
                collapse = ",")
        warning("There are not enough enrichments at: ",
                comp, "\nWe recommend you to perform the fractionation!")
    }


    coverage.df$ProteinCoverage <- as.numeric(coverage.df$ProteinCoverage)
    coverage.df$Compartments <- as.character(coverage.df$Compartments)
    coverage.df$ColorCode <- as.character(coverage.df$ColorCode)

    plot(ggplot(data = coverage.df,
        aes(x = coverage.df$Compartments, y = coverage.df$ProteinCoverage)) +
        geom_bar(stat="identity", fill = coverage.df$ColorCode) +
        scale_x_discrete(limits=c(compartments)) +
        theme_bw() +
        theme(text = element_text(size = 16),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(face = "bold", color="black"),
                axis.text.y = element_text(face = "bold", color="black")) +
        labs(title = "Marker Protein Coverage Compartment-Wise",
            y = "% Protein Coverage",
            x = "Compartment"
        ))

    coverage <- round(length(covered.proteins) / length(markerproteins), 2)
    cat("Overall Coverage of marker proteins : ", coverage)

    return (covered.proteins)
}

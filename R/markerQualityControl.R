#' @title Evaluate the quality of the marker proteins
#' @description Given the proteomics data, quality of the overlapped
#' marker proteins are evaluated by correlating replicates of fractions.
#' @param coveredProteins character; list of marker proteins, gene
#' symbols, that are covered in 3365 marker proteins.
#' @param protein.data data.frame; fractionated proteomics data,
#' rownames are gene symbols associated protein.
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'r.markers <- markerQualityControl(c.prots, df)
#'}
#' @import ggplot2
#' @import gridExtra
#' @importFrom stats cor
#' @return robustMarkers


markerQualityControl <- function(coveredProteins, protein.data){

    # replicate-wise correlation marker QC
    m.prot.df <- protein.data[coveredProteins, ]
    m.prot.df <- m.prot.df[complete.cases(m.prot.df),]

    if( nrow(m.prot.df) < 1)
        stop('Make sure your inputs are correct. Type ?markerQualityControl')

    fracs.df.A <- m.prot.df[, grepl("\\.A\\.", colnames(m.prot.df))]
    fracs.df.B <- m.prot.df[, grepl("\\.B\\.", colnames(m.prot.df))]

    cor.reps.pearson <- vapply(rownames(fracs.df.A),
        function (x) {cor(unlist(fracs.df.A[x, ]), unlist(fracs.df.B[x, ]),
                            method="pearson")},
                            as.numeric(""))

    replicate.df <- data.frame(Protein = names(cor.reps.pearson),
                                Correlation = cor.reps.pearson)

    p1 <- ggplot(replicate.df, aes(x = replicate.df$Correlation)) +
        geom_density(alpha = .7, fill = "deepskyblue") +
        theme_minimal() +
        theme(text = element_text(size = 14),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(face = "bold", color="black"),
                axis.text.y = element_text(face = "bold", color="black")) +
        geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
        labs(title = "Replicate-Wise Correlation Marker QC",
                tag = "A",
                y = "Density",
                x = "Pearson Corr.")

    #remove replicate-wise markerp proteins
    rep.prots <- names(cor.reps.pearson[cor.reps.pearson < 0.8 ])

    message(sprintf("Number of removed replicate-wise proteins: %s",
                    length(rep.prots)))

    # sample-wise correlation marker QC
    prot.names <- setdiff(rownames(m.prot.df), rep.prots)
    markerProteins <- SubCellBarCode::markerProteins[prot.names,][,3:7]
    m.prot.df <- m.prot.df[prot.names,]

    # function for to calculate sample-wise correlation for each protein
    prot.cor <- function(df, marker.df, cor.method = c("spearman", "pearson")){
        unlist(lapply(rownames(m.prot.df), function(x){
            p.cor <- cor(t(df[x,]), t(marker.df[x,]), method = cor.method)
            names(p.cor) <- x
            p.cor
        }))}


    pearson.corA <- prot.cor(df = m.prot.df[grepl("\\.A\\.",
                                                colnames(m.prot.df))],
                                marker.df = markerProteins,
                                cor.method = "pearson")

    pearson.corB <- prot.cor(df = m.prot.df[grepl("\\.B\\.",
                                                colnames(m.prot.df))],
                                marker.df = markerProteins,
                                cor.method = "pearson")

    pearson.cor <- data.frame(RepA = pearson.corA, RepB = pearson.corB)
    pearson.cor$MinP.Cor <- apply(pearson.cor, 1, min)

    spear.corA <- prot.cor(df = m.prot.df[grepl("\\.A\\.",
                                            colnames(m.prot.df))],
                                marker.df = markerProteins,
                                cor.method = "spearman")

    spear.corB <- prot.cor(df = m.prot.df[grepl("\\.B\\.",
                                            colnames(m.prot.df))],
                                marker.df = markerProteins,
                                cor.method = "spearman")

    spear.cor <- data.frame(RepA = spear.corA, RepB = spear.corA)
    spear.cor$MinS.cor <- apply(spear.cor, 1, min)

    df <- data.frame(Protein = rownames(spear.cor),
                        Pearson = pearson.cor$MinP.Cor,
                        Spearman = spear.cor$MinS.cor)

    cols <- SubCellBarCode::markerProteins[prot.names,][8]
    Color <- cols$Colour

    p2 <- ggplot(df, aes(x = df$Pearson, y = df$Spearman)) +
        geom_point(colour = Color, size = 2) +
        geom_hline(yintercept = 0.6, linetype="dashed", color = "red") +
        geom_vline(xintercept = 0.8, linetype="dashed", color = "red") +
        labs(title = "Sample-Wise Correlation Marker QC",
            tag = "B",
            y = "Spearman Corr.",
            x = "Pearson Corr.") +
        theme_minimal() +
        theme(text = element_text(size = 14),
            plot.title = element_text(hjust = 0.5, color = "black", size = 14),
            axis.text.x = element_text(face = "bold", color="black"),
            axis.text.y = element_text(face = "bold", color="black"))

    sample.removed.prot <- df[df$Pearson < 0.8 | df$Spearman < 0.599,]
    sample.removed.prot <- as.character(sample.removed.prot$Protein)

    message(sprintf("Number of removed sample-wise proteins: %s",
                    length(sample.removed.prot)))

    robustMarkerProteins <- setdiff(prot.names, sample.removed.prot)

    message(sprintf("Number of total removed marker proteins: %s",
                    length(sample.removed.prot) + length(rep.prots)))

    grid.arrange(p1, p2, ncol=2)

    # check if there is a depletion after marker qc
    compartment.size <- c(358, 351, 252, 174, 192, 121, 231, 198, 242,
                            132, 220, 215, 341, 69, 269)

    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                        "C1", "C2", "C3", "C4", "C5", "M1", "M2")
    r.marker.df <- SubCellBarCode::markerProteins[robustMarkerProteins, ]

    coverageCompWise <- lapply(seq_len(length(compartments)), function(x){
        temp.df <- r.marker.df[r.marker.df$Compartments == compartments[x], ]
        values <- list(Compartments = compartments[x],
            ProteinCoverage = 100 * ((dim(temp.df)[1]) /compartment.size[x]))
    })

    r.cov.df <- as.data.frame(do.call("rbind", coverageCompWise))

    non.enriched.loc <- r.cov.df[r.cov.df$ProteinCoverage < 20, ]
    if(nrow(non.enriched.loc) == 1){
        warning(sprintf("There is not enough enrichment at %s localization.
                \nWe recommend you to perform the fractionation, again.",
                        as.character(non.enriched.loc$Compartments)))
    }else if(nrow(non.enriched.loc) > 1){
        comp <- paste(as.character(non.enriched.loc$Compartments),
                    collapse = ",")
        warning(sprintf("There are not enough enrichment at %s localizations.
                        \nWe recommend you to perform the fractionation,
    as we describe at the manuscprit.", comp))
    }

    return(robustMarkerProteins)
}

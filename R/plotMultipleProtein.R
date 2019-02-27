#'@title Visualization of multiple protein localizations
#'@description Distributions of subcellular localizations
#'of multiple proteins both ar the
#'compartment and neighborhood level are plotted.
#'@param sampleClassification data.frame; merged classification,
#'combination of compartment and neighborhood classifications
#'per protein.
#'@param proteinList vector; protein gene symbol names.
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'r.markers <- markerQualityControl(c.prots, df)
#'
#'cls <- svmClassification(r.markers, df, markerProteins)
#'
#'test.A <- cls[[1]]$svm.test.prob.out
#'test.B <- cls[[2]]$svm.test.prob.out
#'
#'t.c.df <- computeThresholdCompartment(test.A, test.B)
#'
#'t.n.df <- computeThresholdNeighborhood(test.A, test.B)
#'
#'all.A <- cls[[1]]$all.prot.pred
#'all.B <- cls[[2]]$all.prot.pred
#'
#'c.cls.df <- applyThresholdCompartment(all.A, all.B, t.c.df)
#'
#'n.cls.df <- applyThresholdNeighborhood(all.A, all.B, t.n.df)
#'
#'cls.df  <- mergeCls(c.cls.df, n.cls.df)
#'
#'proteasome26s <- c("PSMA7", "PSMC3", "PSMB1", "PSMA1", "PSMA3","PSMA4",
#'"PSMA5", "PSMB4", "PSMB6", "PSMB5","PSMC2", "PSMC4", "PSMB3", "PSMB2",
#'"PSMD4", "PSMA6", "PSMC1", "PSMC5", "PSMC6", "PSMB7", "PSMD13")
#'
#'multipleProt.df <- plotMultipleProtein(cls.df, "proteasome26s" )
#'
#'}
#'@import ggplot2
#'@import gridExtra
#'@return multipleProt.df


plotMultipleProtein <- function(sampleClassification, proteinList){

    protein.df <- sampleClassification[proteinList, ][, seq_len(3)]

    if( nrow(protein.df) < 1 )
        stop('Thre is no overlap between protein list and data')

    #check proteins that are not in the data
    not.exist.prot <- setdiff(proteinList, rownames(protein.df))

    # export the non-existing proteins
    if(length(not.exist.prot) > 0){
        message(sprintf("%s not exist in data",
            paste0(not.exist.prot, collapse = ",")))
    }

    protein.df <- protein.df[complete.cases(protein.df), ]

    multipleProt.df <- protein.df

    neighborhoods <- c("Secretory", "Nuclear", "Cytosol", "Mitochondria",
                        "Unclassified")
    compartments <- c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4", "C1",
                        "C2", "C3", "C4", "C5", "M1", "M2", "Unclassified")

    neighCls <- lapply(neighborhoods, function(n){

        temp.n.df <- protein.df[protein.df$NeighborhoodCls == n, ]

        temp.df <- data.frame(Neighborhood = n,
                            Count = nrow(temp.n.df))

    })

    neigh.summay.df <- do.call("rbind", neighCls)

    neigh.cols <- c("Secretory" = "goldenrod1", "Nuclear" = "grey80",
                    "Cytosol" = "cadetblue1", "Mitochondria" = "peru",
                    "Unclassified" = "mistyrose")

    #bar plot of the neighborhood level classifications
    p1 <- ggplot(neigh.summay.df,
                aes(x = neigh.summay.df$Neighborhood,
                    y = neigh.summay.df$Count,
                    fill = neigh.summay.df$Neighborhood)) +
            geom_bar(stat = "identity", position = "stack", width = 1) +
            scale_fill_manual(values = c(neigh.cols)) +
            labs(title = " ",
                    tag = "B",
                    y = "Protein Counts",
                    x = "Neighborhoods") +
            theme_minimal() +
            theme(text = element_text(size = 16),
                    legend.position="none",
                    axis.text.x = element_text(face = "bold",
                                                color="black",
                                                angle = 45,
                                                hjust = 1),
                    axis.text.y = element_text(face = "bold", color="black"))


    compCls <- lapply(compartments, function(n){

        temp.n.df <- protein.df[protein.df$CompartmentCls == n, ]

        temp.df <- data.frame(Compartment = n,
                                Count = nrow(temp.n.df))

    })

    comp.summay.df <- do.call("rbind", compCls)

    comp.cols <- cols <- c("S1" = "gold", "S2" = "orange", "S3" = "salmon",
                            "S4" = "tomato2", "N1" = "grey90", "N2" = "grey70",
                            "N3" = "grey50", "N4" = "grey30",
                            "C1" = "lightblue", "C2" = "aquamarine",
                            "C3" = "cyan", "C4" = "deepskyblue2",
                            "C5" = "turquoise3", "M1" = "burlywood4",
                            "M2" = "tan4", "Unclassified" = "mistyrose")

    #bar plot of the compartment level classifications
    p2 <- ggplot(comp.summay.df,
                aes(x = comp.summay.df$Compartment,
                    y = comp.summay.df$Count,
                    fill = comp.summay.df$Compartment)) +
                geom_bar(stat = "identity", position = "stack", width = 1) +
                scale_fill_manual(values = c(comp.cols)) +
                labs(title = " ",
                    tag = "A",
                    y = "Protein Counts",
                    x = "Compartments") +
                theme_minimal() +
                theme(text = element_text(size = 16),
                    legend.position="none",
                    axis.text.x = element_text(face = "bold",
                                                color="black",
                                                angle = 45,
                                                hjust = 1),
                    axis.text.y = element_text(face = "bold",
                                                color="black"))

    grid.arrange(p2, p1, ncol=2)
}

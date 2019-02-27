#'@title Sankey plot for condition-dependent protein relocalization
#'@description Identify candidate condition-dependent relocated
#'proteins by comparing neighborhood classifications.
#'@param sampleCls1 data.frame; merged classification,
#'combination of compartment and neighborhood classification.
#'@param sampleCls2 data.frame; merged classification,
#'combination of compartment and neighborhood classification.
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
#'sankeyData <- sankeyPlot(cls.df, hcc827GEFClass)
#'
#'}
#'@import networkD3
#'@return sankeyData


sankeyPlot <- function(sampleCls1, sampleCls2){

    # just keep neighborhood classification
    df1 <- sampleCls1[, seq_len(2)]
    df2 <- sampleCls2[, seq_len(2)]

    df1 <- df1[intersect(rownames(df1), rownames(df2)), ]
    df2 <- df2[intersect(rownames(df1), rownames(df2)), ]

    df2 <- df2[rownames(df1), ]

    df <- data.frame(Protein = df1$Protein,
                        C.A = df1$NeighborhoodCls,
                        C.B = df2$NeighborhoodCls)

    #remove unclassfied proteins at both conditions
    f.df <- df[!grepl("Unclassified",df$C.A) & ! grepl("Unclassified",df$C.B),]

    neighborhoods <- c("Secretory", "Nuclear", "Cytosol", "Mitochondria")
    label.neigh <- c(0, 1, 2, 3)

    neighborhoods2 <- c("Secretory", "Nuclear", "Cytosol", "Mitochondria")
    label.neigh2 <- c(4, 5, 6, 7)

    # summarize neighborhood to neighborhood localization
    sourceTarget <- lapply(seq_len(4), function(n1){
        temp.n1.df <- df[df$C.A == neighborhoods[n1], ]

        results <- lapply(seq_len(4), function(n2){
                temp.n2.df <- temp.n1.df[temp.n1.df$C.B == neighborhoods2[n2], ]
                values <- list(source = label.neigh[n1],
                                target = label.neigh2[n2],
                                value = nrow(temp.n2.df))
        })
        result.df <- data.frame(do.call(rbind.data.frame, results))
    })

    # source to link for the sankey plot
    link.df <- do.call("rbind", sourceTarget)
    link.df <- link.df[link.df$value > 0, ]
    rownames(link.df) <- seq_len(nrow(link.df))

    # prepare nodes dataframe
    neigh.names <- c("Cond1.Secretory", "Cond1.Nuclear", "Cond1.Cytosol",
                        "Cond1.Mitochondria", "Cond2.Secretory",
                        "Cond2.Nuclear", "Cond2.Cytosol", "Cond2.Mitochondria")

    group <- as.character(c(1, 2, 3, 4, 1, 2, 3, 4))

    id <- c(0, 1, 2, 3, 4, 5, 6, 7)

    sankeyNodes <- data.frame(neigh.names, id, group)

    sankeyData <- list(sankeyNodes, link.df)
    print(sankeyData)

    # plot sankey plot
    networkD3::sankeyNetwork(Links = link.df,
                                Nodes = sankeyNodes,
                                Source = "source",
                                Target = "target",
                                Value = "value",
                                NodeID = "neigh.names",
                                fontSize = 12,
                                nodeWidth = 30,
                                colourScale = JS(
                                    'd3.scaleOrdinal()
                                    .domain(["0", "1", "2", "3", "4",
                                            "5", "6", "7"])
                                    .range(["#FFC125", "#CCCCCC",
                                            "#98F5FF", "#CD853F"])'))

}

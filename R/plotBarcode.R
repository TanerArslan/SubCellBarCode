#'@title Visualize the SubCellBarCode
#'@description Stacked bar plot are plotted for compartment and
#'neighborhood level with respect to classification probabilities.
#'@param sampleClassification data.frame; merged classification,
#'combination of compartment and neighborhood classification.
#'@param protein character; protein gene symbol name
#'@param s1PSM data.frame; minimum PSM count table. Row names should
#'be gene centric protein id.
#'@export
#'@examples {
#'
#'#create mock data
#'plot.df <- data.frame(Protein = "TP53",
#'NeighborhoodCls = "Nuclear",
#'CompartmentCls = "N1",
#'Secretory = as.numeric(0.01),
#'Nuclear = as.numeric(0.95),
#'Cytosol = as.numeric(0.02),
#'Mitochondria = as.numeric(0.02),
#'S1 = as.numeric(0.02),
#'S2 = as.numeric(0.02),
#'S3 = as.numeric(0.02),
#'S4 = as.numeric(0.02),
#'N1 = as.numeric(0.72),
#'N2 = as.numeric(0.02),
#'N3 = as.numeric(0.02),
#'N4 = as.numeric(0.02),
#'C1 = as.numeric(0.02),
#'C2 = as.numeric(0.02),
#'C3 = as.numeric(0.02),
#'C4 = as.numeric(0.02),
#'C5 = as.numeric(0.02),
#'M1 = as.numeric(0.02),
#'M2 = as.numeric(0.02))
#'
#'rownames(plot.df) <- "TP53"
#'
#'psm.df <- data.frame(Protein = "TP53",
#'PSMs.for.quant = as.numeric(31))
#'
#'rownames(psm.df) <- "TP53"
#'
#'proteinPlot <- plotBarcode(plot.df, "TP53", psm.df)
#'}
#'@import ggplot2
#'@importFrom graphics plot
#'@return proteinPlot


plotBarcode <- function(sampleClassification, protein, s1PSM){

    p.df <- sampleClassification[as.character(protein),]

    if( sum(is.na(p.df)) > 0 )
        stop('Invalid protein name.
            Please make sure the protein is in the classification output.')

    comp.df <- data.frame(Locs = as.character(colnames(p.df[8:22])),
                            Probability = as.numeric(unname(t(p.df[1,8:22]))),
                            Level = "Compartment")

    neigh.df <- data.frame(Locs = as.character(colnames(p.df[4:7])),
                            Probability = as.numeric(unname(t(p.df[1,4:7]))),
                            Level = "Neighborhood")

    cols <- c("S1" = "gold", "S2" = "orange", "S3" = "salmon",
                "S4" = "tomato2", "N1" = "grey90", "N2" = "grey70",
                "N3" = "grey50", "N4" = "grey30", "C1" = "lightblue",
                "C2" = "aquamarine", "C3" = "cyan", "C4" = "deepskyblue2",
                "C5" = "turquoise3", "M1" = "burlywood4", "M2" = "tan4",
                "Cytosol" = "cadetblue1", "Secretory" = "goldenrod1",
                "Nuclear" = "grey80", "Mitochondria" = "peru")

    df <- rbind(comp.df, neigh.df)

    df$Locs <- factor(df$Locs,
                    levels = c("Nuclear", "Secretory", "Cytosol",
                                "Mitochondria", " ", "S1", "S2", "S3", "S4",
                                "N1", "N2", "N3", "N4", "C1", "C2", "C3",
                                "C4", "C5", "M1", "M2"))
    #get the PSM count
    psm <- as.numeric(s1PSM[protein,][2])

    if( length(psm) < 1 & !is.numeric(psm))
        stop('PSM count could not obtain properly.
            Please check the PSM input data')

    #calls for classification and probabilities
    c.c <- as.character(p.df[,3])
    if(c.c == "Unclassified"){
        c.p <- round(max(p.df[8:22]), 2)
    }else{
        c.p <- round(p.df[,c.c], 2)
    }


    n.c <- as.character(p.df[,2])
    if(n.c == "Unclassified"){
        n.p <- round(max(p.df[4:7]), 2)
    }else{
        n.p <- round(p.df[,n.c], 2)
    }

    #labels for the plot
    c.l <- sprintf("Call: %s \nProb: %s", c.c, c.p)
    n.l <- sprintf("Call: %s \nProb: %s", n.c, n.p)

    proteinPlot <- ggplot(df,
                        aes(x = df$Level,
                            y = df$Probability,
                            fill = df$Locs,
                            group = df$Probability)) +
                    geom_bar(stat = "identity",
                            position = "stack",
                            width = 0.5) +
                    facet_wrap(.~Level, scales = "free") +
                    scale_fill_manual(values = c(cols[seq_len(4)],
                                        "white",
                                        cols[5:19]),
                                        drop=FALSE) +
                    theme_bw() +
                        theme(text = element_text(size = 14),
                        strip.text.x = element_text(size = 12,
                                        colour = "dodgerblue4"),
                        axis.title.x=element_blank(),
                        axis.text.x = element_text(face = "bold",
                                                    color="black"),
                        plot.title = element_text(hjust = 0.5,
                                        color = "black",
                                        size=14),
                        legend.title = element_blank()) +
                    labs(title = paste(protein, "Classification", sep = " "),
                            y = "Probability",
                            caption = sprintf("#PSM: %s", psm)) +
                    scale_x_discrete(labels=c("Compartment" = c.l,
                                                "Neighborhood" = n.l))

    plot(proteinPlot)
}

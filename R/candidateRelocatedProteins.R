#'@title Identify candidate relocated proteins
#'@description Identify candidate condition-dependent
#'relocated proteins by comparing neighborhood classifications
#'with respect to protein-protein pearson correlation and minumum PSM,
#'peptide spectrum matching, count.
#'@param sampleCls1 data.frame; merged classification,
#'combination of compartment and neighborhood classification.
#'@param s1PSM data.frame; minimum PSM count table
#'across ten TMT channel
#'@param s1Quant data.frame; fractionation quantification data
#'@param sampleCls2 data.frame; merged classification,
#'combination of compartment and neighborhood classification.
#'@param s2PSM data.frame; minimum PSM count table
#'across ten TMT channel
#'@param s2Quant data.frame; fractionation quantification data
#'@param annotation boolean; labeling the selected proteins
#'@param min.psm numeric; minimum psm, peptide spectra matching value
#'@param pearson.cor numeric; pearson correlation threshold
#'@export
#'@examples {
#'
#'candidate.df <- candidateRelocatedProteins(hcc827GEFClass, hcc827GefPSMCount,
#'hcc827GEF, hcc827GEFClass, hcc827GefPSMCount, hcc827GEF,
#'annotation = FALSE)
#'
#'
#'}
#'@import ggplot2
#'@import ggrepel
#'@importFrom graphics plot
#'@importFrom stats complete.cases
#'@return candidate.df

candidateRelocatedProteins <- function(sampleCls1, s1PSM,s1Quant, sampleCls2,
                                s2PSM, s2Quant, annotation = FALSE, min.psm,
                                pearson.cor){

    # just keep neighborhood classification
    df1 <- sampleCls1[, seq_len(2)]
    df2 <- sampleCls2[, seq_len(2)]

    df1 <- df1[intersect(rownames(df1), rownames(df2)), ]
    df2 <- df2[intersect(rownames(df1), rownames(df2)), ]

    df2 <- df2[rownames(df1), ]

    df <- data.frame(Protein = df1$Protein,
                        C.A = df1$NeighborhoodCls,
                        C.B = df2$NeighborhoodCls)

    ###########
    s1Quant <- SubCellBarCode::calRowMean(s1Quant)
    s2Quant <- SubCellBarCode::calRowMean(s2Quant)

    #remove unclassfied proteins at both conditions
    f.df <- df[!grepl("Unclassified", df$C.A) & ! grepl("Unclassified",df$C.B),]
    rownames(f.df) <- f.df$Protein

    #protein pearson correlation between samples
    f.df$Pearson.Corr <- unlist(unname((lapply(rownames(f.df), function(x){
        condA.prot <- s1Quant[x,]
        condB.prot <- s2Quant[x,]
        pearson.cor <- cor(t(condA.prot), t(condB.prot), method = "pearson")
    }))))


    #PSMCount
    df1PSMCount <- s1PSM[rownames(f.df), ]
    df2PSMCount <- s2PSM[rownames(f.df), ]

    #Minumum PSM count between two condition
    psm.df <- data.frame(Condition1 = df1PSMCount$PSMs.for.quant,
                        Condition2 = df2PSMCount$PSMs.for.quant)
    f.df$Min.PSMs <- apply(psm.df, 1, min)

    f.df$Relocated <- ifelse(f.df$C.A == f.df$C.B, "Background", "Relocated")
    f.df$ColorCode <- ifelse(f.df$Relocated == "Relocated",
                                "steelblue1", "Black")
    f.df <- f.df[order(f.df$Relocated , decreasing = FALSE), ]

    Relocalization <- f.df$Relocated

    if( ! annotation){

        plot(ggplot(f.df,
                aes(x = log2(f.df$Min.PSMs),
                    y = f.df$Pearson.Corr,
                    color = Relocalization))+
            geom_point(size = 2) +
            scale_colour_manual(values=c("Black","steelblue1")) +
            geom_hline(yintercept = 0.8, linetype="dashed", color = "red") +
            geom_vline(xintercept = 1.3, linetype="dashed", color = "red") +
            theme_bw() +
            theme(text = element_text(size = 16),
                    axis.text.x = element_text(face = "bold", color="black"),
                    axis.text.y = element_text(face = "bold", color="black")) +
            labs(title = "",
                    y = "Pearson Correlation ",
                    x = "Log2(Min.PSM)"))

        f.df <- f.df[f.df$Relocated == "Relocated", ]

        candidate.df <- f.df[,c(seq_len(5))]
        candidate.df <- candidate.df[candidate.df$Pearson.Corr < 0.8 &
                                        candidate.df$Min.PSMs > 2, ]

    }else{

        text.label.df <- subset(f.df, f.df$Pearson.Corr < pearson.cor &
                                    f.df$Min.PSMs > min.psm &
                                    f.df$Relocated == "Relocated")
        text.label.df$Label <- rownames(text.label.df)

        rem.df <- subset(f.df, !rownames(f.df) %in% rownames(text.label.df))
        rem.df$Label <- ""

        annot.df <- rbind(text.label.df, rem.df)
        annot.df <- annot.df[rownames(f.df),]

        Relocalization <- annot.df$Relocated
        plot(ggplot(annot.df,
                    aes(log2(annot.df$Min.PSMs),
                        annot.df$Pearson.Corr,
                        label = annot.df$Label,
                        color = Relocalization)) +
                geom_point(size = 2) +
                geom_text_repel() +
                scale_colour_manual(values=c("Black","steelblue1")) +
                geom_hline(yintercept = 0.8, linetype="dashed", color = "red") +
                geom_vline(xintercept = 1.3, linetype="dashed", color = "red") +
                theme_bw() +
                theme(text = element_text(size = 16),
                    axis.text.x = element_text(face = "bold", color="black"),
                    axis.text.y = element_text(face = "bold", color="black")) +
                labs(title = "",
                    y = "Pearson Correlation ",
                    x = "Log2(Min.PSM)"))

        annot.df <- annot.df[annot.df$Relocated == "Relocated", ]

        candidate.df <- annot.df[, c(seq_len(5))]
        candidate.df <- candidate.df[candidate.df$Pearson.Corr < 0.8 &
                                        candidate.df$Min.PSMs > 2, ]

    }
}



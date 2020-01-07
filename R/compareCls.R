#'@title Compare exon and gene centric classifications
#'@description Comparison of the gene centric and exon
#'centric classification. Additionally, correlation analysis
#'is performed using quantification data.
#'@param geneCls, data frame gene centric classification output
#'@param exonCls, data frame exon centric classification output
#'@export
#'@examples {
#'
#'exon.cls <- data.frame(Protein = c("ENSE00000331854",
#'                                      "ENSE00000331855",
#'                                      "ENSE00000331859"),
#'                          NeighborhoodCls = c("Cytosol",
#'                                       "Cytosol",
#'                                       "Cytosol"),
#'                          CompartmentCls = c("C1","C1","C1"),
#'                          Secretory = c(0.1, 0.1, 0.1),
#'                          Nuclear = c(0.2, 0.2, 0.2),
#'                          Cytosol = c(0.2, 0.2, 0.2),
#'                          Mitochondria = c(0.2, 0.2, 0.2),
#'                          S1 = c(0.2, 0.2, 0.2),
#'                          S2 = c(0.2, 0.2, 0.2),
#'                          S3 = c(0.2, 0.2, 0.2),
#'                          S4 = c(0.2, 0.2, 0.2),
#'                          N1 = c(0.2, 0.2, 0.2),
#'                          N2 = c(0.2, 0.2, 0.2),
#'                          N3 = c(0.2, 0.2, 0.2),
#'                          N4 = c(0.2, 0.2, 0.2),
#'                          C1 = c(0.2, 0.2, 0.2),
#'                          C2 = c(0.2, 0.2, 0.2),
#'                          C3 = c(0.2, 0.2, 0.2),
#'                          C4 = c(0.2, 0.2, 0.2),
#'                          C5 = c(0.2, 0.2, 0.2),
#'                          M1 = c(0.2, 0.2, 0.2),
#'                          M2 = c(0.2, 0.2, 0.2),
#'                          GeneSymbol = c("COPB1", "COPB1", "COPB1"),
#'                          PeptideCount = c(2, 4, 7))
#'
#'gene.cls <- data.frame(Protein = c("COPB1"),
#'NeighborhoodCls = c("Cytosol"),
#'CompartmentCls = c("C1"),
#'Secretory = c(0.1),
#'Nuclear = c(0.2),
#'Cytosol = c(0.2),
#'Mitochondria = c(0.2),
#'S1 = c(0.2),
#'S2 = c(0.2),
#'S3 = c(0.2),
#'S4 = c(0.2),
#'N1 = c(0.2),
#'N2 = c(0.2),
#'N3 = c(0.2),
#'N4 = c(0.2),
#'C1 = c(0.2),
#'C2 = c(0.2),
#'C3 = c(0.2),
#'C4 = c(0.2),
#'C5 = c(0.2),
#'M1 = c(0.2),
#'M2 = c(0.2))
#'
#'comp.df <- compareCls(gene.cls, exon.cls)
#'
#'}
#'@return c.df

compareCls <- function(geneCls, exonCls){

    #overlap gene symbols between gene and exon centric
    os <- intersect(unique(exonCls$GeneSymbol), geneCls$Protein)

    # get the common symbols
    geneCls <- subset(geneCls, geneCls$Protein %in% os)
    geneCls <- geneCls[, seq_len(3)]

    # get the common symbols
    exonCls <- subset(exonCls, exonCls$GeneSymbol %in% os)
    exonCls <- exonCls[,c(1, 2, 3, 23, 24)]

    # merge classification between gene centric and exon centric
    comp.cls <- lapply(geneCls$Protein, function(g){
        gene.df <- subset(geneCls, geneCls$Protein %in% g)
        exon.df <- subset(exonCls, exonCls$GeneSymbol %in% g)
        c.df <- cbind(gene.df, exon.df, row.names = NULL)
    })

    c.df <- replaced.df <- do.call("rbind", comp.cls)

    #reorder the columns and rename
    c.df <- c.df[,c(4, 7, 2, 5, 3, 6, 8)]
    colnames(c.df) <- c("Exon_id","GeneSymbol",
                            "Gene_Neighborhood", "Exon_Neighborhood",
                            "Gene_Compartment", "Exon_Compartment",
                            "PeptideCount")

    #compute pearson correlation
    c.df$Pearson.Corr <- unlist(unname((lapply(rownames(c.df), function(i){
        e <- SubCellBarCode::hcc827exon[as.character(c.df[i,1]),][2:11]
        g <- SubCellBarCode::hcc827Ctrl[as.character(c.df[i,2]),][1:10]
        pear.cor <- cor(t(e), t(g), method = "pearson")
    }))))

    return(c.df)
}

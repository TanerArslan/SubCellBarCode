#'@title Visualization of marker proteins by t-SNE map
#'@description The marker proteins are visualized in 3D t-SNE map
#'to see the distributions of the marker proteins.
#'@param markerProteins character; robust marker proteins, gene
#'symbols, that are present in the given data and overlapped with package's
#'marker protein list.
#'@param protein.data data.frame; fractionated proteomics data
#'@param dims integer; dimensionality
#'@param theta numeric; Speed/accuracy trade-off ,increase for less accuracy
#'@param perplexity integer; Perplexity parameter
#'@export
#'@examples {
#'
#'df <- loadData(SubCellBarCode::hcc827Ctrl)
#'
#'c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])
#'
#'set.seed(21)
#'tsneMap.df <- tsneVisualization(protein.data = df,
#'markerProteins = c.prots[1:20],
#'dims = 2, theta = c(0.4), perplexity = c(5))
#'}
#'
#'@import Rtsne
#'@import scatterplot3d
#'@importFrom graphics legend
#'@return tsneMap.df

tsneVisualization <- function(protein.data,
                                markerProteins, dims,
                                theta,
                                perplexity){

    tsne.df <- protein.data[markerProteins, ]
    annotation.df <- SubCellBarCode::markerProteins[rownames(tsne.df), ]

    if( ! identical(rownames(tsne.df), rownames(annotation.df)) )
        stop('Make sure your markerProteins are covered in 3365 Proteins and
            there is no extra entry.')

    if(dims == 3){

        message("Optimization process started. This may take some time!")

        lst1 <- unlist(lapply(theta, function(x)
                lapply(perplexity, function(y) sort(x + y))))
        lst2 <- unlist(lapply(theta, function(x)
                lapply(perplexity, function(y){
                tsne.mpa <- Rtsne::Rtsne(tsne.df, dims=3, theta=x, perplexity=y)
                tsne.mpa$itercosts[1]})))

        lst.df <- data.frame(Parameters = lst1, Cost = lst2)
        lst.df <- lst.df[order(lst.df$Cost, decreasing = FALSE), ]


        min.theta.perp <- unlist(strsplit(as.character(lst.df[1, 1]),
                                split = ".", fixed = TRUE))

        message("Optimization was performed.")

        #get the optimization parameters
        theta.val <- as.numeric(paste(0,
                        as.character(min.theta.perp[2]), sep = "."))

        perplexity.val <- as.numeric(min.theta.perp[1])

        cat("Theta value: ", theta.val)
        cat("\nPerplexity value: ", perplexity.val)


        rtsne.map <- Rtsne::Rtsne(tsne.df,
                        dims = 3,
                        theta= theta.val,
                        perplexity = perplexity.val)


        #plot 3D tsne-map
        d <- data.frame(x=rtsne.map$Y[, 1],
                        y=rtsne.map$Y[, 2],
                        z=rtsne.map$Y[, 3])

        scatterplot3d::scatterplot3d(d$y, d$x, d$z,
                                    pch=20,
                                    color=annotation.df$Colour,
                                    grid = TRUE,
                                    box = FALSE)
        legend("topright",
                legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
                            "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
                pch = 16,
                cex = 0.75,
                col = c("gold", "orange", "salmon", "tomato2",
                        "grey90","grey70", "grey50", "grey30",
                        "lightblue", "aquamarine", "cyan", "deepskyblue2",
                        "turquoise3", "burlywood4", "tan4"))

        scatterplot3d::scatterplot3d(d$x, d$y, d$z,
                                    pch=20,
                                    color=annotation.df$Colour,
                                    grid = TRUE,
                                    box = FALSE)
        legend("topright",
                legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
                            "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
                pch = 16,
                cex = 0.75,
                col = c("gold", "orange", "salmon", "tomato2",
                        "grey90", "grey70", "grey50", "grey30",
                        "lightblue","aquamarine", "cyan","deepskyblue2",
                        "turquoise3","burlywood4","tan4"))

        scatterplot3d::scatterplot3d(d$z, d$y, d$x,
                                    pch=20,
                                    color = annotation.df$Colour,
                                    grid = TRUE,
                                    box = FALSE)
        legend("topright",
                legend = c("S1", "S2", "S3", "S4","N1", "N2", "N3", "N4",
                            "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
                pch = 16,
                cex = 0.75,
                col = c("gold", "orange", "salmon", "tomato2",
                        "grey90", "grey70", "grey50", "grey30",
                        "lightblue", "aquamarine", "cyan", "deepskyblue2",
                        "turquoise3", "burlywood4", "tan4"))

        tsneMap.df <- data.frame(Proteins = rownames(tsne.df), d[, seq_len(3)])

    }else if(dims == 2){

        message("Optimization process started. This may take some time!")

        lst1 <- unlist(lapply(theta,
                function(x) lapply(perplexity, function(y) sort(x + y))))
        lst2 <- unlist(lapply(theta,
                function(x) lapply(perplexity, function(y){
            tsne.mpa <- Rtsne::Rtsne(tsne.df, dims = 2, theta=x, perplexity=y)
            tsne.mpa$itercosts[1]})))

        lst.df <- data.frame(Parameters = lst1, Cost = lst2)
        lst.df <- lst.df[order(lst.df$Cost, decreasing = FALSE), ]


        min.theta.perp <- unlist(strsplit(as.character(lst.df[1, 1]),
                                split = ".", fixed = TRUE))

        message("Optimization was performed.")

        #get the optimization parameters
        theta.val <- as.numeric(paste(0,
                        as.character(min.theta.perp[2]), sep = "."))

        perplexity.val <- as.numeric(min.theta.perp[1])

        cat("Theta value: ", theta.val)
        cat("\nPerplexity value: ", perplexity.val)
        rtsne.map <- Rtsne::Rtsne(tsne.df,
                        dims = 2,
                        theta = theta.val,
                        perplexity = perplexity.val)

        #plot 2D tsne-map
        d <- data.frame(x=rtsne.map$Y[, 1],
                        y=rtsne.map$Y[, 2])

        plot(d$x, d$y, col = annotation.df$Colour, type = "p" , pch = 20)
        legend("topright",
                inset=c(-0.04,0),
                legend = c("S1", "S2", "S3", "S4", "N1", "N2", "N3", "N4",
                            "C1", "C2", "C3", "C4", "C5", "M1", "M2"),
                pch = 16,
                cex = 0.75,
                xpd = TRUE,
                col = c("gold", "orange", "salmon", "tomato2",
                        "grey90", "grey70", "grey50", "grey30",
                        "lightblue", "aquamarine", "cyan", "deepskyblue2",
                        "turquoise3", "burlywood4", "tan4"))

        tsneMap.df <- data.frame(Proteins = rownames(tsne.df), d[, seq_len(2)])

    } else{

        stop("Plesae select the correct dimension. It is either 2 or 3.")

    }
}

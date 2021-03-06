% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeProbability.R
\name{mergeProbability}
\alias{mergeProbability}
\title{Merge compartment probabilities to neighborhood probabilities}
\usage{
mergeProbability(df)
}
\arguments{
\item{df}{data.frame; all predictions at the neighborhood level and
probablity vectors for each protein}
}
\value{
merged.df
}
\description{
Compartment levels classifications are summed up to
associated neighborhood levels. It is a helper function.
}
\examples{
{

#create mock data
df <- data.frame(Protein = "TP53",
S1 = as.numeric(0.02),
S2 = as.numeric(0.02),
S3 = as.numeric(0.02),
S4 = as.numeric(0.02),
N1 = as.numeric(0.72),
N2 = as.numeric(0.02),
N3 = as.numeric(0.02),
N4 = as.numeric(0.02),
C1 = as.numeric(0.02),
C2 = as.numeric(0.02),
C3 = as.numeric(0.02),
C4 = as.numeric(0.02),
C5 = as.numeric(0.02),
M1 = as.numeric(0.02),
M2 = as.numeric(0.02))

rownames(df) <- "TP53"

merged.df <- mergeProbability(df)

}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/applyThresholdToCompartment.R
\name{applyThresholdCompartment}
\alias{applyThresholdCompartment}
\title{Apply thresholds to compartments}
\usage{
applyThresholdCompartment(all.repA, all.repB, threshold.df)
}
\arguments{
\item{all.repA}{data.frame; all predictions and probablity vectors
for each protein in replicate A}

\item{all.repB}{data.frame; all predictions and probablity vectors
for each protein in replicate B}

\item{threshold.df}{data.frame; collection od precision and recall
values for each compaartment}
}
\value{
c.cls.df
}
\description{
Apply thresholds for all predictions to increase the
true positive rate and remove poor classification.
}
\examples{
{

df <- loadData(SubCellBarCode::hcc827Ctrl)

c.prots <- calculateCoveredProtein(rownames(df), markerProteins[,1])

set.seed(7)
c.prots <- sample(c.prots, 550)
cls <- svmClassification(c.prots, df, markerProteins)

test.A <- cls[[1]]$svm.test.prob.out
test.B <- cls[[2]]$svm.test.prob.out

t.c.df <- computeThresholdCompartment(test.A, test.B)

all.A <- cls[[1]]$all.prot.pred
all.B <- cls[[2]]$all.prot.pred

c.cls.df <- applyThresholdCompartment(all.A, all.B, t.c.df)
}
}

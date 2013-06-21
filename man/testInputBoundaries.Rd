\name{testInputBoundaries}
\alias{testInputBoundaries}
\docType{data}
\title{
A table of chromosome boundary positions for DNA copy number analysis
%%   ~~ data name/kind ... ~~
}
\description{
The entire length of the genome was divided into 50009 bins, with chromosomes
laid out in the usual genomic order: 1,...,22,X,Y. Each observation in the
table provides the \code{start} and \code{end} bin numbers of each chromosome
(\code{chrom}). 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(testInputBoundaries)}
\format{
  A data frame with 24 observations on the following 3 variables.
  \describe{
    \item{\code{chrom}}{a numeric vector}
    \item{\code{start}}{a numeric vector}
    \item{\code{end}}{a numeric vector}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Navin N, Kendall J, Troge J, Andrews P, Rodgers L, McIndoo J, Cook K,
Stepansky A, Levy D, Esposito D et al. 2011.
Tumour evolution inferred by single-cell sequencing. Nature 472(7341): 90-U119.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(testInputBoundaries)
## maybe str(testInputBoundaries) ; plot(testInputBoundaries) ...
}
\keyword{datasets}

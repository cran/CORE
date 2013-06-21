\name{testInputCORE}
\alias{testInputCORE}
\docType{data}
\title{
A table of DNA copy number gain events observed in 100 individual tumor cells
%%   ~~ data name/kind ... ~~
}
\description{
Each observation in the table corresponds to a DNA copy number gain event in
one of 100 individual breast cancer cells. The entire length of the genome was
divided into 50009 bins. An event is an interval in chromosome \code{chrom}
whose start and end bin numbers are given by \code{start} and \code{end}.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(testInputCORE)}
\format{
  A data frame with 2490 observations on the following 3 variables.
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
data(testInputCORE)
## maybe str(testInputCORE) ; plot(testInputCORE) ...
}
\keyword{datasets}

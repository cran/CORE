\name{CORE}
\alias{CORE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CORE(dataIn, keep = NULL, startcol = "start", endcol = "end", chromcol = "chrom", weightcol = "weight", maxmark = 1, minscore = 0, pow = 1, assoc = c("I", "J", "P"), nshuffle = 0, boundaries = NULL, seedme = sample(1e+08, 1), shufflemethod = c("SIMPLE", "RESCALE"), tiny = -1, distrib = c("vanilla", "multicore", "sge"), njobs = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataIn}{
%%     ~~Describe \code{dataIn} here~~
}
  \item{keep}{
%%     ~~Describe \code{keep} here~~
}
  \item{startcol}{
%%     ~~Describe \code{startcol} here~~
}
  \item{endcol}{
%%     ~~Describe \code{endcol} here~~
}
  \item{chromcol}{
%%     ~~Describe \code{chromcol} here~~
}
  \item{weightcol}{
%%     ~~Describe \code{weightcol} here~~
}
  \item{maxmark}{
%%     ~~Describe \code{maxmark} here~~
}
  \item{minscore}{
%%     ~~Describe \code{minscore} here~~
}
  \item{pow}{
%%     ~~Describe \code{pow} here~~
}
  \item{assoc}{
%%     ~~Describe \code{assoc} here~~
}
  \item{nshuffle}{
%%     ~~Describe \code{nshuffle} here~~
}
  \item{boundaries}{
%%     ~~Describe \code{boundaries} here~~
}
  \item{seedme}{
%%     ~~Describe \code{seedme} here~~
}
  \item{shufflemethod}{
%%     ~~Describe \code{shufflemethod} here~~
}
  \item{tiny}{
%%     ~~Describe \code{tiny} here~~
}
  \item{distrib}{
%%     ~~Describe \code{distrib} here~~
}
  \item{njobs}{
%%     ~~Describe \code{njobs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataIn, keep = NULL, startcol = "start", endcol = "end", 
    chromcol = "chrom", weightcol = "weight", maxmark = 1, minscore = 0, 
    pow = 1, assoc = c("I", "J", "P"), nshuffle = 0, boundaries = NULL, 
    seedme = sample(1e+08, 1), shufflemethod = c("SIMPLE", "RESCALE"), 
    tiny = -1, distrib = c("vanilla", "multicore", "sge"), njobs = 1) 
{
    doshuffles <- "NO"
    shufflemethod <- match.arg(shufflemethod)
    assoc <- match.arg(assoc)
    if (!(class(dataIn) \%in\% c("CORE", "matrix", "data.frame"))) 
        stop("invalid class of argument dataIn")
    if (class(dataIn) == "CORE") 
        if ((!("assoc" \%in\% keep) & dataIn$assoc != assoc) | 
            !(("pow" \%in\% keep) & dataIn$pow != pow)) {
            cankeep <- c("maxmark", "minscore", "nshuffle", "boundaries", 
                "seedme", "shufflemethod", "tiny")
            for (item in intersect(keep, cankeep)) assign(item, 
                dataIn[[item]])
            dataIn <- dataIn$input
            startcol <- "start"
            endcol <- "end"
            chromcol <- "chrom"
            weightcol <- "weight"
        }
    if (class(dataIn) != "CORE") {
        if (!(startcol \%in\% dimnames(dataIn)[[2]])) 
            stop("start column missing in input data")
        if (!(endcol \%in\% dimnames(dataIn)[[2]])) 
            stop("end column missing in input data")
        if (!(chromcol \%in\% dimnames(dataIn)[[2]])) {
            dataIn <- cbind(dataIn, rep(1, nrow(dataIn)))
            dimnames(dataIn)[[2]][ncol(dataIn)] <- chromcol
        }
        if (!(weightcol \%in\% dimnames(dataIn)[[2]])) {
            dataIn <- cbind(dataIn, rep(1, nrow(dataIn)))
            dimnames(dataIn)[[2]][ncol(dataIn)] <- weightcol
        }
        z <- as.matrix(dataIn[, c(chromcol, startcol, endcol, 
            weightcol), drop = F])
        rm(dataIn)
        dimnames(z)[[2]] <- c("chrom", "start", "end", "weight")
        if (!is.null(boundaries)) {
            if (!(startcol \%in\% dimnames(boundaries)[[2]])) 
                stop("start column missing in boundary table")
            if (!(endcol \%in\% dimnames(boundaries)[[2]])) 
                stop("end column missing in boundary table")
            if (!(chromcol \%in\% dimnames(boundaries)[[2]])) 
                stop("chrom column missing in boundary table")
            boundaries <- as.matrix(boundaries[, c(chromcol, 
                startcol, endcol), drop = F])
        }
        result <- switch(assoc, I = ICOREiteration(z, maxmark, 
            pow, minscore), J = JCOREiteration(z, maxmark, pow, 
            minscore), P = PCOREiteration(z, maxmark, minscore))
        returnme <- list(input = z, call = match.call(), minscore = minscore, 
            maxmark = maxmark, pow = pow, assoc = assoc, coreTable = result, 
            seedme = seedme, boundaries = boundaries, shufflemethod = shufflemethod, 
            nshuffle = nshuffle, tiny = tiny)
        if (nshuffle > 0) 
            doshuffles <- "FROMSCRATCH"
    }
    else {
        assoc <- dataIn$assoc
        pow <- dataIn$pow
        cankeep <- c("maxmark", "minscore", "nshuffle", "boundaries", 
            "seedme", "shufflemethod", "tiny")
        for (item in intersect(keep, cankeep)) assign(item, dataIn[[item]])
        returnme <- list(input = dataIn$input, call = match.call(), 
            minscore = minscore, maxmark = maxmark, pow = dataIn$pow, 
            assoc = dataIn$assoc, shufflemethod = shufflemethod, 
            seedme = seedme, nshuffle = nshuffle, tiny = tiny, 
            boundaries = boundaries)
        if (dataIn$coreTable[nrow(dataIn$coreTable), "score"] >= 
            minscore & nrow(dataIn$coreTable) < maxmark) {
            z <- dataIn$input
            if (assoc == "I") {
                for (i in 1:nrow(dataIn$coreTable)) {
                  fixus <- z[, "chrom"] == dataIn$coreTable[i, 
                    "chrom"] & z[, "start"] <= dataIn$coreTable[i, 
                    "start"] & z[, "end"] >= dataIn$coreTable[i, 
                    "end"]
                  z[fixus, "weight"] <- z[fixus, "weight"] * 
                    (1 - ((dataIn$coreTable[i, "end"] - dataIn$coreTable[i, 
                      "start"] + 1)/(z[fixus, "end"] - z[fixus, 
                      "start"] + 1))^dataIn$pow)
                }
                z <- z[z[, "weight"] >= tiny, , drop = F]
                result <- ICOREiteration(z, maxmark - nrow(dataIn$coreTable), 
                  dataIn$pow, minscore)
            }
            if (assoc == "J") {
                for (i in 1:nrow(dataIn$coreTable)) {
                  fixus <- z[, "chrom"] == dataIn$coreTable[i, 
                    "chrom"] & pmax(z[, "start"], dataIn$coreTable[i, 
                    "start"]) <= pmin(z[, "end"], dataIn$coreTable[i, 
                    "end"])
                  z[fixus, "weight"] <- z[fixus, "weight"] * 
                    (1 - ((pmin(dataIn$coreTable[i, "end"], z[fixus, 
                      "end"]) - pmax(dataIn$coreTable[i, "start"], 
                      z[fixus, "start"]) + 1)/(pmax(dataIn$coreTable[i, 
                      "end"], z[fixus, "end"]) - pmin(dataIn$coreTable[i, 
                      "start"], z[fixus, "start"]) + 1))^dataIn$pow)
                }
                z <- z[z[, "weight"] >= tiny, , drop = F]
                result <- JCOREiteration(z, maxmark - nrow(dataIn$coreTable), 
                  dataIn$pow, minscore)
            }
            if (assoc == "P") {
                for (i in 1:nrow(dataIn$coreTable)) {
                  fixus <- z[, "chrom"] == dataIn$coreTable[i, 
                    "chrom"] & z[, "start"] <= dataIn$coreTable[i, 
                    "start"] & z[, "end"] >= dataIn$coreTable[i, 
                    "end"]
                  z <- z[-fixus, , drop = F]
                }
                result <- PCOREiteration(z, maxmark - nrow(dataIn$coreTable), 
                  minscore)
            }
            returnme$coreTable <- rbind(dataIn$coreTable, result)
        }
        else returnme$coreTable <- dataIn$coreTable[which(dataIn$coreTable[1:min(maxmark, 
            nrow(dataIn$coreTable)), "score"] >= minscore), , 
            drop = F]
        if (nshuffle > 0) {
            if (dataIn$shufflemethod != shufflemethod | dataIn$seedme != 
                seedme | nrow(dataIn$coreTable) < nrow(returnme$coreTable) | 
                dataIn$nshuffle == 0) 
                doshuffles <- "FROMSCRATCH"
            else if (dataIn$nshuffle < nshuffle) 
                doshuffles <- "ADD"
            else returnme$simscores <- dataIn$simscores[1:nrow(returnme$coreTable), 
                1:nshuffle, drop = F]
        }
    }
    class(returnme) <- "CORE"
    if (is.null(boundaries)) {
        y <- returnme$input[order(returnme$input[, "chrom"]), 
            , drop = F]
        boundaries <- cbind(unique(y[, "chrom"]), tapply(X = y[, 
            "start"], INDEX = y[, "chrom"], FUN = min), tapply(X = y[, 
            "end"], INDEX = y[, "chrom"], FUN = max))
    }
    dimnames(boundaries)[[2]] <- c("chrom", "start", "end")
    randfun <- match.fun(paste(assoc, "CORErandomized", sep = ""))
    distrib <- match.arg(distrib)
    exportus <- c("advanceRNG", "containment.indicator", "fill.values", 
        "makepairs", "randomEventMoves", "envelope.indicator", 
        "interval.schedule", "overlap.indicator", "randomRescaledEventMoves")
    if (distrib == "multicore" & doshuffles != "NO") {
        ncores <- min(njobs, switch(doshuffles == "FROMSCRATCH", 
            nshuffle, nshuffle - dataIn$nshuffle), detectCores())
        cl <- makeCluster(getOption("cl.cores", ncores))
        clusterExport(cl = cl, varlist = exportus)
    }
    if (doshuffles == "FROMSCRATCH") {
        returnme$simscores <- switch(distrib, vanilla = randfun(COREobj = returnme, 
            boundaries = boundaries), multicore = matrix(nrow = nrow(returnme$coreTable), 
            data = unlist(parSapply(cl = cl, X = 1:ncores, FUN = randfun, 
                COREobj = returnme, boundaries = boundaries, 
                nprocs = ncores))), sge = matrix(nrow = nrow(returnme$coreTable), 
            data = unlist(sge.parSapply(X = 1:njobs, FUN = randfun, 
                COREobj = returnme, boundaries = boundaries, 
                nprocs = njobs, njobs = njobs, global.savelist = exportus, 
                simplify = F))))
    }
    else if (doshuffles == "ADD") {
        returnme$simscores <- cbind(dataIn$simscores[1:nrow(returnme$coreTable), 
            , drop = F], switch(distrib, vanilla = randfun(COREobj = returnme, 
            boundaries = boundaries, rngoffset = dataIn$nshuffle), 
            multicore = matrix(nrow = nrow(returnme$coreTable), 
                data = unlist(parSapply(cl = cl, X = 1:ncores, 
                  FUN = randfun, COREobj = returnme, boundaries = boundaries, 
                  rngoffset = dataIn$nshuffle, nprocs = ncores))), 
            sge = matrix(nrow = nrow(returnme$coreTable), data = unlist(sge.parSapply(X = 1:njobs, 
                FUN = randfun, COREobj = returnme, boundaries = boundaries, 
                rngoffset = dataIn$nshuffle, nprocs = njobs, 
                njobs = njobs, global.savelist = exportus, simplify = F)))))
    }
    if ("simscores" \%in\% names(returnme)) 
        returnme$p <- (rowSums(returnme$simscores > returnme$coreTable[, 
            "score"]) + 1)/(ncol(returnme$simscores) + 2)
    if (exists("cl")) 
        stopCluster(cl)
    return(returnme)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

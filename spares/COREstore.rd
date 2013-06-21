\name{CORE}
\alias{CORE}
\title{Cores of Recurrent Events}

\description{
Given a collection of intervals s_1,...,s_N, find K intervals c_1,...,c_K which 
approximately minimize Sum_i Prod_k (1-E(s_i,c_k)), where E(s_i,c_k) is a 
geometric measure of association between s_i and c_k. Perform permutation 
tests to estimate the significance of finding.
}
\usage{
CORE(dataIn, keep = NULL, startcol = "start", endcol = "end", 
chromcol = "chrom", weightcol = "weight", maxmark = 1, minscore = 0, pow = 1, 
assoc = c("I", "J", "P"), nshuffle = 0, boundaries = NULL, 
seedme = sample(1e+08, 1), shufflemethod = c("SIMPLE", "RESCALE"), tiny = -1, 
distrib = c("vanilla", "multicore", "sge"), njobs = 1)
}
\arguments{
\item{datain}{
A matrix, a data frame or an object of class "CORE". If \code{dataIn} 
is a matrix or a data frame, it should have columns with names specified by 
the \code{startcol} and \code{endcol} arguments, otherwise the function exits 
with an error.
}
\item{keep}{
A character vector. If \code{dataIn} is of class "CORE", \code{keep} specifies 
the names of items of \keep{dataIn} to be kept at their input values. These 
values take precedence over the corresponding argument values as specified in 
the function call. \code{keep} is ignored if \code{dataIn} is not of class 
"CORE".
}
\item{startcol}{
A character string. If \code{dataIn} is a matrix or a data frame, 
\code{startcol} specifies the name of the column containing start coordinates 
of the input intervals. Otherwise \code{startcol} is ignored.
}
\item{endcol}{
A character string. If \code{dataIn} is a matrix or a data frame, 
\code{endcol} specifies the name of the column containing end coordinates 
of the input intervals. Otherwise \code{endcol} is ignored.
}
\item{chromcol}{
A character string. If \code{dataIn} is a matrix or a data frame, 
\code{chromcol} specifies the name of the column containing chromosome numbers 
of the input intervals. Otherwise \code{chromcol} is ignored.
}
\item{weightcol}{
A character string. If \code{dataIn} is a matrix or a data frame, 
\code{weightcol} specifies the name of the column containing initial weights 
of the input intervals. Otherwise \code{weightcol} is ignored.
}
\item{maxmark}{
An integer for the maximal number of cores to be computed. The actual number 
of cores to be computed is the smaller of \code{maxmark} and the number of 
cores with scores exceeding \code{minscore}.
}
\item{minscore}
{
A single numeric value for the minimal allowed score of the cores to be 
reported.
}
\item{pow}{
A single numeric value of at least 1 for the power parameter used in 
computing the association measure beween the cores and the input intervals
(see Details).
}
\item{assoc}{
A character specifying the type of association measure to be used (see Details).
}
\item{nshuffle}{
An integer specifying the number of randomizations to be performed for 
estimating significance.
}
\item{boundaries}{
A matrix or a data frame that must have three columns whose names are given by 
\code{chromcol}, \code{startcol} and \code{endcol}. These specify the 
chromosome numbers and their start and end positions (see Details).
}
\item{seedme}{
An integer specifying the random number generator seed (see Details).
}
\item{shufflemethod}{
A character string specifying the event randomization method used for 
estimation of significance. If "SIMPLE" (default), each event is placed at 
random with equal probability for any position where it can fit within 
chromosome boundaries. If "RESCALE", each event is placed at random in a 
randomly chosen chromosome, and the event length is multiplied by the length 
ratio of the new to the original chromosome.
}
\item{tiny}{
A single numeric value specifying the weight below which events are removed 
from the input event set.
}
\item{distrib}{
A character string specifying the method of distributed computing used for 
estimation of significance. If "vanilla" (default), no distributed computing 
is performed. If "multicore", parallel computation is performed using 
functions from CRAN core package parallel, with the number of worker processes 
being the least of \code{njobs}, \code{nshuffle} and the value returned by the 
detectCores function of parallel. If "sge", parallel computation is performed 
using functions from the CRAN package Rsge, with the number of worker 
processes being the least of \code{njobs} and \code{nshuffle}.
}
\item{njobs}{
If distributed computing is used for estimation of significance, a single 
integer specifying the desired number of worker processes.
}
}
\details{
The three measures of association specified by \code{associ}re defined as 
follows (|| denotes the length of an interval). For "I" (inclusion) 
E(s_i,c_k) = (|c_k|/|s_i|)^pow if c_k is contained in s_i and 0 otherwise. 
For "J" (Jaccard) E(s_i,c_k) = J(s_i,c_k)^pow, where J is the Jaccard index. 
For "P" (piercing) E(s_i,c_k) = 1 if c_k is contained and 0 otherwise. 
In all cases the left (right) boundary of an optimal c_k is one of the left 
(right) boundaries in the set of input interval events. In addition, there 
are no event interval boundaries in the interior of an optimal c_k in case "P".

The \code{boundaries} argument is used for assessing statistical significance 
of the solution. If \code{boundaries} is not specified, the chromosome 
boundaries for each chromosome are taken to be the leftmost left and the 
rightmost right boundaries of all events in the chromosome.

If significance of finding is estimated, the random number generator stream, 
and hence the resultant estimate, only depends on \code{seedme} and is 
independent of the parallelization option chosen.
}
\value{
An object of class "CORE" with the following items.
\item{input}{
A matrix with four columns called "chrom", "start", "end" and "weight", 
specifying the input interval events.
}
\item{call}{
A character string specifying the function call.
}
\item{coreTable}{
A matrix with columns named "start", "end" and "score", for start and end 
positions and CORE scores of the cores found by the algorithm.
}
\item{seedme}{
If significance estimate was performed, the random number generator seed.
}
\item{assoc}{
One of "I", "J" or "P", indicating the geometric measure of association used.
}
\item{shufflemethod}{
One of "SIMPLE" or "RESCALE", indicating the randomization method used.
}
\item{p}{
A numeric vector of the length equal to the row dimension of \code{coreTable} 
containing estimated p-values for the cores.
}
\item{simscores}{
A matrix with the row dimension equal to that of \code{coreTable} and 
\code{nshuffle} columns, containing core scores computed for \code{nshuffle} 
sets of randomized events.
}
\item{minscore}{
A single numeric value for the minimal score of the reported cores.
}
\item{maxmark}{
A single numeric value for the requested maximal number of cores to be 
computed.
}
\item{tiny}{
A single numeric value for the weight below which events were removed from 
the input set.
}
\item{pow}{
A single numeric value for the power used in computing the association 
measures.
}
\item{boundaries}{
A matrix with three columns named "chrom", "start" and "end", indicating 
chromosome numbers and boundary positions used for estimation of significance.
}
}
\author{Alex Krasnitz}
\examples{
\dontrun{
#Compute 10 cores and perform a very small number of randomizations (not enough
#for a meaningful estmate of significance).
myCOREobj<-CORE(dataIn=testInputCORE,maxmark=10,nshuffle=3,
boundaries=testInputBoundaries,seedme=123)
#Extend this computation to a much larger number of randomizations, using 8 
#cores of a host computer.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
nshuffle=100,distrib="multicore",njobs=8)
}
}

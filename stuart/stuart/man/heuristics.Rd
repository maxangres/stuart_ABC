% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heuristics.R
\name{heuristics}
\alias{heuristics}
\title{Generating heuristics for the use in STUART subtest construction}
\usage{
heuristics(
  data,
  factor.structure,
  capacity = NULL,
  repeated.measures = NULL,
  mtmm = NULL,
  grouping = NULL,
  localization = "nodes",
  ...
)
}
\arguments{
\item{data}{A data.frame containing all relevant data.}

\item{factor.structure}{A list linking factors to items. The names of the list elements correspond to the factor names. Each list element must contain a character-vector of item names that are indicators of this factor.}

\item{capacity}{A list containing the number of items per subtest. This must be in the same order as the \code{factor.structure} provided. If a single number, it is applied to all subtests. If \code{NULL} all items are evenly distributed among the subtests.}

\item{repeated.measures}{A list linking factors that are repeated measures of each other. Repeated factors must be in one element of the list - other sets of factors in other elements of the list. When this is \code{NULL} (the default) a cross-sectional model is estimated.}

\item{mtmm}{A list linking factors that are measurements of the same construct with different methods. Measurements of the same construct must be in one element of the list - other sets of methods in other elements of the list. When this is \code{NULL} (the default) a single method model is estimated.}

\item{grouping}{The name of the grouping variable. The grouping variable must be part of \code{data} provided and must be a numeric variable.}

\item{localization}{Which parameterization to use when depositing pheromones. Can be either 'nodes' (the default) for depositing pheromones on selected nodes or 'arcs' for depositing on selection arcs.}

\item{...}{Other arguments normally provided to \code{\link{mmas}}, which will be ignored.}
}
\value{
Returns a list of the same length as the \code{factor.structure} argument provided.
}
\description{
Creates uninformative heuristic matrices for the use in \code{\link{mmas}}.
}
\details{
This function generates a list of matrices which can be used as heuristics for all STUART constructions. This is mainly intended to write the structure of the heuristic matrices to an object, change components in line with theoretically derived heuristics and feed them back into \code{\link{mmas}} via the \code{heuristics} argument. The generated heuristics will contain only 1s and 0s, making it no heuristic information. Selection probabilities can be altered by manipulating the contents of the object created by \code{heuristics}. Setting a value to 0 will result in prohibiting a certain choice to be made. Please note, that it will lead to unpredictable behavior if the diagonal elements of the matrices produced in the arcs parameterization are set to values other than 0.
}
\examples{

# heuristics for node localization
data(fairplayer)
fs <- list(si = names(fairplayer)[83:92])

(heu <- heuristics(fairplayer, fs, 4))

# Define anchor-item
heu$si[1] <- 10000
heu

# heuristics for arc localization
data(fairplayer)
fs <- list(si = names(fairplayer)[83:92])

(heu <- heuristics(fairplayer, fs, 4, localization = 'arcs'))

# Define equal selection of odd and even items
heu$si[1:10,] <- c(rep(c(0, 1), 5), rep(c(1, 0), 5))
heu

}
\seealso{
\code{\link{mmas}}
}
\author{
Martin Schultze
}
\concept{ACO subtests}

\name{wbgt_df}
\alias{wbgt_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the web bulb globe temperature (wbgt)}
\description{
    An R wrapper for Liljegren's C implemenation of his method for computing the wbgt.
}
\usage{
wbgt_df(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A data frame containing columns with names equal to the inputs into the \code{wbgt} function in this package.
  }
}
\value{
    A copy of the data frame input with a new column containing the wet bulb globe temperature in degrees C 
    for each observation. In contrast to the function `wbgt` in this package, this function does not 
    retain the additional calculations (`Tpsy`, `Tnwb`, `Tg`, `status`).
}
\author{
    Max Lieblich, June Spector
}

%\examples{
%wbgt(2015, 1, 2, 4, 4, -8, 15, 46.2, -120.0, 1000, 1017.8, 12, 45, 1.24, 2.1, -0.1, 0)
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}

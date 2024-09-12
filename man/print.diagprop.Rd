\name{print.diagprop}
\alias{print.diagprop}
\title{
Print diagprop object
}
\description{
Prints the p-values of testing proportionality of the subdistribution hazards for each covariate and the p-value of overall proportionality. 
}
\usage{
\method{print}{diagprop}(x, ...)
}

\arguments{
  \item{x}{a diagprop object}
  \item{...}{ additional arguments to print()}
}


\seealso{
\code{\link{diag_prop}}, \code{\link{plot.diagprop}}
}

\examples{
library(crskdiag)
data(dat2)
out2 <- diag_crr(Crsk(time,cause)~z1+z2,data=dat2,test="prop",seed=1234)
print(out2)
}

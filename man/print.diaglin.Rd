\name{print.diaglin}
\alias{print.diaglin}
\title{
Print diaglin object
}
\description{
Prints the p-values of testing linear functional form for each covariate and the p-value of testing link function. 
}
\usage{
\method{print}{diaglin}(x, ...)
}
\arguments{
  \item{x}{a diaglin object}
  \item{...}{ additional arguments to print()}
}


\seealso{
\code{\link{diag_lin}}, \code{\link{plot.diaglin}}
}

\examples{
library(crskdiag)
data(dat1)
out1 <- diag_crr(Crsk(time,cause)~z1+z2,data=dat1,test="lin",seed=1234)
print(out1)
}

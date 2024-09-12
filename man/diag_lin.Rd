\name{diag_lin}
\alias{diag_lin}
\title{
Diagnostics for the linear functional form and the link function 
}
\description{
Checks the linear functional form and the log-log link function assumptions of the Fine and Gray model. For internal use only.
}
\usage{
diag_lin(t, ic, z, n.total, Nit, n.sim, n.plot, seed, minor_included)
}
\arguments{
  \item{t}{a vector of observed times.}
  \item{ic}{a vector of cause indicators.}
  \item{z}{a matrix of covariates for fitting the model.}
  \item{n.total}{total number of samples}
  \item{Nit}{number of iterations for model estimation. Default is 20.}
  \item{n.sim}{number of iterations in the analytical diagnostics. Default is 1000.}
  \item{n.plot}{number of plots of simulated processes used in the graphical diagnostics. Default is 10.}
  \item{seed}{a seed used in generating the simulated processes.}
  \item{minor_included}{include the minor term in FG model if TRUE. Default is TRUE.}
}
\details{
More details can be found in Li, Scheike and Zhang (2015).
}
\value{
a "diaglin" class object.
}
\references{
Li, Jianing, Thomas H. Scheike, and Mei-Jie Zhang. "Checking Fine and Gray subdistribution hazards model with cumulative sums of residuals." 
Lifetime data analysis 21.2 (2015): 197-217.
}
\author{
Jianing Li
}

\seealso{
\code{\link{diag_crr}}, \code{\link{plot.diaglin}}, \code{\link{print.diaglin}}
}

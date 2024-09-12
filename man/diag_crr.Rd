\name{diag_crr}
\alias{diag_crr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking Fine and Gray subdistribution hazards model with cumulative sums of residuals
}
\description{
Provides with a class of analytical methods and graphical approaches for checking the assumptions of the Fine and Gray subdistribution hazards model based on the cumulative sums of residuals. It validates the model in three aspects: 
proportionality of hazard ratio, the linear functional form, and the link function.
}
\usage{
diag_crr(formula, data, test = c("lin", "prop"), Nit = 20, n.sim = 1000, n.plot = 10, 
		seed = NULL, minor_included = 1)
}
\arguments{
  \item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a competing risks object as returned by the Crsk function.}
  \item{data}{a dataset contains the time, causes of event, and the covariates. For the causes, 0 is censoring, 1 is the cause of interest, 2 is all other causes.}
  \item{test}{type of diagnostic method used in the test. One test at a time. "lin" generates the results of testing linear functional form for each covariate and link function. "prop" generates the results of testing proportionality for each covariate and overall proportionality.}
  \item{Nit}{number of iterations in model estimation. Default is 20.}
  \item{n.sim}{number of iterations in the analytical diagnostics. Default is 1000.}
  \item{n.plot}{number of plots of simulated processes used in the graphical diagnostics. Default is 10.}
  \item{seed}{a seed used in generating the simulated processes.}
  \item{minor_included}{include the minor term in FG model if TRUE. Default is TRUE.}
}
\details{
More details about the model diagnostic procedures can be found in Li, Scheike and Zhang (2015).
}
\value{
 \item{test}{type of diagnostic method used in the test.}
 \item{varname}{covariate names.}
 \item{used}{total sample size used in analysis. Subject with missing value will be excluded from analysis.}
 \item{NJP}{total number of unique time points for the cause of interest.}
 \item{TJP}{unique time points for the cause of interest leading by a zero.}
 \item{n.plot}{number of plots of simulated processes used in the graphical diagnostics.}
 \item{beta}{estimated regression coefficients.}
 \item{beta.se}{estimated standard errors for regression coefficients.}
 \item{dlamb0}{a vector of estimated cumulative baseline subdistribution hazards at observed failure times.}
 \item{pval}{p-values of the model diagnostic test.}
 \item{mav}{maximum value of the cumulative residual process (lin) or standardized residual process (prop).}
 \item{B}{cumulative residual process (lin) or standardized cumulative residual process (prop). Used for plot functions.}
 \item{uniX}{unique values of the covariates. Used for plot and available only when test="lin".}
 \item{TC}{number of unique covariates. Available only when test="lin".}
}
\references{
Li, Jianing, Thomas H. Scheike, and Mei-Jie Zhang. "Checking Fine and Gray subdistribution hazards model with cumulative sums of residuals." 
Lifetime data analysis 21.2 (2015): 197-217.

Fine, Jason P., and Robert J. Gray. "A proportional hazards model for the subdistribution of a competing risk." 
Journal of the American statistical association 94.446 (1999): 496-509.
}
\author{
Jianing Li
}

\seealso{
\code{\link{diag_lin}}, \code{\link{diag_prop}}, \code{\link{plot.diaglin}}, \code{\link{plot.diagprop}}, \code{\link{print.diaglin}}, \code{\link{print.diagprop}}
}
\examples{
data(dat1)
out1 <- diag_crr(Crsk(time,cause)~z1+z2,data=dat1,test="lin",seed=1234)
print(out1)
plot(out1)

data(dat2)
out2 <- diag_crr(Crsk(time,cause)~z1+z2,data=dat2,test="prop",seed=1234)
print(out2)
plot(out2)
}

\keyword{subdistribution}
\keyword{diagnostics}
\keyword{linear functional form}
\keyword{link function}
\keyword{proportionality}

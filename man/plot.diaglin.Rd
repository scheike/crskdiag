\name{plot.diaglin}
\alias{plot.diaglin}

\title{
Plot diaglin object
}

\description{
Plots the observed and simulated cumulative sums of residuals in testing linear functional form for each covariate and link function for the fitted Fine and Gray model.
}

\usage{
\method{plot}{diaglin}(x, col = c("red", "black"), lty = c(1, 2), lwd = c(2, 1), txt.pos=c(0.85,0.1), 
		lgd.pos="topright", xlim=NULL, ylim=NULL, select=NULL, ...)
}

\arguments{
  \item{x}{a diaglin object.}
  \item{col}{colors of the plots. The first element is used for observed process, and the second element is used for simulated processes.}
  \item{lty}{line types of the plots. The first element is used for observed process, and the second element is used for simulated processes.}
  \item{lwd}{line widths of the plots. The first element is used for observed process, and the second element is used for simulated processes.}
  \item{txt.pos}{the horizontal and vertical positions of the p-value. The assigned numbers represents the percentage in the x-y panel.}
  \item{lgd.pos}{the position of the legend.}
  \item{xlim}{the x limitation of the plots}
  \item{ylim}{the y limitation of the plots}
  \item{select}{the covariate needs to be plotted. Valid value is an integer from 1 to the number of covariates plus 1. 
		For example, if the number of covariates is 3, then 1, 2, 3 indicates the plots of checking corresponding covariates in the formula, and 4 indicates the plot of checking the link function. 
		Multiple values are allowed}
  \item{...}{ additional arguments to plot()}
}

\seealso{
\code{\link{diag_lin}}, \code{\link{print.diaglin}}
}

\examples{
library(crskdiag)
data(dat1)
out1 <- diag_crr(Crsk(time,cause)~z1+z2,data=dat1,test="lin",seed=1234)
plot(out1)
}

\name{n.wilcox.ord}
\alias{n.wilcox.ord}
\title{Sample size for Wilcoxon-Mann-Whitney for ordinal data}
\description{
  Function computes sample size for the two-sided Wilcoxon test when applied to two independent samples with ordered categorical responses.
}
\usage{
n.wilcox.ord(power = 0.8, alpha = 0.05, t, p, q)
}

\arguments{
  \item{power}{required Power}
  \item{alpha}{required two-sided Type-I-error level}
  \item{t}{sample size fraction n/N, where n is sample size of group B and N is the total sample size}
  \item{p}{vector of expected proportions of the categories in group A, should sum to 1}
  \item{q}{vector of expected proportions of the categories in group B, should be of equal length as p and should sum to 1}
}

\details{
This function approximates the total sample size, N, needed for the two-sided Wilcoxon test when comparing two independent samples, A and B, when data are ordered categorical according to Equation 12 in Zhao et al.(2008). Assuming that the response consists of D ordered categories \eqn{C_1 ,..., C_D}. The expected proportions of these categories in two treatments A and B must be specified as numeric vectors \eqn{p_1,...,p_D} and \eqn{q_1,...,q_D}, respectively. The argument t allows to compute power for an unbalanced design, where \eqn{t=n_B/N} is the proportion of sample size in treatment B.
}

\value{
  \item{total sample size}{Total sample size}
  \item{m}{Sample size group 1}
  \item{n}{Sample size group 2}
}
\references{Zhao YD, Rahardja D, Qu Yongming. Sample size calculation for the Wilcoxon-Mann-Whitney test adjsuting for ties. Statistics in Medicine 2008; 27:462-468}
\author{Ralph Scherer}
\examples{
## example out of:
## Zhao YD, Rahardja D, Qu Yongming. 
## Sample size calculation for the Wilcoxon-Mann-Whitney test adjsuting for ties. 
## Statistics in Medicine 2008; 27:462-468
n.wilcox.ord(power = 0.8, alpha = 0.05, t = 0.53, p = c(0.66, 0.15, 0.19), q = c(0.61, 0.23, 0.16))
}
\keyword{htest}
\keyword{nonparametric}
\concept{sample size}
\concept{ordered categorical}

\name{n.ttest}
\alias{n.ttest}
\title{
  n.ttest computes sample size for paired and unpaired t-tests.
}
\description{
  n.ttest computes sample size for paired and unpaired t-tests. Design
  may be balanced or unbalanced. Homogeneous and heterogeneous variances
  are allowed.
}
\usage{
n.ttest(power = 0.8, alpha = 0.05, mean.diff = 0.8, sd1 = 0.83, sd2 = 2.65, 
        k = 1, design = "unpaired", fraction = "balanced", variance = "equal")
}
\arguments{
  \item{power}{
    Power (1 - Type-II-error)
  }
  \item{alpha}{
    Two-sided Type-I-error}
  \item{mean.diff}{
    Expected mean difference
  }
  \item{sd1}{
    Standard deviation in group 1
  }
  \item{sd2}{
    Standard deviation in group 2
  }
  \item{k}{
    Sample fraction k
  }
  \item{design}{
    Type of design. May be paired or unpaired
  }
  \item{fraction}{
    Type of fraction. May be balanced or unbalanced
  }
  \item{variance}{
    Type of variance. May be homo- or heterogeneous
  }
}
\value{
  \item{Total sample size }{Sample size for both groups together}
  \item{Sample size group 1 }{Sample size in group 1}
  \item{Sample size group 2 }{Sample size in group 2}
}
\references{
  Bock J., Bestimmung des Stichprobenumfangs fuer biologische
  Experimente und kontrollierte klinische Studien. Oldenbourg 1998
}
\author{
  Ralph Scherer
}
\examples{
n.ttest(power = 0.8, alpha = 0.05, mean.diff = 0.80, sd1 = 0.83, k = 1,
design = "unpaired", fraction = "balanced", variance = "equal")

n.ttest(power = 0.8, alpha = 0.05, mean.diff = 0.80, sd1 = 0.83, sd2 =
2.65, k = 0.7, design = "unpaired", fraction = "unbalanced", variance =
"equal")
}
\keyword{ htest }

\name{meshr-package}
\alias{meshr-package}
\alias{meshr}
\docType{package}
\title{
Enrichment analysis for MeSH terms. 
}
\description{
  meshr package conducts a MeSH enrichment analysis employing gene-MeSH 
  annotation data. A hypergeometric test accounting for a multiple
  tesing correction is used to find significantly enriched MeSH terms. 
}
\details{
\tabular{ll}{
Package: \tab meshr\cr
Version: \tab 0.99.0\cr
Date: \tab 12-16-2012\cr
biocViews: \tab Bioinformatics, Statistics, Annotation,
MultipleComparisons\cr
Imports: \tab MeSH.db, fdrtool, Category, methods\cr
Depends: \tab R (>= 2.10), gendoo.Hs.db, MeSH.db, fdrtool, Category, methods\cr
License: \tab Artistic-2.0\cr
}

Index:
\preformatted{
meshHyperGTest        Performs a hypergeometric statistical test. 
}

Further information is available in the vignettes. 

}
\author{
Gota Morota, Koki Tsuyuzaki, Takeru Nakazato, Itoshi Nikaido

Maintainer: Gota Morota <morota@ansci.wisc.edu>
}

\examples{
ls("package:meshr")
}

\seealso{
  \code{\link{MeSHHyperGParams-class}},  \code{\link{MeSHHyperGResult-class}},
  \code{\link{meshHyperGTest}}
}

\keyword{ package }

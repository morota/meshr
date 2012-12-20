\name{meshHyperGTest}
\alias{meshHyperGTest}
\alias{meshHyperGTest,MeSHHyperGParams-method}

\title{Hypergeometric Tests for MeSH term association}

\description{ Given a \code{MeSHHyperGParams} object 
  containing a set of selected and background gene IDs,
  and  gene-MeSH annotation data of interest, \code{meshHyperGTest}
  performs Hypergeomtric test for over-representation of each MeSH term
  accouting for the multiple testing correction.  }

\arguments{
  \item{p}{A \code{MeSHHyperGParams} object }
}

\details{

  For details on creating \code{MeSHHyperGParams} object, please read
  the documentation in the \link[meshr]{MeSHHyperGParams-class}. 

}

\value{
  A \code{MeSHHyperGResult} object.
}


\author{
Gota Morota, Koki Tsuyuzaki, Takeru Nakazato, Itoshi Nikaido

Maintainer: Gota Morota <morota@ansci.wisc.edu>
}


\examples{
## create selected and universe gene vectors 
gene.table <- toTable(org.Hs.egALIAS2EG)
gene.list <- unique(gene.table[,1])

set.seed(139949)
n <- 100
t <- 999

# sample gene index 
genes.index  <- sample(1:length(gene.list),t)
# sample selected gene index 
selected.genes <- gene.list[genes.index[1:n]]
# sample universe gene index 
universe.genes <- gene.list[genes.index[1:t]]

## create MeSHHyperGParams instance  
meshParams <- new("MeSHHyperGParams",
                   geneIds=selected.genes,
                   universeGeneIds=universe.genes,
                   annotation="GendooMeSHA",
                   pvalueCutoff=0.05,
                   pAdjust="none")

## hypergeometric test 
meshR <- meshHyperGTest(meshParams) 
summary(meshR)

## hypergeometric test using different annotation data 
meshParams@annotation <- "GendooMeSHB" 
meshR <- meshHyperGTest(meshParams) 
summary(meshR)
}


\seealso{
  \code{\link{meshr-package}}, \code{\link{MeSHHyperGParams-class}},
  \code{\link{MeSHHyperGResult-class}}
}


\keyword{models}


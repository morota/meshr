##
## Accessor methods for MeSHHyperGResult class 
##

setMethod("meshIds", "MeSHHyperGResult", function(r) r@meshIds)

setMethod("meshTerms", "MeSHHyperGResult", function(r) r@meshTerms)

setMethod("meshCategory", "MeSHHyperGResult", function(r) r@meshCategory)

setMethod("annotation", "MeSHHyperGResult", function(object) object@annotation)

## generic is defined in Category/R/AllGenerics.R
setMethod("pvalues", "MeSHHyperGResult", function(r) r@pvalues)

## generic is defined in methods/R/AllGenerics.R
setMethod("show", "MeSHHyperGResult", function(object){
cat("MeSH enrichment analysis for category", object@meshCategory, '\n')
cat("Annotation used: ", object@annotation, '\n')
cat("Number of MeSH terms identified: ", length(object@meshIds), '\n')
})

## generic is defined in base/R/AllGenerics.R
setMethod("summary", "MeSHHyperGResult", function(object){
  mesh.df <- data.frame(MESHID = object@meshIds, MESHTERM = object@meshTerms, PVALUE = object@pvalues)
  mesh.df
})


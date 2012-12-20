setMethod("meshHyperGTest", signature(p="MeSHHyperGParams"),
          function(p) .meshHyperGTestInternal(p) )

.meshHyperGTestInternal <- function(p) {

  ##
  ## MeSH enrichment analysis 
  ##

  ## Map gene ID to MeSH ID through annotation data
  my.keytype <- c("gene_id")
  my.cols <- c("gene_id", "mesh")
  selected.mesh <- select(get(p@annotation), keys = p@geneIds, cols = my.cols, keytype = my.keytype)
  universe.mesh <- select(get(p@annotation), keys = p@universeGeneIds, cols = my.cols, keytype = my.keytype)
  selected.table <- table(selected.mesh[,2])
  universe.table <- table(universe.mesh[,2])

  ## Hypergeometric test 
  pvals <- array()
  for (i in 1:length(p@geneIds)) {
    numWdrawn <- selected.table[i]
    mesh.index <- which(names(selected.table[i])==names(universe.table))
    numW <- universe.table[mesh.index]
    numB <- length(p@universeGeneIds) - numW
    numDrawn <- length(p@geneIds)
    pvals[i] <- phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=FALSE)
  }
  
  ## Multiple testing correction  
  stats <- switch(p@pAdjust, BH = {
    p.adjust(pvals, "BH")
  }, QV = {
    suppressWarnings(fdrtool(pvals, statistic="pvalue", plot=FALSE, verbose=FALSE)$qval)
  }, lFDR = {
    suppressWarnings(fdrtool(pvals, statistic="pvalue", plot=FALSE, verbose=FALSE)$lfdr)
  }, none = pvals)

  ## Choose siginificantly enriched MeSH terms 
  mesh.list <- names(selected.table[which(stats < p@pvalueCutoff)])
  sort.index <- unlist(sort(stats[which(stats < p@pvalueCutoff)], index.return=TRUE)[2])
  pval.vec <- unlist(sort(stats[which(stats < p@pvalueCutoff)], index.return=TRUE)[1])
  mesh.list <- mesh.list[sort.index]

  ## Retrieve full name of MeSH category 
  tmp.mesh.cat <- unlist(strsplit(p@annotation, ""))
  mesh.cat <- tmp.mesh.cat[length(tmp.mesh.cat)]
  switch(mesh.cat, A = {
    mesh.full.cat <- "Anatomy"
  }, B = {
    mesh.full.cat <- "Organisms"
  }, C = {
    mesh.full.cat <- "Diseases"
  }, D = {
    mesh.full.cat <- "Drugs and chemicals"
  }, G = {
    mesh.full.cat <- "Phenomena and Processes"
  }, S = {
    mesh.full.cat <- "Substance Names"
  })

  ## Mapping  MeSH ID to MeSH term 
  mesh.df <- select(MeSHTERM, keys=mesh.list, cols=c("MESHID","MESHTERM","MESHCATEGORY"), keytype="MESHID")
  # remove categories not specified 
  mesh.df <- mesh.df[mesh.cat==mesh.df[,3],][,1:2]
  # remove MeSH terms appearing multiple times within same category 
  mesh.df <- mesh.df[!duplicated(mesh.df[,1]),]
  mesh.df <- mesh.df[sort.index,]
  mesh.df$PVALUE <- pval.vec
  
  new("MeSHHyperGResult",
      meshCategory=mesh.full.cat,
      annotation=p@annotation, 
      meshIds=mesh.df[,1],
      meshTerms=mesh.df[,2],
      pvalues=mesh.df[,3])
}


# QTL_summary function
#
# written by Brian Ritchey
#
# Returns gene name, their chromosome location, their position (in Mb), their protein description, distance from peak,... 
# and number of pubmed hits for term(s) of interest 
# First run QTL_genes...required for this function

QTL_summary <- 
  function(bayes, chr, pos1, pos2, search_term){ 
  require(rentrez)
  require(pubmed_count)
  
    o <- QTL_genes(bayes, chr, pos1, pos2) 
      if(!missing(bayes)){
        o <- o[,c(3,1,2,4,5)]   # Rearrange columns so gene name is first
}
else{
  o <- o[,c(3,1,2,4)]
}
  o <- o[order(o$Position),]  # Order by gene position
  o$Position <- round((o$Position/10^6), 
                        digits = 2)
  rownames(o) <- NULL
  eQTL_mathces <- match(o$Gene_Name, eQTL$Gene_Name)
  o$eQTL_LOD <- 0
  o$eQTL_LOD <- eQTL[eQTL_mathces, "LOD"]
  o$eQTL_LOD <- round(o$eQTL_LOD, digits = 1)  
  rownames(o) <- NULL
o
  gene <- as.character(o$Gene_Name)
    if(!missing(search_term)){
        pmed_count <- pubmed_count(gene, search_term)
        o <- cbind(o, pmed_count)
o
}

o 
} 



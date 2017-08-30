# flank_LOD function.

# written by Brian Ritchey.
# Used in conjunction with R/qtl package.
# Returns poisiton the based on a gene location LOD score of nearest flanking marker after a scanone analysis has been performed. 
# Utilizes R/qtl's "find.flanking" funciton.
#
# Aguments: x - should contain a data frame of gene chromosome and Mb position.
#           cross - the object of class cross used for R/qtl analysis.
#           results - the result of a scanone analysis stored as an object in your workspace.



flank_LOD <- 
  function(x, cross, result){
    require(qtl)
    output <- numeric(length = nrow(x))

    chr <- x$Chromosome[1]
    for(i in 1:nrow(x)){
  	    z <- find.flanking(cross, chr = chr, pos = x$Position[i])
      z[] <- lapply(z, as.character)
    	  z <- unname(z)                 
   	    z <- unlist(z)   # perform those operations to convert factor to char
    	  z <- z[3]  
    match <-  which(rownames(result) == z)
    output[i] <- round(result[match,3], digits = 1)
  }
  output
}





# flank_LOD example:
#
# flank_LOD_chr1 <- flank_LOD(load_CE_chr1_genes, 
#                             cross = peggy_chol_all_data, 
#                             result = load_CE.hk)



# pubmed_count function
#
# written by Brian Ritchey
#
# Arguments: x - list of genes (as characters)
#            search_term - term or terms of interest as character(s).  For multiple terms use c("term1", "term2")...
#
# Returns the number of PubMed hits for queries of "gene_name" and "search_term"


pubmed_count <- 
  function(x, search_term){
    require(rentrez)
    
    result <- matrix(nrow = length(x), ncol = length(search_term))
      for(i in 1:length(search_term)){
        for(j in 1:length(x)){
          result[j,i] <- entrez_search(db = "pubmed", 
                                      term = paste0(x[j], " AND ", search_term[i]), 
                                      retmax = 0)$count 
    }
  }
colnames(result) <- search_term
result 
} 






# gwas_catalog function
#
# written by Brian Ritchey


gwas_catalog <-

  function(gene){
  html_1 <- "https://www.ebi.ac.uk/gwas/api/search/downloads?q=text:%22"
  html_2 <- gene
  html_3 <- "%22&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&traitfilter[]=&dateaddedfilter=&facet=association"

  html <- paste0(html_1, html_2, html_3)
  results <- read.delim(html, sep = "\t", stringsAsFactors = F)

  results
}


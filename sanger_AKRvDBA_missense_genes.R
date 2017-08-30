# sanger_AKRvDBA_missense_genes function
#
# written by Brian Ritchey
#
# 

sanger_AKRvDBA_missense_genes <- 

  function(chr, pos1, pos2){
  require(rvest)
  require(rentrez)
  
  pos1 <- pos1 * 10^6
  pos2 <- pos2 * 10^6
  
  query_1 <- "https://www.sanger.ac.uk/sanger/Mouse_SnpViewer/rel-1211?gene=&context=0&loc="
  query_2 <- "%3A"
  query_3 <- "-"
  query_4 <- "&release=rel-1211&sn=missense_variant&st=akr_j&st=dba_2j"
  
  sanger_query <- paste0(query_1, chr, query_2, pos1, query_3, pos2, query_4)  
  sanger_html <- read_html(sanger_query) 
  gene_table <- sanger_html %>%  
    html_node("#t_snps_0 > div.scrollable > table")%>%
    html_table(header = T)
  
  # Exclude rows without dbSNP or where AKR allele == DBA allele
  exclude <- which(gene_table$dbSNP == "-" | gene_table[,6] == gene_table[,7])
  gene_table <- gene_table[-exclude,]
  rownames(gene_table) <- NULL
  gene_table
  genes <- unique(gene_table$Gene)
  rs_list <- vector("list", length(unique(gene_table$Gene)))
  
  for(i in 1:length(rs_list)){
    rs_list[[i]] <- gene_table$dbSNP[gene_table$Gene == (unique(gene_table$Gene)[i])]
  }
  rs_list
  names(rs_list) <- genes
rs_list

  # Calls and returns result from missense_for_provean function.
  # The code above filters the region and decides which gene rs to feed missense_for_provean function (rs_list)

  missense_for_provean(rs_list = rs_list, gene_table = gene_table)
}

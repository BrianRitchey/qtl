# missense_for_provean function
#
# written by Brian Ritchey
#


function(rs_list= rs_list, gene_table = gene_table){
  require(rentrez) 
  final <- vector("list")
  no_match <- vector("list")
  genes <- names(rs_list)
  
  for(rs_i in 1:length(rs_list)){
    missense_SNPs <- rs_list[[rs_i]]
    for(i in 1:length(missense_SNPs)){
     rs <- entrez_fetch(db = "snp", 
                        id = missense_SNPs[i], 
                        rettype = "flt")
      
      rs <- strsplit(rs, "residue")
      rs_split <- strsplit(rs[[1]][3], split = "\\|")
      
      pro_acc <- gsub(pattern = " prot_acc=", 
                                    replacement = "", 
                                    x = rs_split[[1]][4])
      
      pro_acc <- strsplit(pro_acc, split = "\n")
      pro_acc <- pro_acc[[1]][1]
      
      # break loop once an accession "NP_" is found
        if(length(grep(pattern = "NP_", x = pro_acc)) == 1){
             break
      }   
    }
    
 # Get protein sequence and accession name
 # If pro_seq doesn't exist, "NP_033256.2"... 
    # a soat1 acc which is the same in AKR + DBA/2 will be forced in
    pro_seq <- tryCatch({
      entrez_fetch(db = "protein", 
                               id = pro_acc, 
                               rettype = "fasta")}, 
      
      error = function(e){
        entrez_fetch(db = "protein", 
                                 id = "NP_033256.2", 
                                 rettype = "fasta")
      })
    
    pro_seq <- gsub(pro_seq, pattern = "\n", replacement = "")
    pro_seq_split <- strsplit(pro_seq, split = "]")
    pro_name <- pro_seq_split[[1]][1]
    pro_seq <- pro_seq_split[[1]][2]
    pro_seq_string <- strsplit(pro_seq, "")
    
    # Get amino acid position of residue change for each rs (AKR vs. DBA)
    # Also get referece and mutant alleles
    
    aa_positions <- character(length = length(missense_SNPs))
    ref_alleles <- character(length = length(missense_SNPs))
    mut_alleles <- character(length = length(missense_SNPs))
    
    
    # "rs30510802" is rs number that matches pro_acc above
    for(i in 1:length(missense_SNPs)){
      rs <- tryCatch({
        entrez_fetch(db = "snp", 
                     id = missense_SNPs[i], 
                     rettype = "flt")},
        error = function(e){  
          entrez_fetch(db = "snp", 
                                   id = "rs30510802", 
                                   rettype = "flt")
        })       
   
      rs <- strsplit(rs, "residue")
      rs_split <- strsplit(rs[[1]][3], split = "\\|")
      
      aa_positions[i] <-  trimws(
        gsub(pattern = "aa_position=", 
             replacement = "", 
             x = rs_split[[1]][2]))
      
      ref_alleles[i] <- trimws(
        gsub(pattern = "=", 
             replacement = "", 
             x = rs_split[[1]][1]))
      
      mutant_allele <- gsub(pattern = "=", 
                            replacement = "", 
                            x = rs[[1]][2])
      
      mutant_allele <- strsplit(mutant_allele, split = "\\|")
      
      mut_alleles[i] <- trimws(mutant_allele[[1]][1])
      
    }
    
    # Concatenate into Provean variant syntax (i.e. Y284H)
    # Account for no match conditions (NA) and orientation of ref/mutant in paste0
    
    ref_vec <- ref_alleles == pro_seq_string[[1]][as.integer(aa_positions)]
    mut_vec <- mut_alleles == pro_seq_string[[1]][as.integer(aa_positions)]
    
    pro_variants <- character(length = length(missense_SNPs))
    none <- character(length = length(missense_SNPs))
    
    for(i in 1:length(missense_SNPs)){
      
      if(is.na(ref_vec[i])){
        pro_variants[i] <- NA
        none[i] <- missense_SNPs[i] 
      }else if(ref_vec[i] == T){
        pro_variants[i] <- paste0(ref_alleles[i], aa_positions[i], mut_alleles[i])
      }else if (mut_vec[i] == T){
        pro_variants[i] <- paste0(mut_alleles[i], aa_positions[i], ref_alleles[i])
      }else{
        pro_variants[i] <- NA
        none[i] <- missense_SNPs[i] 
      }
    }
    pro_variants <- na.exclude(pro_variants)
    attributes(pro_variants) <- NULL
    no_match[[rs_i]] <- none
    pro_name <- paste0(pro_name, "]")
    output <- character(length = (3 + length(pro_variants)))
    output <- c(genes[rs_i], pro_name, pro_seq, pro_variants)
    final[[rs_i]] <- output
  }
  names(no_match) <- genes
  no_match

  names(final) <- genes
  final_df <- list_to_df(final)
  rownames(final_df) <- NULL
  final_rm <- which(is.na(final_df[,4]))

  if(length(final_rm) > 0 ){
    final <- final[-final_rm] 
  }
  no_match <- unlist(no_match, use.names = F)
  nm_rm <- which(gene_table$dbSNP == no_match, useNames = F)

  # In the case everything matches
  if(length(nm_rm) > 0){
    gene_table <- gene_table[-nm_rm,] 
  }
  rownames(gene_table) <- NULL
   final[[length(final) + 1]] <- gene_table
  # final[[length(final) + 1]] <- no_match
  # names(final) <- c(genes, "gene_table", "no_rs_match")
  
  final
}


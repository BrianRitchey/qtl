


QTL_genes <-
function(bayes, chr, pos1, pos2){
  if (!missing(bayes)){
    chr_subset <- subset(Mouse_Genome, 
                                          subset = Chromosome == as.integer(bayes$chr[1]))  
    final <- subset(chr_subset, 
                               subset = Position > (bayes$pos[1] * 10^6)  & Position < (bayes$pos[3] * 10^6))
    distance_to_peak <- abs(((bayes$pos[2] * 10^6) - final$Position) / 10^6)
    distance_to_peak <- round(distance_to_peak, digits = 2)
    output <- data.frame(final, distance_to_peak)
    rownames(final) <- NULL
    output
  }
  else{
    chr_subset <- subset(Mouse_Genome, 
                         subset = Chromosome == chr)
    final <- subset(chr_subset, 
                    subset = Position > (pos1 * 10^6) & Position < (pos2*10^6))
    output <- data.frame(final)
    output
  }  
output
}





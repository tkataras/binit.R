
#fuction itself
#datain file should be the "****_allValidPairs" file, this R script gives nice col.names, bins, creates data.frame with only useful info and then returns head
binit <- function(datain, sbin){
  datain <- read.table(datain)
  colnames(datain) <-  c("read_name", "ch_1_locus", "rd_pos_1", "strand1", "ch_2_locus", "rd_pos_2", "strand2", "frag_lngth","res_frag_1", "res_frag_2", "map_qual_1", "map_qual_2", "chrom_assign")
  
  bin_rd_pos_1_no_round <- datain$rd_pos_1/sbin
  bin_rd_pos_1 <- floor(bin_rd_pos_1_no_round) #using floor rather than ceiling to start at bin 0
    
  bin_rd_pos_2_no_round <- datain$rd_pos_2/sbin
  bin_rd_pos_2 <- floor(bin_rd_pos_2_no_round)
  binit_out <-data.frame(datain$read_name, datain$ch_1_locus, datain$rd_pos_1, bin_rd_pos_1, datain$ch_2_locus, datain$rd_pos_2, bin_rd_pos_2, datain$chrom_assign, row.names=TRUE)
  
  n <-head(binit_out)
  return(n)
  sink("binit_output", append=FALSE, split=FALSE )
  
  }



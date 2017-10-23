
#fuction itself
#datain file should be the "****_allValidPairs" file, this R script gives nice col.names, bins, creates data.frame with only useful info and then returns head
binit <- function(datain="~/bigdata/hic/data/SRR2240738/mm10/output/hic_results/data/mNPe/mNPe_allValidPairs", sbin=10000, sav.loc="./", ch.sizes="~/bigdata/hic/data/SRR2240738/mm10/mm10.chrom.sizes"){
 
   ###### ERROR CHECK, adding "/"
  if(!is.character(sav.loc)) stop("sav.loc must be a character string")
  split.string = unlist(strsplit(sav.loc,""))
  l = length(split.string)
  lastchar = split.string[l]
  if(lastchar!="/") datain = paste0(sav.loc,"/")
  
  ##### DATA MANIP
  # datain file should be the "****_allValidPairs"
  inputlist = strsplit(datain,"/")
  inputsplit = inputlist[[1]] # could use unlist command
  l = length(inputsplit)
  outfile_name = paste0(inputsplit[l],"_binit_out")
  ######
  
  
 

  
  #### Chromosome size
  #getting chrom sizes
  chrom.size<-read.table(ch.sizes)
  
  numbins <- ceiling(chrom.size[,2]/sbin)
  
  chromes = as.character(chrom.size$V1[-1])
  
  cumbins<-cumsum(numbins)
  cumbins = cumbins[-length(cumbins)]
  
  cumbins= as.data.frame(cumbins,ncol=1,row.names = chromes)
  
  
  
  
  
  datain <- read.table(datain)
  colnames(datain) <-  c("read_name", "rd1_ch", "rd_pos_1", "strand1", "rd2_ch", "rd_pos_2", "strand2", "frag_lngth","res_frag_1", "res_frag_2", "map_qual_1", "map_qual_2", "chrom_assign")
  
  bin_rd_pos_1<- floor(datain$rd_pos_1/sbin)
  bin_rd_pos_2<- floor(datain$rd_pos_2/sbin)
  
  #### spliting 0-0 to 0,0
  binit_in_progress<-data.frame(datain$read_name, datain$rd1_ch, datain$rd_pos_1, bin_rd_pos_1, datain$rd2_ch, datain$rd_pos_2, bin_rd_pos_2, datain$chrom_assign, row.names=TRUE)
  
  chrm_assign_col <-(binit_in_progress$datain.chrom_assign)
  char_chrm_assign_col <-as.character(chrm_assign_col) #needs to be char for strsplit
  splt_chrm_assign_col <- strsplit(char_chrm_assign_col,"-")
  
  unlist_chrm_assign_col <- unlist(splt_chrm_assign_col, recursive = TRUE, use.names = TRUE)
  
  df_chrom_assign_col <- matrix(unlist_chrm_assign_col,ncol=2, byrow=TRUE)
  
  type_both <- as.data.frame(df_chrom_assign_col)
  colnames(type_both) = c("type1","type2")
  #########
  
  binit_out<-data.frame(datain$read_name, datain$rd1_ch, datain$rd_pos_1, bin_rd_pos_1, datain$rd2_ch, datain$rd_pos_2, bin_rd_pos_2, datain$chrom_assign, type_both$type1, type_both$type2, row.names=TRUE)
  
  binit_preprocessed <- binit_out
  
  
  #####actual binning
  for(chr in chromes){
    add_index = cumbins[chr,] #adding the cumulative sum of bins to make the bin number outputs all sequential
    
    temp1 = binit_out$bin_rd_pos_1[ binit_out$datain.rd1_ch == chr ]
    binit_out$bin_rd_pos_1[ binit_out$datain.rd1_ch == chr ] = temp1 + add_index
    
    temp2 = binit_out$bin_rd_pos_2[ binit_out$datain.rd2_ch == chr ]
    binit_out$bin_rd_pos_2[ binit_out$datain.rd2_ch == chr ] = temp2 + add_index
  }
  
  #####making the reads sequential as well
  chrom.sizer = chrom.size[-length(chrom.size),]
  chrom.sizer <-as.data.frame(chrom.sizer,ncol=1,row.names = chromes) 
  #chrom.sizer <- colnames(chrom.sizer$V1)
  for(chr in chromes){
    add_indexr = chrom.sizer[chr,2] #adding the cumulative sum of reads to make the reads all sequential, not indexed on each chrom
    
    tempr1 = binit_out$datain.rd_pos_1[ binit_out$datain.rd1_ch == chr ]
    binit_out$datain.rd_pos_1[ binit_out$datain.rd1_ch == chr ] = tempr1 + add_indexr
    
    tempr2 = binit_out$datain.rd_pos_2[ binit_out$datain.rd2_ch == chr ]
    binit_out$datain.rd_pos_2[ binit_out$datain.rd2_ch == chr ] = tempr2 + add_indexr
  }
  
  
  # ##check to see if the sequntial binning worked
  # length(binit_out$bin_rd_pos_2)
  # n= 200000:300000
  # 
  # test <- binit_out$datain.rd_pos_1[n]>binit_out$datain.rd_pos_2[n+1]
  # class(test)
  # sum(test)
  # length(test)
  # length(n)
  
  
 # unique_bin_combos <-unique(binit_out[c("bin_rd_pos_1", "bin_rd_pos_2")])
##### making sure all bins are larger or = on the "left"
   
  verify <- binit_out$bin_rd_pos_1 > binit_out$bin_rd_pos_2
  #   head(verify) 
  # sum(verify)
  #   which(verify) 
  # max
#binit_out[458087,]  #here bin_rd_pos_1 is greater than bin_rd_pos_2
   
     
     holdingbin1 <-binit_out$bin_rd_pos_1
     holdingbin2 <-binit_out$bin_rd_pos_2
     
     holdingchr1 <-binit_out$datain.rd1_ch
     holdingchr2 <-binit_out$datain.rd2_ch
     
     holdingrdpos1 <- binit_out$datain.rd_pos_1
     holdingrdpos2 <- binit_out$datain.rd_pos_2
     
     binit_out$bin_rd_pos_1[verify] <- holdingbin2[verify] 
     binit_out$bin_rd_pos_2[verify] <- holdingbin1[verify]
     
     binit_out$datain.rd1_ch[verify] <- holdingchr2[verify]
     binit_out$datain.rd2_ch[verify] <- holdingchr1[verify]
     
     binit_out$datain.rd_pos_1[verify] <- holdingrdpos2[verify]
     binit_out$datain.rd_pos_2[verify] <- holdingrdpos1[verify]
     
    
     
     
  

    
   ##### selecting all rows in binit_out within unique pair of bins
    
     
      
     unique_bins <- unique(cbind(binit_out$bin_rd_pos_1, binit_out$bin_rd_pos_2)) #2 col matrix with all unique bin combos
    colnames(unique_bins) <- c("bin_1", "bin_2")
    
    empty_frame <- data.frame(matrix(0,ncol=9,nrow=nrow(unique_bins)))
    colnames(empty_frame) = c("a-a", "a-b", "b-a", "b-b", "a-x", "b-x", "x-a", "x-b", "x-x" )
    
    
    final_frame <- cbind.data.frame(unique_bins, empty_frame)#initializing
    
    
    
    tot_unique_bins = 1:(nrow(unique_bins))  
    for(n in tot_unique_bins){
        
        this_loop_combo_1 <- unique_bins[n,1] ###need to get all rows in binit_out where binit_out$bin_rd_pos_1 = this_loop_combo$bin_1 and binit_out$bin_rd_pos_2 = this_loop_combo$bin_2
        this_loop_combo_2 <- unique_bins[n,2]
        
        this_loop_rows <-binit_out[binit_out$bin_rd_pos_1 == this_loop_combo_1,]
        this_loop_rows_final <- this_loop_rows[this_loop_rows$bin_rd_pos_2 == this_loop_combo_2,]
    
       count_00 = sum(this_loop_rows_final$type_both.type1 == 0 &  this_loop_rows_final$type_both.type2 == 0 )
       count_01 = sum(this_loop_rows_final$type_both.type1 == 0 &  this_loop_rows_final$type_both.type2 == 1 )
       count_02 = sum(this_loop_rows_final$type_both.type1 == 0 &  this_loop_rows_final$type_both.type2 == 2 )
       count_10 = sum(this_loop_rows_final$type_both.type1 == 1 &  this_loop_rows_final$type_both.type2 == 0 )
       count_20 = sum(this_loop_rows_final$type_both.type1 == 2 &  this_loop_rows_final$type_both.type2 == 0 )
       count_12 = sum(this_loop_rows_final$type_both.type1 == 1 &  this_loop_rows_final$type_both.type2 == 2 ) 
       count_21 = sum(this_loop_rows_final$type_both.type1 == 2 &  this_loop_rows_final$type_both.type2 == 1 )
       count_11 = sum(this_loop_rows_final$type_both.type1 == 1 &  this_loop_rows_final$type_both.type2 == 1 )
       count_22 = sum(this_loop_rows_final$type_both.type1 == 2 &  this_loop_rows_final$type_both.type2 == 2 )
       
       #moving the "3s" or "conficting reads into unassigned
       count_00 = count_00 + sum(this_loop_rows_final$type_both.type1 == 3 &  this_loop_rows_final$type_both.type2 == 3)
       count_00 = count_00 + sum(this_loop_rows_final$type_both.type1 == 3 &  this_loop_rows_final$type_both.type2 == 0)
       count_00 = count_00 + sum(this_loop_rows_final$type_both.type1 == 0 &  this_loop_rows_final$type_both.type2 == 3)
       count_10 = count_10 + sum(this_loop_rows_final$type_both.type1 == 1 &  this_loop_rows_final$type_both.type2 == 3)
       count_20 = count_20 + sum(this_loop_rows_final$type_both.type1 == 2 &  this_loop_rows_final$type_both.type2 == 3)
       count_02 = count_02 + sum(this_loop_rows_final$type_both.type1 == 3 &  this_loop_rows_final$type_both.type2 == 2)
       count_01 = count_01 + sum(this_loop_rows_final$type_both.type1 == 3 &  this_loop_rows_final$type_both.type2 == 1)
       
     final_frame$`x-x`[n] = count_00
     final_frame$`x-a`[n] = count_01 
     final_frame$`x-b`[n] = count_02
     final_frame$`a-x`[n] = count_10
     final_frame$`b-x`[n] = count_20
     final_frame$`a-a`[n] = count_11
     final_frame$`b-b`[n] = count_22
     final_frame$`a-b`[n] = count_12
     final_frame$`b-a`[n] = count_21
       
        
        
        }
     
     
     
     write.table(final_frame, outfile_name , append=FALSE, sep="\t", quote=FALSE)

  
  
  
}




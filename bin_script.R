
#fuction itself
#datain file should be the "****_allValidPairs" file, this R script gives nice col.names, bins, creates data.frame with only useful info and then returns head
binit <- function(input_data="~/bigdata/hic/data/SRR2240738/mm10/output/hic_results/data/mNPe/mNPe_allValidPairs", sbin=10000, sav.loc="./", ch.sizes="~/bigdata/hic/data/SRR2240738/mm10/mm10.chrom.sizes"){
  
  ###### ERROR CHECK, adding "/"
  if(!is.character(sav.loc)) stop("sav.loc must be a character string")
  split.string = unlist(strsplit(sav.loc,""))
  l = length(split.string)
  lastchar = split.string[l]
  if(lastchar!="/") input_data = paste0(sav.loc,"/")
  
  ##### DATA MANIP, making sure input has "/" and creating outfile name based on the infile name
  # input_data file should be the "****_allValidPairs"
  inputlist = strsplit(input_data,"/") 
  inputsplit = inputlist[[1]] # could use unlist command
  l = length(inputsplit)
  outfile_name = paste0(inputsplit[l],"_binit_out_11_15")
  ######
  
  
  
  #### Chromosome size
  #getting chrom sizes and creating cumulative numbers of bins for each chromosome
  chrom.size<-read.table(ch.sizes)
  
  numbins <- ceiling(chrom.size[,2]/sbin)
  
  chromes = as.character(chrom.size$V1[-1])
  
  cumbins<-cumsum(numbins)
  cumbins = cumbins[-length(cumbins)]
  
  cumbins= as.data.frame(cumbins,ncol=1,row.names = chromes)
  
  
  
  ###introduce input data
  datain <- read.table(input_data)
  colnames(datain) <-  c("read_name", "rd1_ch", "rd_pos_1", "strand1", "rd2_ch", "rd_pos_2", "strand2", "frag_lngth","res_frag_1", "res_frag_2", "map_qual_1", "map_qual_2", "chrom_assign")
  
  
  
  ##### removing any incomplete cases or levels that appear in read 2 but not read 1 (this means the parent with most chormosmomes needs to be "read 2")
  datain[,5] <-factor(datain[,5], levels=c(levels(datain$rd1_ch)))#overwriting levels of rd2_ch with rd1_ch net change is rd2_ch loses chrY which is what i was trying to accomplish
  datain <- datain[complete.cases(datain),]#removing all rows with <NA> cases (in test dataset only 83 of 859703)
  
  
  ##### removing all inter-chromosomal reads
  datain_no_inter <- datain[datain$rd1_ch == datain$rd2_ch,]
  
  
  ### preliminary bining, the bin numbers this produces will be indexed on each chromosome
  bin_rd_pos_1<- floor(datain_no_inter$rd_pos_1/sbin) #gives chromosome indexed bin number
  bin_rd_pos_2<- floor(datain_no_inter$rd_pos_2/sbin)
  
  # #### spliting 0-0 to 0,0 as two columns
  # binit_in_progress<-data.frame(datain_no_inter$read_name, datain_no_inter$rd1_ch, datain_no_inter$rd_pos_1, bin_rd_pos_1, datain_no_inter$rd2_ch, datain_no_inter$rd_pos_2, bin_rd_pos_2, datain_no_inter$chrom_assign, row.names=TRUE)
  # 
  # chrm_assign_col <-(binit_in_progress$datain_no_inter.chrom_assign)
  # char_chrm_assign_col <-as.character(chrm_assign_col) #needs to be char for strsplit
  # splt_chrm_assign_col <- strsplit(char_chrm_assign_col,"-")
  # 
  # unlist_chrm_assign_col <- unlist(splt_chrm_assign_col, recursive = TRUE, use.names = TRUE)
  # 
  # df_chrom_assign_col <- matrix(unlist_chrm_assign_col,ncol=2, byrow=TRUE)
  # 
  # type_both <- as.data.frame(df_chrom_assign_col)
  # colnames(type_both) = c("type1","type2")
  #########
  
  ##### data processsing cont.
  binit_out<-data.frame(datain_no_inter$read_name, datain_no_inter$rd1_ch, datain_no_inter$rd_pos_1, bin_rd_pos_1, datain_no_inter$rd2_ch, datain_no_inter$rd_pos_2, bin_rd_pos_2, datain_no_inter$chrom_assign, row.names=NULL)
  
  #### creating placeholder
  binit_preprocessed <- binit_out
  
  
  ##### making all bins seqential along entire genome
  for(chr in chromes){
    add_index = cumbins[chr,] #adding the cumulative sum of bins to make the bin number outputs all sequential
    #chr <- as.factor(chr)
    
    temp1 = binit_out$bin_rd_pos_1[ binit_out$datain_no_inter.rd1_ch == chr ]
    binit_out$bin_rd_pos_1[ binit_out$datain_no_inter.rd1_ch == chr ] = temp1 + add_index
    
    temp2 = binit_out$bin_rd_pos_2[ binit_out$datain_no_inter.rd2_ch == chr ]
    binit_out$bin_rd_pos_2[ binit_out$datain_no_inter.rd2_ch == chr ] = temp2 + add_index
  }
  
  #####making the reads sequential as well
  chrom.sizer = chrom.size[-length(chrom.size),]
  chrom.sizer <-as.data.frame(chrom.sizer,ncol=1,row.names = chromes) 
  #chrom.sizer <- colnames(chrom.sizer$V1)
  for(chr in chromes){
    add_indexr = chrom.sizer[chr,2] #adding the cumulative sum of reads to make the reads all sequential, not indexed on each chrom
    
    tempr1 = binit_out$datain.rd_pos_1[ binit_out$datain_no_inter.rd1_ch == chr ]
    binit_out$datain.rd_pos_1[ binit_out$datain_no_inter.rd1_ch == chr ] = tempr1 + add_indexr
    
    tempr2 = binit_out$datain.rd_pos_2[ binit_out$datain_no_inter.rd2_ch == chr ]
    binit_out$datain.rd_pos_2[ binit_out$datain_no_inter.rd2_ch == chr ] = tempr2 + add_indexr
  }
  
 
  
  #####  SORTING then selecting all rows in binit_out within unique pair of bins
  sorted_binit_out <- binit_out[order(binit_out$bin_rd_pos_1,binit_out$bin_rd_pos_2),]
  
  unique_bins <- unique(cbind(sorted_binit_out$bin_rd_pos_1, sorted_binit_out$bin_rd_pos_2)) #2 col matrix with all unique bin combos
  colnames(unique_bins) <- c("bin_1", "bin_2")
  
  ####initializing the output data.frame
  empty_frame <- data.frame(matrix(0,ncol=9,nrow=nrow(unique_bins)))
  colnames(empty_frame) = c("a-a", "a-b", "b-a", "b-b", "a-x", "b-x", "x-a", "x-b", "x-x" )
  final_frame <- cbind.data.frame(unique_bins, empty_frame)#initializing
  
  
  
  ####New counting scheme (turning the factor column "chrom.assign" into counts in several columns)
  unique_bin_row <-1
  
  length_binit_out <- 1:(length(sorted_binit_out$datain_no_inter.chrom_assign))
  for(row in length_binit_out){
    
    ##if the current row bin combo doesnt match the unique bin combo in final_frame, this moves to the next unique bin line and adds counts as if it had started this loop with the new value of unique_bin_row
    if(sorted_binit_out$bin_rd_pos_1[row] != final_frame$bin_1[unique_bin_row]||(sorted_binit_out$bin_rd_pos_2[row] != final_frame$bin_2[unique_bin_row])){ 
      unique_bin_row <- sum(unique_bin_row + 1)
      if(sorted_binit_out$bin_rd_pos_1[row] == final_frame[unique_bin_row,1]&&(sorted_binit_out$bin_rd_pos_2[row] == final_frame[unique_bin_row,2])){
        
        if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-0"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1) 
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-1"){final_frame$`x-a`[unique_bin_row]= sum(final_frame$`x-a`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-2"){final_frame$`x-b`[unique_bin_row]= sum(final_frame$`x-b`[unique_bin_row]+1) 
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-0"){final_frame$`a-x`[unique_bin_row]= sum(final_frame$`a-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-0"){final_frame$`b-x`[unique_bin_row]= sum(final_frame$`b-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-2"){final_frame$`a-b`[unique_bin_row]= sum(final_frame$`a-b`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-1"){final_frame$`b-a`[unique_bin_row]= sum(final_frame$`b-a`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-1"){final_frame$`a-a`[unique_bin_row]= sum(final_frame$`a-a`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-2"){final_frame$`b-b`[unique_bin_row]= sum(final_frame$`b-b`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-3"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-0"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-3"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-3"){final_frame$`a-x`[unique_bin_row]= sum(final_frame$`a-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-3"){final_frame$`b-x`[unique_bin_row]= sum(final_frame$`b-x`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-2"){final_frame$`x-b`[unique_bin_row]= sum(final_frame$`x-b`[unique_bin_row]+1)
        }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-1"){final_frame$`x-a`[unique_bin_row]= sum(final_frame$`x-a`[unique_bin_row]+1)
        }else{stop("unexpected data in chrom_assign. see line 212")}
        
        
      }
      #this if statement determines if the current row(read) falls within the first pair of unique bins
    }else if(sorted_binit_out$bin_rd_pos_1[row] == final_frame[unique_bin_row,1]&&(sorted_binit_out$bin_rd_pos_2[row] == final_frame[unique_bin_row,2])){
      
      #this looks terrible, but it just goes through each possible chom assignment and says for the current row, add one to right col in final_frame, if that cell in final frame already has some, add one more
      if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-0"){final_frame$`x-x`[unique_bin_row]<- sum(final_frame$`x-x`[unique_bin_row]+1) 
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-1"){final_frame$`x-a`[unique_bin_row]= sum(final_frame$`x-a`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-2"){final_frame$`x-b`[unique_bin_row]= sum(final_frame$`x-b`[unique_bin_row]+1) 
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-0"){final_frame$`a-x`[unique_bin_row]= sum(final_frame$`a-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-0"){final_frame$`b-x`[unique_bin_row]= sum(final_frame$`b-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-2"){final_frame$`a-b`[unique_bin_row]= sum(final_frame$`a-b`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-1"){final_frame$`b-a`[unique_bin_row]= sum(final_frame$`b-a`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-1"){final_frame$`a-a`[unique_bin_row]= sum(final_frame$`a-a`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-2"){final_frame$`b-b`[unique_bin_row]= sum(final_frame$`b-b`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-3"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-0"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "0-3"){final_frame$`x-x`[unique_bin_row]= sum(final_frame$`x-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "1-3"){final_frame$`a-x`[unique_bin_row]= sum(final_frame$`a-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "2-3"){final_frame$`b-x`[unique_bin_row]= sum(final_frame$`b-x`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-2"){final_frame$`x-b`[unique_bin_row]= sum(final_frame$`x-b`[unique_bin_row]+1)
      }else if(sorted_binit_out$datain_no_inter.chrom_assign[row] == "3-1"){final_frame$`x-a`[unique_bin_row]= sum(final_frame$`x-a`[unique_bin_row]+1)
      }else{stop("unexpected data in chrom_assign. see line 255")}
      
      
      
    }
    
  }
  
 
  
  ###adding back in chromosomes from sorted_binit_out, only using the unique bin rows
  correct_chr_1 <- data.frame(matrix(0,ncol=1,nrow=nrow(unique_bins)))
  correct_chr_2 <- data.frame(matrix(0,ncol=1,nrow=nrow(unique_bins)))
  
  unique_bin_row2 <-1
  
  length_binit_out2 <- 1:(length(sorted_binit_out$datain_no_inter.chrom_assign))
  
  for(row2 in length_binit_out2){
    
    if(sorted_binit_out$bin_rd_pos_1[row2] != final_frame$bin_1[unique_bin_row2]||(sorted_binit_out$bin_rd_pos_2[row2] != final_frame$bin_2[unique_bin_row2])){
      unique_bin_row2 = unique_bin_row2 + 1
      correct_chr_1[unique_bin_row2,1] <- as.character(sorted_binit_out$datain_no_inter.rd1_ch[row2])
      correct_chr_2[unique_bin_row2,1] <- as.character(sorted_binit_out$datain_no_inter.rd2_ch[row2])
    }
    #this if statement determines if the current row(read) falls within the first pair of unique bins
    else if(sorted_binit_out$bin_rd_pos_1[row2] == final_frame[unique_bin_row2,1]&&(sorted_binit_out$bin_rd_pos_2[row2] == final_frame[unique_bin_row2,2])){
      correct_chr_1[unique_bin_row2,1] <- as.character(sorted_binit_out$datain_no_inter.rd1_ch[row2])
      correct_chr_2[unique_bin_row2,1] <- as.character(sorted_binit_out$datain_no_inter.rd2_ch[row2])
      
    }}
 
  #### output
  final_frame_out <- cbind.data.frame(final_frame$bin_1, final_frame$bin_2,correct_chr_1,correct_chr_2, final_frame$`a-a`,final_frame$`b-b`, final_frame$`a-b`,final_frame$`b-a`,final_frame$`a-x`,final_frame$`b-x`,final_frame$`x-a`,final_frame$`x-b`,final_frame$`x-x`)
  colnames(final_frame_out) = c("bin_1", "bin_2","chr_1", "chr_2", "a-a","b-b", "a-b", "b-a", "a-x", "b-x", "x-a", "x-b", "x-x" )
  
  write.table(final_frame_out, outfile_name , append=FALSE, sep="\t", quote=FALSE)
  
  
}




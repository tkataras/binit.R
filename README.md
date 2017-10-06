# binit.R
adds column for bin of each read after giving col names to whole file
output has 8 columns

to use, example:
output <- binit(~/bigdata/hic/etc..., 10000)

currently reciving this error code:

> write(binit_out,"binit_out", ncolumns=8, append=FALSE, sep="\t")
Error in cat(x, file = file, sep = c(rep.int(sep, ncolumns - 1), "\n"),  : 
  argument 1 (type 'list') cannot be handled by 'cat'

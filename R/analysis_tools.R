
#function for CPM normalising read counts, requires data.frame from classify_windows() and a vector with the total sample read count in the same order as the df as input
#returns the input dataframe with a new noramlised read count column 

normalise_CPM <- function(df, nreadcount_vec){
  
  out_df <- data.frame()
  for (i in seq(length(unique(df$samplename)))){
    sample_df <- df[df$samplename == unique(df$samplename)[i],]
    nreads <- nreadcount_vec[i]
    sample_df$CPM <- as.numeric(sample_df$nReads)/(as.numeric(nreads)/1000000)
    sample_df$log2CPM <- log2(sample_df$CPM+1)
    out_df <- rbind(out_df,sample_df)
  }
  return(out_df)
}
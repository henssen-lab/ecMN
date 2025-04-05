##
#counts the number of positive windows per sample and returns a dataframe with samplename and number of positive windows
count_positive_windows <- function(df){
  sample_vec <- c()
  n_positive_windows_vec <- c()
  
  for (sample in unique(df$samplename)){
    sample_df <- df[df$samplename==sample,]
    n_positive_windows <- dim(sample_df[sample_df$window_call_final==1,])[1]
    sample_vec <- c(sample_vec,sample)
    n_positive_windows_vec <- c(n_positive_windows_vec,n_positive_windows)
  }
  
  out_df <- data.frame(sample_vec,n_positive_windows_vec)
  colnames(out_df) <- c("sample","n_positive_windows_vec")
  return(out_df)
}

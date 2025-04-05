#function to compute mean, variance and to determine thresholds for window calling
#fits a nagative binomial for reads per bin and determines coverage % threshold by the nth perentile
#takes as input an dataframe in the following format: chr start end n_reads_in_bin n_bases_covered bin_size perc_bases_covered sample_name
#additionally, takes as input propability P, and the nth percentile
#returns a vector in following format: sample_name mean_reads variance_reads threshold_reads mean_perc_covered variance_perc_covered threshold_perc_covered
compute_thresholds <- function(df,P,percentile, plot, out_dir_plots) {
  
  colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename")
  
  #compute mean and variance
  mean_reads <- mean(df$nReads)
  variance_reads <- var(df$nReads)
  mean_perc_covered <- mean(df$percent_covered)
  variance_perc_covered <- var(df$percent_covered)
  
  ####
  #read threshold
  ###
  #fit negative binomial 
  fit_params <- fitdist(df$nReads, "nbinom")
  fit_PMF <- dnbinom(0:30000, size=fit_params$estimate[1], mu=fit_params$estimate[2]) #probability mass function
  #compute cumulative mass function for the first 
  fit_CMF <- pnbinom(order(df$nReads)[1:5000], size=fit_params$estimate[1], mu=fit_params$estimate[2]) 
  
  index_read_theshold <- Position(function(x) x > 1-P,fit_CMF) #determine the first index for which 1-P of the CMF is greater than the predefined propability P (the probability of n reads in a bin is smaller than P)
  read_thresh <- order(df$nReads)[1:5000][index_read_theshold]
  
  ####
  #coverage threshold
  ###
  #as for now, the the threshold for the coverage is defined as the 0.999 percentile
  coverage_thresh <- quantile(df$percent_covered, probs = percentile)
  #returns a vector in following format: sample_name mean_reads variance_reads threshold_reads mean_perc_covered variance_perc_covered threshold_perc_covered
  
  out_vec <- c(df$samplename[1],mean_reads,variance_reads,read_thresh,mean_perc_covered,variance_perc_covered,coverage_thresh)
  
  ##plotting
  #plots histogram with fitted probability mass function and cumulative mass function
  if (plot == T) {
    pdf(paste(out_dir_plots,"/",df$samplename[1],"_PMF_CMF.pdf", sep = ""), width = 10, height=5)
    plot1 = ggplot(df, aes(x = nReads)) + 
      geom_histogram(aes(y=..density..),colour = "black", fill = "white", binwidth = 1) +
      geom_line(data=as.data.frame(fit_PMF), aes(1:length(fit_PMF),fit_PMF), color="blue", linewidth=0.5) + xlim(0,100)
    plot2 = ggplot(as.data.frame(fit_CMF[1:2500])) + 
      geom_point(aes(seq(1,length(fit_CMF[1:2500])), fit_CMF[1:2500])) + xlab("n reads") + ylab("cumulative probability")
    ggarrange(plot1, plot2, ncol = 2, nrow = 1, align = 'h')
    print(ggarrange(plot1, plot2, ncol = 2, nrow = 1, align = 'h'))
    dev.off()
  }
  
  
  return (out_vec)
}


#takes the list of dataframes from load_bed_files(), returns a table with thresholds
#takes list of dataframes, P probability for read count, percentile for coverage, plot, plot_dir; as input
generate_threshold_df <- function(bed_file_list,P=0.001,percentile=0.999,plot=F,out_dir_plots=""){
  
  sample_thresholds_list <- lapply(bed_file_list, function(x) {
    thresh_vec <- compute_thresholds(x,P=P,percentile=percentile,plot=plot,out_dir_plots= out_dir_plots)
  } )
  sample_thresholds <- data.frame(do.call(rbind, sample_thresholds_list))
  colnames(sample_thresholds) <- c("Sample","mean_reads","variance_reads","threshold_reads","mean_perc_covered","variance_perc_covered","threshold_perc_covered")

  return(sample_thresholds)
  
}

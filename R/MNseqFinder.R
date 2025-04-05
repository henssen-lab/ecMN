##script containing functions for sMNseqfinder

#requires tidyr ggplot2 fitdistrplus MASS ggpubr



###############################################################################
#function to compute mean, variance and to determine thresholds for window calling
#fits a nagative binomial for reads per bin and determines coverage % threshold by the nth perentile
#takes as input an dataframe in the following format:
# chr start end n_reads_in_bin n_bases_covered bin_size perc_bases_covered sample_name
#additionally, takes as input propability P, and the nth percentile
#returns a vector in following format: sample_name mean_reads variance_reads threshold_reads mean_perc_covered variance_perc_covered threshold_perc_covered
compute_thresholds <- function(df,P=0.001,percentile=0.999, plot=F, our_dir_plots="") {
  
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
  if (plot == T) {
    pdf(paste(our_dir_plots,"/",df$samplename[1],"_PMF_CMF.pdf", sep = ""), width = 10, height=5)
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


##this function counts adjacent windows and classifices windows as Edge or body segment, it ultimately call if a window is positive
count_adjacent_windows <- function(df_thresholded_windows, window_proximity_n=4){
  df_thresholded_windows <- df_thresholded_windows[df_thresholded_windows$window_class ==1,]
  colnames(df_thresholded_windows) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename","window_class")
  colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename")
  df_thresholded_windows_index <- as.numeric(rownames(df_thresholded_windows))
  df_thresholded_windows_index_neighbour_distance <- data.frame()
  #find the nearest up/downstream positive window
  for (i in seq(1:length(df_thresholded_windows_index))){
    if (i == 1) {
      neighbour_upstream <- 0
      neighbour_downstream <- (df_thresholded_windows_index[i]-df_thresholded_windows_index[i+1])*-1
    } else if (i == length(df_thresholded_windows_index)) {
      neighbour_upstream <- df_thresholded_windows_index[i]-df_thresholded_windows_index[i-1]
      neighbour_downstream <- 0
    } else {
      neighbour_upstream <- df_thresholded_windows_index[i]-df_thresholded_windows_index[i-1]
      neighbour_downstream <- (df_thresholded_windows_index[i]-df_thresholded_windows_index[i+1])*-1
    }
    neighbour_vec <- c(neighbour_upstream,neighbour_downstream)
    df_thresholded_windows_index_neighbour_distance <- rbind(df_thresholded_windows_index_neighbour_distance,neighbour_vec)
  }
  df_thresholded_windows_index_neighbour_distance$index <- df_thresholded_windows_index
  colnames(df_thresholded_windows_index_neighbour_distance) <- c("upstream_nearest_window","downstream_nearest_window","index")
  #define edges based on the nearest window. I.e. if the nearest upstreamwindow is far away, but downstream window is in proximity it is classified as left_edge
  #if there is a window in proximity up and dowstream it is classified as body (middle)
  #if there is a near window upstream, but not downstream, it is classified as right_edge
  #0:NA 1:edge 2:body
  
  edge_body_class_vec <- c()
  for (i in seq(1:length(df_thresholded_windows_index_neighbour_distance$index))) {
    left_window <- df_thresholded_windows_index_neighbour_distance[i,]$upstream_nearest_window
    right_window <- df_thresholded_windows_index_neighbour_distance[i,]$downstream_nearest_window
    if (left_window >= window_proximity_n & right_window <= window_proximity_n){
      edge_body_class_vec <- c(edge_body_class_vec,1)
    }
    else if (left_window <= window_proximity_n & right_window >= window_proximity_n){
      edge_body_class_vec <- c(edge_body_class_vec,1)
    }
    else if (left_window <= window_proximity_n & right_window <= window_proximity_n){
      edge_body_class_vec <- c(edge_body_class_vec,2)
    }
    else {
      edge_body_class_vec <- c(edge_body_class_vec,0)
    }
  }
  df_thresholded_windows_index_neighbour_distance$edge_body_class <- edge_body_class_vec
  #compute edge_body_class sum by summing up left and right facing window scores; 
  #sum of >2 denotes postive called region
  
  left_sum_vec <- c()
  right_sum_vec <- c()
  final_window_call <- c()
  for (i in seq(1:length(df_thresholded_windows_index_neighbour_distance$index))) {
    if (i==1){
      left_sum <- 0
      right_sum <- df_thresholded_windows_index_neighbour_distance[i,]$edge_body_class + df_thresholded_windows_index_neighbour_distance[i+1,]$edge_body_class
    }
    else if (i==length(df_thresholded_windows_index_neighbour_distance$index)){
      left_sum <- df_thresholded_windows_index_neighbour_distance[i,]$edge_body_class + df_thresholded_windows_index_neighbour_distance[i-1,]$edge_body_class
      right_sum <- 0
    }
    else {
      left_sum <- df_thresholded_windows_index_neighbour_distance[i,]$edge_body_class + df_thresholded_windows_index_neighbour_distance[i-1,]$edge_body_class
      right_sum <- df_thresholded_windows_index_neighbour_distance[i,]$edge_body_class + df_thresholded_windows_index_neighbour_distance[i+1,]$edge_body_class
    }
    left_sum_vec <- c(left_sum_vec, left_sum)
    right_sum_vec <- c(right_sum_vec,right_sum)
    ##if left or right_sum is > 2, then the window is called as postive
    ##I know that this is sloppy code, but its late and I am tired, at least it works
    if (left_sum>2|right_sum>2) {
      final_window_call <- c(final_window_call,1)
    }
    else {
      final_window_call <- c(final_window_call,0)
    }
  }
  df_thresholded_windows_index_neighbour_distance$window_call_final <- final_window_call
  df_thresholded_windows_index_neighbour_distance$left_sum <- left_sum_vec
  df_thresholded_windows_index_neighbour_distance$right_sum <- right_sum_vec
  return(df_thresholded_windows_index_neighbour_distance)
}
#function to handle samples without called positive windows, returns negative dummy data
failed_window_calling_df <- function(){
  df <- data.frame(list(0,0,1,0,0,0,0))
  colnames(df) <- c("upstream_nearest_window","downstream_nearest_window","index","edge_body_class","window_call_final","left_sum","right_sum")
  return(df)
}


##this function samples all windows which meet the threshold conditions and classifies windows into edge/body/non-call
#calls count_adjacent_windows()
classify_window <- function(df,threshold_df, plot=F, plot_outdir=""){
  colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename")
  sample <- df$samplename[1]
  threshold <- threshold_df[threshold_df$Sample == sample,]
  read_threshold <- as.numeric(threshold$threshold_reads)
  coverage_treshold <- as.numeric(threshold$threshold_perc_covered)
  df$window_class <- 0 #init column for window classificaiton
  #set window_class to 1 for all windows satisfying filter conditions (>= thresholds), trycatch in case there is an error (0 called windows)
  tryCatch(df[(df$nReads>=read_threshold&df$percent_covered>=coverage_treshold),]$window_class <- 1, error=function(e) df=df)
  window_calls_df <- tryCatch(count_adjacent_windows(df[df$window_class == 1,],window_proximity_n=4), error=function(e) failed_window_calling_df())
  
  ##matching
  df$window_call_final <- 0
  
  tryCatch(df[window_calls_df[window_calls_df$window_call_final>0,]$index,]$window_call_final <- 1, error=function(e) df=df)
  if (plot==T){
    plot_scatter(df,threshold_df,plot_outdir)
  }
  return(df)
}



#function to plot scatters for each sample (nReads vs. percent covered) and to also mark positive windows

plot_scatter <- function(df, threshold_df, out_dir_scatter){
  sample <- df$samplename[1]
  threshold <- threshold_df[threshold_df$Sample == sample,]
  y_intercept <- as.numeric(threshold$threshold_reads)
  x_intercept <- as.numeric(threshold$threshold_perc_covered)
  
  df$Window <- ""
  tryCatch(df[df$window_call_final==1,]$Window <- "Pass", error=function(e) df=df)
  df[df$window_call_final==0,]$Window <- "Fail"
  
  pdf(paste(out_dir_scatter,"/",sample,"_scatter.pdf", sep = ""), width = 10, height=10)
  plot_scatter <- ggplot(data=df, aes(x=percent_covered, nReads, color = Window)) +
    geom_point() + 
    scale_color_manual(values = c("Fail" = "black", "Pass" = "steelblue")) +
    geom_rug(col="steelblue",alpha=0.1, size=1.5)+
    geom_hline(yintercept=y_intercept) + geom_vline(xintercept = x_intercept)+
    ylab("n(Reads)") + xlab("% window covered") + ggtitle(sample)  + theme_bw()
  print(plot_scatter)
  dev.off()

##functions for normalisation and qc
normalise_CPM <- function(df, read_count_df){
  colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename","window_class","window_call_final")
  sample <- df$samplename[1]
  nread <- read_count_df[read_count_df$well == sample,]$read_count #total librarysize
  df$CPM <- df$nReads/(nread/1000000)
  df$log2CPM <- log2(df$CPM+1)
  
  return(df)
  
}

count_positive_windows <- function(df){
  sample <- df$samplename[1]
  n_positive_windows <- dim(df[df$window_call_final==1,])[1]
  return(c(sample,n_positive_windows))
  
}


#function for loading bed files
#must be in the following format: chr start end n(Reads) n(bases covered) bin_length %(bin_covered)
#meta file must be in the following format: path_to_file sample_name class 
load_bed_files_lead_meta <- function(path_bed_dir, path_meta){
  files_bed <- list.files(path=path, pattern="*.bed", full.names=TRUE, recursive=FALSE)
  beds <- data.frame()
  bed_file_list <- lapply(files_bed, function(x) {
    sample <- read.table(x, header=F) 
    sample$V8 <- str_split_1(x, "_")[7]
    beds <- rbind(beds, sample)
  })
  
}



files_bed
sample_names <- c()
for (path in files_bed){
  print(str_split_1(path,"_")[7])
  sample_names <- c(sample_names,str_split_1(path,"_")[7])
}
c(files_bed)
data.frame(files_bed,sample_names)
meta_df <- data.frame(files_bed,sample_names)
colnames(meta_df) <- c("path","well")
meta_df_df <- merge(meta_df, meta, by='well')

meta_df_df[,c("path","well","class")]
write.table(meta_df_df[,c("path","well","class")], "/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/meta/meta_pilot.txt", sep= "\t", row.names=FALSE, col.names=FALSE, quote=FALSE) 

files_bed <- list.files(path="/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/25kb_bins_q20", pattern="*.bed", full.names=TRUE, recursive=FALSE)
name_sample <- str_split_1("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/25kb_bins_q20/sorted_P3469_scMN_F5_S35_dedup_25kb_bin_q20.bed", "_")
meta_df


name_sample[length(name_sample)]
beds <- data.frame()
bed_file_list <- lapply(files_bed, function(x) {
  sample <- read.table(x, header=F) 
  sample$V8 <- str_split_1(x, "_")[7]
  beds <- rbind(beds, sample)
})


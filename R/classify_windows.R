#functions to classify windows


##this function counts adjacent windows and classifies windows as Edge or body segment, it ultimately calls if a window is positive
#window_proximity_n defines the search radius for a adjacent positive window defined by the user-set thresholds: i.e. for n=4: at each thresholded window, we search for an adjacent window within n bin sizes.
#classification as body or edge segment is dependent on the search radius
count_adjacent_windows <- function(df_thresholded_windows, window_proximity_n){
  df_thresholded_windows <- df_thresholded_windows[df_thresholded_windows$window_class ==1,]
  colnames(df_thresholded_windows) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename","window_class")
  #colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename")
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
  #sum of >2 denotes postive called region (final)
  
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
call_windows <- function(df,threshold_df, plot, plot_outdir,window_proximity_n){
  colnames(df) <- c("chr","start","end","nReads","nBases_covered","binsize","percent_covered","samplename")
  sample <- df$samplename[1]
  threshold <- threshold_df[threshold_df$Sample == sample,]
  read_threshold <- as.numeric(threshold$threshold_reads)
  coverage_treshold <- as.numeric(threshold$threshold_perc_covered)
  df$window_class <- 0 #init column for window classificaiton
  #set window_class to 1 for all windows satisfying filter conditions (>= thresholds), trycatch in case there is an error (in the case of 0 called windows)
  tryCatch(df[(df$nReads>=read_threshold&df$percent_covered>=coverage_treshold),]$window_class <- 1, error=function(e) df=df)
  window_calls_df <- tryCatch(count_adjacent_windows(df[df$window_class == 1,],window_proximity_n), error=function(e) failed_window_calling_df())
  
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
}

#calls all the above functions
#takes the bed_file_list from load_bed_files() and the thresholdtable from generate_threshold_df() as input, as well as the number radius search size
classify_windows <- function(bed_file_list,sample_thresholds,plot = F, plot_outdir = "",window_proximity_n=4){
  
  window_calls <- lapply(bed_file_list, function(x) {
    call_windows(x,sample_thresholds,plot, plot_outdir,window_proximity_n) #samples, which dont have any positive windows will throw an error, tryCatch to continue the function, returns df with window_class = 0 window_call_final = 0
  } )
  called_window_bed <- data.frame(do.call(rbind, window_calls))
  return(called_window_bed)
}

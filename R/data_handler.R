##data_handling


#function read in bed files from a meta table with structure: path  sample_name class
#returns a list of each bed file as data.frame
load_bed_files <- function(path_meta){
  
  #read files
  meta_df <- read.csv(path_meta, header = F, sep = "\t", col.names = c("path","sample_name","class"))
  
  
  beds <- data.frame()
  bed_file_list <- lapply(meta_df$path, function(x) {
    print(paste("Reading: ",x))
    sample <- read.table(x, header=F) 
    sample$V8 <- meta_df[meta_df["path"]==x,]$sample_name
    colnames(sample) <- c("chr","start","end","n_reads","n_bases","bin_size","perc_covered","sample")
    beds <- rbind(beds, sample)
  })
  return(bed_file_list)
}
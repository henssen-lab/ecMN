library(tidyr)
library(stringr)
library(ggplot2)
library(fitdistrplus)
library(MASS)
library(ggpubr)
library(ggsignif)
library(DescTools)


source("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/R/data_handler.R")
source("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/R/compute_thresholds.R")
source("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/R/classify_windows.R")
source("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/R/QC_tools.R")
source("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/R/analysis_tools.R")



#0. load meta table and reconstructions
meta <- read.csv("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/meta/meta_pilot.txt", sep="\t", header = F, col.names = c("path","sample","class"))
reconstrucion_bed <- read.csv("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/meta/TR14_hg38_reconstruction_decoil_with_CDK4_shasta_assembly.bed.txt", header = 1, sep = "\t")
#1. read files
bed_file_list <- load_bed_files("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/meta/meta_pilot.txt")
#2. compute thresholds
sample_threshold_df <- generate_threshold_df(bed_file_list,P=0.01,percentile=0.99,plot=F,out_dir_plots="/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/Tests/CMF")
#3. classify windows
window_calls <- classify_windows(bed_file_list,sample_threshold_df,plot = F, plot_outdir = "/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/MNseqFinder/Tests/plots",window_proximity_n=4)
#4. add read counts from external files and add window counts to the meta table
window_calls_counts <- count_positive_windows(window_calls)

files <- list.files(path="/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/uniq_read_counts", pattern="*.count", full.names=TRUE, recursive=FALSE)
read_counts <- data.frame()
for (file in files){
  sample <- read.table(file, header=F)
  sample$V2 <- str_split_1(file, "_")[7]
  print(sample)
  read_counts <- rbind(read_counts,sample)
}
colnames(read_counts) <- c("n_unique_reads","sample")

meta_qc <- merge(meta, window_calls_counts, by='sample')
meta_qc <- merge(meta_qc,read_counts , by='sample')

###
## analysis
###

library(GenomicRanges)
library(regioneR)
library(UpSetR)
library(karyoploteR)

#5. CPM normalise read counts
window_calls_counts_normalised <- normalise_CPM(window_calls,meta_qc$n_unique_reads)

#make different granges for the different amplicons
SMC6_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="SMC6",],keep.extra.columns = TRUE)
MYCN_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="MYCN",],keep.extra.columns = TRUE)
ODC1_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="ODC1",],keep.extra.columns = TRUE)
MDM2_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="MDM2",],keep.extra.columns = TRUE)
CDK4_grange<- makeGRangesFromDataFrame(reconstrucion_bed[reconstrucion_bed$circ_id=="CDK4",],keep.extra.columns = TRUE)

#6. first filter
MN_pass_vec1 <- meta_qc[as.numeric(meta_qc$n_unique_reads) >= 100000 & as.numeric(meta_qc$n_positive_windows_vec) >=20 & meta_qc$class == "MN",]

#7. subset window_calls_counts_normalised df and determine consensus sequence of circular DNA in MN
window_calls_counts_normalised_MN <- window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec1$sample,]
MN_consensus_calls_mat <- unstack(window_calls_counts_normalised_MN,window_call_final~samplename)
MN_consensus_calls_mat$rowsum <- rowSums(MN_consensus_calls_mat)
MN_consensus_calls_vec <- rownames(MN_consensus_calls_mat[MN_consensus_calls_mat$rowsum>=5,])
MN_consensus_regions<- window_calls_counts_normalised_MN[window_calls_counts_normalised_MN$samplename == "A1",][MN_consensus_calls_vec,c(1,2,3)] ##########
MN_consensus_regions <- MN_consensus_regions[MN_consensus_regions$chr!="chrM",] #remove chrM

#8. overlap identified consensus regions with reference
reconstrucion_bed_grange <- makeGRangesFromDataFrame(reconstrucion_bed,keep.extra.columns = TRUE)
consensus_regions_grange <- makeGRangesFromDataFrame(MN_consensus_regions) 

MNsequences_on_ecDNA <- subsetByOverlaps(consensus_regions_grange,reconstrucion_bed_grange) #extract indices of all regions in the consensus, which are overlapping the ecDNA reconstruction
index_MNsequences_on_ecDNA <-rownames(as.data.frame(MNsequences_on_ecDNA)) # for later quantitaive analysis

#9 add bin counts
#assign amplicons to each bin
amplicon_assignment_vec <- c()
for(i in seq_along(MNsequences_on_ecDNA)) {
  grange_MN <- MNsequences_on_ecDNA[i]
  current_gene <-  "nope"
  SMC6_sum <- sum(countOverlaps(SMC6_grange,grange_MN))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,grange_MN))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,grange_MN))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,grange_MN))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,grange_MN))
  
  if(SMC6_sum != 0){current_gene <- "SMC6"}
  if(MYCN_sum != 0){current_gene <- "MYCN"}
  if(ODC1_sum != 0){current_gene <- "ODC1"}
  if(MDM2_sum != 0){current_gene <- "MDM2"}
  if(CDK4_sum != 0){current_gene <- "CDK4"}
  
  amplicon_assignment_vec <- c(amplicon_assignment_vec,current_gene)
}

amplicon_index_df <- cbind(as.data.frame(MNsequences_on_ecDNA)[,1:3],amplicon_assignment_vec)
window_calls_counts_normalised_amplicon_anno <- window_calls_counts_normalised
A1_df <- window_calls_counts_normalised[window_calls_counts_normalised$samplename == "A1",]
A1_df$amplicon <- "z"
A1_df$amplicon[match(row.names(amplicon_index_df), row.names(A1_df))] <- amplicon_index_df$amplicon_assignment_vec
A1_df[row.names(A1_df) %in% row.names(amplicon_index_df), ]

window_calls_counts_normalised_amplicon_anno$amplicon_anno <- rep(A1_df$amplicon, length(unique(window_calls_counts_normalised$samplename)))

#bin counts
sample_count_vec = c()
bp_covered_vec_vec = c()
n_short_fragments_vec = c()
n_long_fragments_vec = c()
for(sample in unique(window_calls_counts_normalised_amplicon_anno$samplename)){
  window_calls_counts_normalised_amplicon_anno_sample_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename == sample & window_calls_counts_normalised_amplicon_anno$window_call_final == 1,]
  fragment_i <- 1
  current_index <- 0
  prev_index <- 0
  fragment_vec <- c()
  if (dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1] != 0){
    for (i in seq(1:dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1])){
      sample_row <- window_calls_counts_normalised_amplicon_anno_sample_df[i,]
      current_index <- as.numeric(row.names(sample_row))
      diff <- current_index - prev_index 
      if (i == 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff <= 3 & i != 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff > 3 & i != 1){
        fragment_i <- fragment_i + 1
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
    }
    
    window_calls_counts_normalised_amplicon_anno_sample_df$fragment <- fragment_vec
    bp_covered <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)) * 25000
    n_short_fragments <- sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment) == 3)
    n_long_fragments <- dim(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)) - sum(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment) == 3)
    
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,bp_covered)
    n_short_fragments_vec <- c(n_short_fragments_vec,n_short_fragments)
    n_long_fragments_vec <- c(n_long_fragments_vec,n_long_fragments)
    
  }
  else {
    sample_count_vec <- c(sample_count_vec,sample)
    bp_covered_vec_vec <- c(bp_covered_vec_vec,0)
    n_short_fragments_vec <- c(n_short_fragments_vec,0)
    n_long_fragments_vec <- c(n_long_fragments_vec,0)
  }
}
#long fragments >3 bins
bin_count_meta <- data.frame(sample = sample_count_vec, bp_covered = bp_covered_vec_vec, n_short_fragments = n_short_fragments_vec, n_long_fragments = n_long_fragments_vec)

meta_qc_ext <- merge(meta_qc,bin_count_meta, by  = "sample") #extended meta
meta_qc_ext[meta_qc_ext$class == "MN", ]

#9. second filter: filter samples based on exclusion cirteria (here: MN, <= 20 windows, <= 100000 reads, and exclude if only short fragments)
MN_pass_vec <- meta_qc_ext[as.numeric(meta_qc_ext$n_unique_reads) >= 100000 & as.numeric(meta_qc_ext$n_positive_windows_vec) >=20 & meta_qc_ext$class == "MN" & meta_qc_ext$n_long_fragments != 0,]
pMN_pass_vec <- meta_qc_ext[as.numeric(meta_qc_ext$n_unique_reads) >= 100000 & as.numeric(meta_qc_ext$n_positive_windows_vec) >=20 & meta_qc_ext$class == "mMN" & meta_qc_ext$n_long_fragments != 0,]
PN_pass_vec <- meta_qc[meta_qc$class == "PN" & as.numeric(meta_qc$n_positive_windows_vec) >= 50,]

###
##upsets


#make a vectors for different amplicons, to which positve sample will be added
SMC6 <- c()
MYCN <- c()
ODC1 <- c()
MDM2 <- c()
CDK4 <- c()

for(x in PN_pass_vec$sample){
  
  sample_coords <- window_calls[window_calls$samplename == x & window_calls$window_call_final > 0, ][1:3] #subset sample and positive windows
  #test for amplicon overlap
  sample_coord_grange <-  makeGRangesFromDataFrame(sample_coords)
  SMC6_sum <- sum(countOverlaps(SMC6_grange,sample_coord_grange))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,sample_coord_grange))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,sample_coord_grange))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,sample_coord_grange))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,sample_coord_grange))
  
  if(SMC6_sum > 0){
    SMC6 <- c(SMC6,x)}
  if(MYCN_sum > 0){
    MYCN <- c(MYCN,x)}
  if(ODC1_sum > 0){
    ODC1 <- c(ODC1,x)}
  if(MDM2_sum > 0){
    MDM2 <- c(MDM2,x)}
  if(CDK4_sum > 0){
    CDK4 <- c(MDM2,x)}
}

Amplicon_list <- list(MDM2 = MDM2, SMC6 = SMC6, MYCN = MYCN, ODC1 = ODC1, CDK4=CDK4 )
pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2//upset_4windows_PN_prov_final.pdf", width = 5, height = 5)

upset(fromList(Amplicon_list), order.by = "freq", mainbar.y.max = 12, empty.intersections = "on", sets.bar.color = "#8B0000", matrix.color = "#8B0000", main.bar.color = "#8B0000", nintersects = 5, keep.order = T, sets = c("MYCN","CDK4","MDM2","ODC1","SMC6"))
dev.off()

##MN Upset

SMC6 <- c()
MYCN <- c()
ODC1 <- c()
MDM2 <- c()
CDK4 <- c()

for(x in MN_pass_vec$sample){
  
  sample_coords <- window_calls[window_calls$samplename == x & window_calls$window_call_final > 0, ][1:3] #subset sample and positive windows
  #test for amplicon overlap
  sample_coord_grange <-  makeGRangesFromDataFrame(sample_coords)
  SMC6_sum <- sum(countOverlaps(SMC6_grange,sample_coord_grange))
  MYCN_sum <- sum(countOverlaps(MYCN_grange,sample_coord_grange))
  ODC1_sum <- sum(countOverlaps(ODC1_grange,sample_coord_grange))
  MDM2_sum <- sum(countOverlaps(MDM2_grange,sample_coord_grange))
  CDK4_sum <- sum(countOverlaps(CDK4_grange,sample_coord_grange))
  
  if(SMC6_sum > 0){
    SMC6 <- c(SMC6,x)}
  if(MYCN_sum > 0){
    MYCN <- c(MYCN,x)}
  if(ODC1_sum > 0){
    ODC1 <- c(ODC1,x)}
  if(MDM2_sum > 0){
    MDM2 <- c(MDM2,x)}
  if(CDK4_sum > 0){
    CDK4 <- c(MDM2,x)}
}

Amplicon_list <- list(MDM2 = MDM2, SMC6 = SMC6, MYCN = MYCN, ODC1 = ODC1, CDK4=CDK4 )
pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/upset_4windows_MN_prov_final.pdf", width = 5, height = 5 )

upset(fromList(Amplicon_list), order.by = "freq", mainbar.y.max = 12,empty.intersections = "on", nintersects = 5, keep.order = T)
dev.off()

##


##
#logCPM mean
#filter PN 
PN_pass_vec <- meta_qc[meta_qc$class == "PN" & as.numeric(meta_qc$n_positive_windows_vec) >= 50,]

#PN
PN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_called_region_mat_ecDNA <- t(PN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(PN_logCPM_called_region_mat_ecDNA)
mean(PN_logCPM_called_region_mat_ecDNA_mean)

PN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_lin_mat <- t(PN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_lin_mat_mean <- rowMeans(PN_logCPM_lin_mat)
mean(PN_logCPM_lin_mat_mean)


#MN
MN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA <- t(MN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(MN_logCPM_called_region_mat_ecDNA)
mean(MN_logCPM_called_region_mat_ecDNA_mean)

MN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_lin_mat <- t(MN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_lin_mat_mean <- rowMeans(MN_logCPM_lin_mat)
mean(MN_logCPM_lin_mat_mean)



#compute (circCPM/linCPM)
MN_CPM_ratio <- log2(as.vector(MN_logCPM_called_region_mat_ecDNA_mean)/MN_logCPM_lin_mat_mean)

PN_CPM_ratio <- log2(as.vector(PN_logCPM_called_region_mat_ecDNA_mean)/PN_logCPM_lin_mat_mean)


##combine MN and PN for fig
MN_CPM_ratio_df <- data.frame(MN_CPM_ratio)
colnames(MN_CPM_ratio_df) <- "CPM_ratio"
MN_CPM_ratio_df$sample <- names(MN_CPM_ratio)
MN_CPM_ratio_df$grp <- "Micronucleus"


PN_CPM_ratio_df <- data.frame(PN_CPM_ratio)
colnames(PN_CPM_ratio_df) <- "CPM_ratio"
PN_CPM_ratio_df$sample <- names(PN_CPM_ratio)
PN_CPM_ratio_df$grp <- "Primary Nucleus"

mean(MN_CPM_ratio_df$CPM_ratio)
mean(PN_CPM_ratio_df$CPM_ratio)

CPM_ratio_df <- rbind(MN_CPM_ratio_df,PN_CPM_ratio_df)
anno <- t.test(MN_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/plots_prov/circ_lin_CPM_ratio_box.pdf", width =6 , height = 6)

ggplot(data=CPM_ratio_df, aes(x=grp, y=CPM_ratio, fill=grp)) +
  geom_boxplot(width = 1,outlier.shape = NA) +
  ggtitle("") +
  xlab("") + ylab("log2 (circular/linear DNA mean(CPM))") +
  geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.1)+
  geom_signif(comparisons = list(c("Primary Nucleus", "Micronucleus")), 
              map_signif_level=TRUE, 
              annotation = formatC(anno, digits = 1), ) + theme_classic()
dev.off()



t.test(CPM_ratio~grp, data = CPM_ratio_df, 
       alternative = "two.sided")


##
#ecDNA ratios
#bulk
ODC1_ratio_bulk <- log2(456.548275527328/95.44642818592865) #3
MDM2_ratio_bulk <- log2(456.548275527328/121.75477264406456) #4
SMC6_ratio_bulk <- log2(456.548275527328/20.464356546540927) #1
CDK4_ratio_bulk <- log2(456.548275527328/70) #5
##

ODC1_ratio_vec <- c()
MDM2_ratio_vec <- c()
SMC6_ratio_vec <- c()
CDK4_ratio_vec <- c()
sample_name_vec <- c()

for(x in MN_pass_vec$sample){
  sample_df <- window_calls_counts_normalised[window_calls_counts_normalised$samplename == x,]
  sample_grange <- makeGRangesFromDataFrame(sample_df,keep.extra.columns = TRUE)
  #Use MYCN as anchor
  MYCN_anchor_grange <- subsetByOverlaps(sample_grange,MYCN_grange)
  
  SMC6_sample_grange <- subsetByOverlaps(sample_grange,SMC6_grange)
  SMC6_MYCN_ratio<- mean(MYCN_anchor_grange$CPM)/mean(SMC6_sample_grange$CPM)
  
  MDM2_sample_grange <- subsetByOverlaps(sample_grange,MDM2_grange)
  MDM2_MYCN_ratio<- mean(MYCN_anchor_grange$CPM)/mean(MDM2_sample_grange$CPM)
  
  ODC1_sample_grange <- subsetByOverlaps(sample_grange,ODC1_grange)
  ODC1_MYCN_ratio<- mean(MYCN_anchor_grange$CPM)/mean(ODC1_sample_grange$CPM)
  
  CDK4_sample_grange <- subsetByOverlaps(sample_grange,CDK4_grange)
  CDK4_MYCN_ratio<- mean(MYCN_anchor_grange$CPM)/mean(CDK4_sample_grange$CPM)
  
  ODC1_ratio_vec <- c(ODC1_ratio_vec,ODC1_ratio_bulk/log2(ODC1_MYCN_ratio))
  MDM2_ratio_vec <- c(MDM2_ratio_vec,MDM2_ratio_bulk/log2(MDM2_MYCN_ratio))
  SMC6_ratio_vec <- c(SMC6_ratio_vec,SMC6_ratio_bulk/log2(SMC6_MYCN_ratio))
  CDK4_ratio_vec <- c(CDK4_ratio_vec,CDK4_ratio_bulk/log2(CDK4_MYCN_ratio))
  
  sample_name_vec <- c(sample_name_vec,x)
  
}



median(ODC1_ratio_vec)
median(MDM2_ratio_vec)
median(SMC6_ratio_vec)

MN_amplicon_ratios <- data.frame(MDM2_ratio_vec,ODC1_ratio_vec,SMC6_ratio_vec,CDK4_ratio_vec,sample_name_vec)
MN_amplicon_ratios_df<- reshape(MN_amplicon_ratios, 
                                direction = "long",
                                varying = list(names(MN_amplicon_ratios)[1:4]),
                                v.names = "CPM",
                                idvar = c("sample_name_vec"),
                                timevar = "Amplicon",
                                times = colnames(MN_amplicon_ratios)[1:4])

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/scratch/ratios_box_ratio_bulk.pdf", width =6 , height = 6)

ggplot(data=MN_amplicon_ratios_df, aes(x=Amplicon, y=CPM, fill=Amplicon)) +
  geom_boxplot(width = 1) +
  ggtitle("") +
  xlab("") + ylab("log2(mean(CPM) ratio)") +
  geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.1) + ylim(-10, 10)

#+
  #geom_signif(comparisons = list(c("MN circular DNA", "MN linear DNA")), 
    #          map_signif_level=TRUE, 
   #           test = "wilcox.test") + theme_classic()
dev.off()



###heatmap version 2
library(pheatmap)
library(ComplexHeatmap)
#make seperate heatmaps and apply row clustering (hclust)
##MN

##
MN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% MN_pass_vec$sample,]

MN_mat_heatmap <- unstack(MN_mat_heatmap_df,log2CPM~samplename)
MN_mat_heatmap_sub <- MN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]

window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename == "A1",][c(10000,10700,80000,80600),]




col_anno <- data.frame(MN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- MN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/heatmapt_MN.pdf", width =10 , height = 3)

heatmap_MN <- pheatmap(data.frame(t(MN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_MN
dev.off()

MN_heatmap_row_order <- rownames(data.frame(t(MN_mat_heatmap_sub))[row_order(heatmap_MN),])


###PN


PN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% PN_pass_vec$sample,]

PN_mat_heatmap <- unstack(PN_mat_heatmap_df,log2CPM~samplename)
PN_mat_heatmap_sub <- PN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


col_anno <- data.frame(PN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- PN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/heatmapt_PN.pdf", width =10 , height = 3)

heatmap_PN <- pheatmap(as.matrix(t(PN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_PN
dev.off()

PN_heatmap_row_order <- rownames(data.frame(t(PN_mat_heatmap_sub))[row_order(heatmap_PN),])



###pMN


pMN_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% pMN_pass_vec$sample,]

pMN_mat_heatmap <- unstack(pMN_mat_heatmap_df,log2CPM~samplename)
pMN_mat_heatmap_sub <- pMN_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


col_anno <- data.frame(pMN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- pMN_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/heatmapt_pMN.pdf", width =10 , height = 3)

heatmap_pMN <- pheatmap(as.matrix(t(pMN_mat_heatmap_sub)), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = TRUE,gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_pMN
dev.off()

pMN_heatmap_row_order <- rownames(data.frame(t(pMN_mat_heatmap_sub))[row_order(heatmap_pMN),])

##Master Heatmap

heatmap_order <- c(pMN_heatmap_row_order,MN_heatmap_row_order,PN_heatmap_row_order)
breaks_rows <- c(length(pMN_heatmap_row_order),length(pMN_heatmap_row_order)+length(MN_heatmap_row_order))
Master_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% heatmap_order,]

Master_mat_heatmap <- unstack(Master_mat_heatmap_df,log2CPM~samplename)
Master_mat_heatmap_sub <- Master_mat_heatmap[c(seq(10000,10700),seq(80000,80600)),]


#row anno
row_anno <- meta_qc_ext[meta_qc_ext$sample %in% heatmap_order, ]
rownames(row_anno) <- row_anno$sample
row_anno$placeholder <- "a"
row_anno_sub <- row_anno[heatmap_order,c("class","placeholder")]

col_anno <- data.frame(Master_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$chr)
row.names(col_anno) <- c(seq(10000,10700),seq(80000,80600))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- Master_mat_heatmap_df[c(seq(10000,10700),seq(80000,80600)),]$amplicon_anno
count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/heatmap_main_v2.pdf", width =10 , height = 3)

heatmap_main <- pheatmap(as.matrix(t(Master_mat_heatmap_sub[,heatmap_order])), color=colorRampPalette(c("white","red","darkred"))(200), annotation_colors = ann_colors,cluster_rows = F,gaps_row = breaks_rows, gaps_col = 701 ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_row = row_anno_sub,annotation_col = col_anno[,c("Amplicon","Region")])
heatmap_main
dev.off()


#make amplicon specific annotations
Master_mat_heatmap_df[Master_mat_heatmap_df$samplename=="A1",][seq(10044,10064),]
#MYCN

amplicon_vec <- c("MYCN","MYCN","MYCN","MYCN","CDK4","CDK4","SMC6","ODC1","ODC1","MDM2")
amplicon_start <- c(10052,10095,10506,10592,80070,82170,10638,10366,80041,80509)
amplicon_end <-c(10056,10097,10508,10613,80076,82174,10677,10389,80045,80549)

amplicon_anno_df <- data.frame(amplicon = amplicon_vec, start = amplicon_start, end = amplicon_end )


#extend 8 bins
amplicon_anno_df$start <- amplicon_anno_df$start - 4
amplicon_anno_df$end <- amplicon_anno_df$end + 4
#make annotation

seq_vec <- c()
region_vec <- c()
sub_temp <- Master_mat_heatmap_df[Master_mat_heatmap_df$samplename=="A1",]
for (i in 1:nrow(amplicon_anno_df)){
  seq_vec <- c(seq_vec, seq(amplicon_anno_df[i,"start"],amplicon_anno_df[i,"end"]))
  region_vec <- c(region_vec,rep(paste(format(round(sub_temp[amplicon_anno_df[i,"start"],"start"] / 1e6, 1), trim = TRUE),"Mb",format(round(sub_temp[amplicon_anno_df[i,"end"],"end"] / 1e6, 1), trim = TRUE),"Mb"),length(seq(amplicon_anno_df[i,"start"],amplicon_anno_df[i,"end"]))))
}

amplicon_extend_anno_df <- data.frame(index = seq_vec, region_vec = region_vec)
rownames(amplicon_extend_anno_df) <- amplicon_extend_anno_df$index

##heatmap_reduced_ final

##Master Heatmap

heatmap_order <- c(pMN_heatmap_row_order,MN_heatmap_row_order,PN_heatmap_row_order)
breaks_rows <- c(length(pMN_heatmap_row_order),length(pMN_heatmap_row_order)+length(MN_heatmap_row_order))

Master_mat_heatmap_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% heatmap_order,]

Master_mat_heatmap <- unstack(Master_mat_heatmap_df,log2CPM~samplename)
Master_mat_heatmap_sub <- Master_mat_heatmap[sort(amplicon_extend_anno_df$index),]


#row anno
row_anno <- meta_qc_ext[meta_qc_ext$sample %in% heatmap_order, ]
rownames(row_anno) <- row_anno$sample
row_anno$placeholder <- "a"
row_anno_sub <- row_anno[heatmap_order,c("class","placeholder")]

col_anno <- data.frame(Master_mat_heatmap_df[sort(amplicon_extend_anno_df$index),]$chr)
row.names(col_anno) <- c(sort(amplicon_extend_anno_df$index))
colnames(col_anno) <- "Region"
col_anno$Amplicon <- Master_mat_heatmap_df[sort(amplicon_extend_anno_df$index),]$amplicon_anno
col_anno$Region_coord <- amplicon_extend_anno_df[paste(sort(amplicon_extend_anno_df$index)),]$region_vec
count(col_anno$Region) #determine number of chr bins
#colour annotation
ann_colors = list(
  Amplicon = c(z = "white", MYCN = "#e41a1c", CDK4 = "#377eb8", ODC1 = "#4daf4a", MDM2 = "#984ea3", SMC6 = "#ff7f00"),
  Region = c(chr2 = "grey", chr12 = "black")
)


pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/heatmap_main_cleaned_new_v1_color_BYR.pdf", width =10 , height = 3)

heatmap_main <- pheatmap(as.matrix(t(Master_mat_heatmap_sub[,heatmap_order])), color=colorRampPalette(c("#2c7bb6","#ffffbf","#d7191c"))(200), annotation_colors = ann_colors,cluster_rows = F,gaps_row = breaks_rows,  gaps_col = c(21-8,40-(2*8),80-(3*8),99-(4*8),137-(5*8),193-(6*8),214-(7*8),237-(8*8),294-(9*8)) ,cluster_cols = FALSE, show_rownames = T, show_colnames = FALSE, annotation_row = row_anno_sub,annotation_col = col_anno[,c("Amplicon","Region","Region_coord")])
heatmap_main
dev.off()

table(col_anno$Region_coord)

gaps_col = c(21,40,80,99,137,193,214,237,294)+1
table(col_anno$Region_coord)



##copy number:
#uses winsorized mean for lin
library("DescTools")

#logCPM mean


#PN
PN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_called_region_mat_ecDNA <- t(PN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(PN_logCPM_called_region_mat_ecDNA)
mean(PN_logCPM_called_region_mat_ecDNA_mean)

PN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_lin_mat <- t(PN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(PN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))

mean(PN_logCPM_lin_mat_mean_wins)


#MN
MN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA <- t(MN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(MN_logCPM_called_region_mat_ecDNA)
mean(MN_logCPM_called_region_mat_ecDNA_mean)

MN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_lin_mat <- t(MN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(MN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))

mean(MN_logCPM_lin_mat_mean_wins)

#amplicon specific copy number

window_calls_counts_normalised_amplicon_anno

MN_logCPM_called_region_mat_SMC6 <- unstack(window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% MN_pass_vec$sample & window_calls_counts_normalised_amplicon_anno$amplicon_anno == "SMC6",],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA_SMC6 <- t(MN_logCPM_called_region_mat_SMC6)
MN_logCPM_called_region_mat_ecDNA_mean_SMC6 <- rowMeans(MN_logCPM_called_region_mat_ecDNA_SMC6)
mean(MN_logCPM_called_region_mat_ecDNA_mean_SMC6)


MN_logCPM_called_region_mat_MYCN <- unstack(window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% MN_pass_vec$sample & window_calls_counts_normalised_amplicon_anno$amplicon_anno == "MYCN",],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA_MYCN <- t(MN_logCPM_called_region_mat_MYCN)
MN_logCPM_called_region_mat_ecDNA_mean_MYCN <- rowMeans(MN_logCPM_called_region_mat_ecDNA_MYCN)
mean(MN_logCPM_called_region_mat_ecDNA_mean_MYCN)

log2(MN_logCPM_called_region_mat_ecDNA_mean_MYCN/mean(MN_logCPM_lin_mat_mean_wins))
ll1 <- scale(log2(MN_logCPM_called_region_mat_ecDNA_mean_MYCN/mean(MN_logCPM_lin_mat_mean_wins)), center = TRUE, scale = TRUE)



PN_logCPM_called_region_mat_MYCN <- unstack(window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% PN_pass_vec$sample & window_calls_counts_normalised_amplicon_anno$amplicon_anno == "MYCN",],CPM~samplename) #CPM matrix
PN_logCPM_called_region_mat_ecDNA_MYCN <- t(PN_logCPM_called_region_mat_MYCN)
PN_logCPM_called_region_mat_ecDNA_mean_MYCN <- rowMeans(PN_logCPM_called_region_mat_ecDNA_MYCN)
mean(PN_logCPM_called_region_mat_ecDNA_mean_MYCN)

log2(PN_logCPM_called_region_mat_ecDNA_mean_MYCN/mean(PN_logCPM_lin_mat_mean_wins))
ll2 <- scale(log2(PN_logCPM_called_region_mat_ecDNA_mean_MYCN/mean(PN_logCPM_lin_mat_mean_wins)), center = TRUE, scale = TRUE)

as.numeric(ll2)


ks.test(ll1,ll2,alternative = "two.sided",
        exact = NULL)


##perc amplicon occupancy


###PN
sample_vec <- c()
MYCN_frac_vec <- c()
MDM2_frac_vec <- c()
CDK4_frac_vec <- c()
ODC1_frac_vec <- c()
SMC6_frac_vec <- c()

for(sample in PN_pass_vec$sample){
  read_count_df <- window_calls_counts_normalised_amplicon_anno
  MYCN_length_norm <- sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$amplicon_anno)
  MDM2_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$amplicon_anno)
  CDK4_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$amplicon_anno)
  ODC1_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$amplicon_anno)
  SMC6_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$amplicon_anno)
  
  length_norm_sum <- sum(c(MYCN_length_norm,MDM2_length_norm,CDK4_length_norm,ODC1_length_norm,SMC6_length_norm))
  
  sample_vec <- c(sample_vec,sample)
  MYCN_frac_vec <- c(MYCN_frac_vec,(MYCN_length_norm / length_norm_sum))
  MDM2_frac_vec <- c(MDM2_frac_vec,(MDM2_length_norm / length_norm_sum))
  CDK4_frac_vec <- c(CDK4_frac_vec,(CDK4_length_norm / length_norm_sum))
  ODC1_frac_vec <- c(ODC1_frac_vec,(ODC1_length_norm / length_norm_sum))
  SMC6_frac_vec <- c(SMC6_frac_vec,(SMC6_length_norm / length_norm_sum))
  
}
gene <- c(rep("MYCN",length(sample_vec)),rep("MDM2",length(sample_vec)),rep("CDK4",length(sample_vec)),rep("ODC1",length(sample_vec)),rep("SMC6",length(sample_vec)))
fraction <- c(MYCN_frac_vec,MDM2_frac_vec,CDK4_frac_vec,ODC1_frac_vec,SMC6_frac_vec)
PN_fraction_df <- data.frame(gene,fraction)
PN_fraction_df$grp <- "Primary Nucleus"


###MN
sample_vec <- c()
MYCN_frac_vec <- c()
MDM2_frac_vec <- c()
CDK4_frac_vec <- c()
ODC1_frac_vec <- c()
SMC6_frac_vec <- c()

for(sample in MN_pass_vec$sample){
  read_count_df <- window_calls_counts_normalised_amplicon_anno
  MYCN_length_norm <- sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MYCN",]$amplicon_anno)
  MDM2_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "MDM2",]$amplicon_anno)
  CDK4_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "CDK4",]$amplicon_anno)
  ODC1_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "ODC1",]$amplicon_anno)
  SMC6_length_norm <-sum(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$nReads)/length(read_count_df[read_count_df$samplename == sample & read_count_df$amplicon_anno == "SMC6",]$amplicon_anno)
  
  length_norm_sum <- sum(c(MYCN_length_norm,MDM2_length_norm,CDK4_length_norm,ODC1_length_norm,SMC6_length_norm))
  
  sample_vec <- c(sample_vec,sample)
  MYCN_frac_vec <- c(MYCN_frac_vec,(MYCN_length_norm / length_norm_sum))
  MDM2_frac_vec <- c(MDM2_frac_vec,(MDM2_length_norm / length_norm_sum))
  CDK4_frac_vec <- c(CDK4_frac_vec,(CDK4_length_norm / length_norm_sum))
  ODC1_frac_vec <- c(ODC1_frac_vec,(ODC1_length_norm / length_norm_sum))
  SMC6_frac_vec <- c(SMC6_frac_vec,(SMC6_length_norm / length_norm_sum))
  
}

gene <- c(rep("MYCN",length(sample_vec)),rep("MDM2",length(sample_vec)),rep("CDK4",length(sample_vec)),rep("ODC1",length(sample_vec)),rep("SMC6",length(sample_vec)))
fraction <- c(MYCN_frac_vec,MDM2_frac_vec,CDK4_frac_vec,ODC1_frac_vec,SMC6_frac_vec)
MN_fraction_df <- data.frame(gene,fraction)
MN_fraction_df$grp <- "Micronucleus"


fraction_df <- rbind(MN_fraction_df,PN_fraction_df)

fraction_df$log2fraction <- log2(fraction_df$fraction+1)



##plot
pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/fractions_box.pdf", width =10 , height = 6)

ggplot(data=fraction_df, aes(x=gene, y=fraction, fill=grp)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("") +
  xlab("") + ylab("normalized amplicon fraction") +
 #geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.1) + 
  geom_point(color="black", size=0.3, alpha=0.9,position=position_jitterdodge(jitter.width = 0.1))+
  theme_classic() + 
  stat_compare_means(aes(group = grp), method = "t.test",
                                       method.args = list(var.equal = F, alternative = "two.sided"),label = "p.signif")


dev.off()

saveRDS(fraction_df, file = "/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/L_fractions_box.rds")


####winsorised
#PN
PN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_called_region_mat_ecDNA <- t(PN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(PN_logCPM_called_region_mat_ecDNA)
mean(PN_logCPM_called_region_mat_ecDNA_mean)

PN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% PN_pass_vec$sample ,],CPM~samplename) #CPM matrix
PN_logCPM_lin_mat <- t(PN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
PN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(PN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))

mean(PN_logCPM_lin_mat_mean_wins)


#MN
MN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_called_region_mat_ecDNA <- t(MN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(MN_logCPM_called_region_mat_ecDNA)
mean(MN_logCPM_called_region_mat_ecDNA_mean)

MN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% MN_pass_vec$sample ,],CPM~samplename) #CPM matrix
MN_logCPM_lin_mat <- t(MN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
MN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(MN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))


#pooled MN
pMN_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% pMN_pass_vec$sample ,],CPM~samplename) #CPM matrix
pMN_logCPM_called_region_mat_ecDNA <- t(pMN_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)]
pMN_logCPM_called_region_mat_ecDNA_mean <- rowMeans(pMN_logCPM_called_region_mat_ecDNA)
mean(pMN_logCPM_called_region_mat_ecDNA_mean)

pMN_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% pMN_pass_vec$sample ,],CPM~samplename) #CPM matrix
pMN_logCPM_lin_mat <- t(pMN_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
pMN_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(pMN_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))


#empty
empty_logCPM_called_region_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% meta[meta$class == "Empty",]$sample ,],CPM~samplename) #CPM matrix
empty_logCPM_called_region_mat_ecDNA <- t(empty_logCPM_called_region_mat)[,as.numeric(index_MNsequences_on_ecDNA)] 
empty_logCPM_called_region_mat_ecDNA_mean <- rowMeans(empty_logCPM_called_region_mat_ecDNA)
mean(empty_logCPM_called_region_mat_ecDNA_mean)

empty_logCPM_lin_mat <- unstack(window_calls_counts_normalised[window_calls_counts_normalised$samplename %in% meta[meta$class == "Empty",]$sample ,],CPM~samplename) #CPM matrix
empty_logCPM_lin_mat <- t(empty_logCPM_lin_mat)[,-as.numeric(index_MNsequences_on_ecDNA)]
empty_logCPM_lin_mat_mean_wins <- rowMeans(t(apply(empty_logCPM_lin_mat, 1, function(x) Winsorize(x, val = quantile(x, probs = c(0.05, 0.95))))))


#compute (circCPM/linCPM)
MN_CPM_ratio <- log2(as.vector(MN_logCPM_called_region_mat_ecDNA_mean)/MN_logCPM_lin_mat_mean_wins)

PN_CPM_ratio <- log2(as.vector(PN_logCPM_called_region_mat_ecDNA_mean)/PN_logCPM_lin_mat_mean_wins)

pMN_CPM_ratio <- log2(as.vector(pMN_logCPM_called_region_mat_ecDNA_mean)/pMN_logCPM_lin_mat_mean_wins)

empty_CPM_ratio <- log2(as.vector(empty_logCPM_called_region_mat_ecDNA_mean)/empty_logCPM_lin_mat_mean_wins)

##combine MN and PN for fig
MN_CPM_ratio_df <- data.frame(MN_CPM_ratio)
colnames(MN_CPM_ratio_df) <- "CPM_ratio"
MN_CPM_ratio_df$sample <- names(MN_CPM_ratio)
MN_CPM_ratio_df$grp <- "single micronucleus"


PN_CPM_ratio_df <- data.frame(PN_CPM_ratio)
colnames(PN_CPM_ratio_df) <- "CPM_ratio"
PN_CPM_ratio_df$sample <- names(PN_CPM_ratio)
PN_CPM_ratio_df$grp <- "Primary Nucleus"

pMN_CPM_ratio_df <- data.frame(pMN_CPM_ratio)
colnames(pMN_CPM_ratio_df) <- "CPM_ratio"
pMN_CPM_ratio_df$sample <- names(pMN_CPM_ratio)
pMN_CPM_ratio_df$grp <- "pooled micronuclei"

mean(MN_CPM_ratio_df$CPM_ratio)
mean(PN_CPM_ratio_df$CPM_ratio)
mean(pMN_CPM_ratio_df$CPM_ratio)

CPM_ratio_df <- rbind(pMN_CPM_ratio_df,MN_CPM_ratio_df,PN_CPM_ratio_df)

#empty
empty_CPM_ratio_df <- data.frame(empty_CPM_ratio)
colnames(empty_CPM_ratio_df) <- "CPM_ratio"
empty_CPM_ratio_df$sample <- names(empty_CPM_ratio)
empty_CPM_ratio_df$grp <- "empty"
CPM_ratio_df_empty <- rbind(pMN_CPM_ratio_df,MN_CPM_ratio_df,PN_CPM_ratio_df,empty_CPM_ratio_df)



#testing

anno <- t.test(MN_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno1 <- t.test(MN_CPM_ratio, pMN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno2 <- t.test(pMN_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
test_vec <- c(anno,anno1,anno2)

pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/circ_lin_CPM_ratio_box_winsorised.pdf", width =8 , height = 6)

ggplot(data=CPM_ratio_df, aes(x=factor(grp, level=c("pooled micronuclei","single micronucleus","Primary Nucleus")), y=CPM_ratio, fill=grp)) +
  geom_boxplot(width = 1,outlier.shape = NA) +
  ggtitle("") +
  xlab("") + ylab("log2 (circular/linear DNA mean(CPM))") +
  geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.1)+
  geom_signif(comparisons = list(c("Primary Nucleus", "single micronucleus"),c("single micronucleus", "pooled micronuclei"),c("pooled micronuclei", "Primary Nucleus")), 
              map_signif_level=TRUE, 
              annotation = formatC(test_vec, digits = 1),step_increase = 0.05 ) + theme_classic()
dev.off()





saveRDS(CPM_ratio_df, file = "/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/K_circ_lin_CPM_ratio_box_winsorised_raw_data.rds")

saveRDS(test_vec, file = "/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/K_circ_lin_CPM_ratio_box_winsorised_raw_welch_test_p.rds")


##for reviewers MNseq method description
anno <- t.test(MN_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno1 <- t.test(MN_CPM_ratio, pMN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno2 <- t.test(pMN_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno3a <- t.test(empty_CPM_ratio, pMN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno3b <- t.test(empty_CPM_ratio, PN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
anno3c <- t.test(empty_CPM_ratio, MN_CPM_ratio, alternative = "two.sided", var.equal = FALSE)$p.value
test_vec <- c(anno,anno1,anno2,anno3a,anno3b,anno3c)



pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/circ_lin_CPM_ratio_box_winsorised.pdf", width =8 , height = 6)

ggplot(data=CPM_ratio_df_empty, aes(x=factor(grp, level=c("pooled micronuclei","single micronucleus","Primary Nucleus","empty")), y=CPM_ratio, fill=grp)) +
  geom_boxplot(width = 1,outlier.shape = NA) +
  ggtitle("") +
  xlab("") + ylab("log2 (circular/linear DNA mean(CPM))") +
  geom_jitter(color="black", size=0.3, alpha=0.9, width = 0.1)+
  geom_signif(comparisons = list(c("Primary Nucleus", "single micronucleus"),c("single micronucleus", "pooled micronuclei"),c("pooled micronuclei", "Primary Nucleus"),c("empty", "pooled micronuclei"),c("empty", "Primary Nucleus"),c("empty", "single micronucleus")), 
              map_signif_level=TRUE, 
              annotation = formatC(test_vec, digits = 1),step_increase = 0.05 ) + theme_classic()
dev.off()



###SIGNAL Plot samples
PN_compute_mat_vec <- c()
MN_compute_mat_vec <- c()
pMN_compute_mat_vec <- c()

PN_samples_signal_sample_vec <- PN_pass_vec$sample
MN_samples_signal_sample_vec <- MN_pass_vec$sample
pMN_samples_signal_sample_vec <- pMN_pass_vec$sample


length(PN_samples_signal_sample_vec)
length(MN_samples_signal_sample_vec)
length(pMN_samples_signal_sample_vec)

for (i in 1:length(PN_samples_signal_sample_vec)){PN_compute_mat_vec <- c(PN_compute_mat_vec,paste("/data/cephfs-1/work/groups/henssen/users/xuro_c/sMNseq/pilot/samples/*",PN_samples_signal_sample_vec[i],"*/coverage/*",PN_samples_signal_sample_vec[i],"*25bp_CPM.bw", sep = ""))}
for (i in 1:length(MN_samples_signal_sample_vec)){MN_compute_mat_vec <- c(MN_compute_mat_vec,paste("/data/cephfs-1/work/groups/henssen/users/xuro_c/sMNseq/pilot/samples/*",MN_samples_signal_sample_vec[i],"*/coverage/*",MN_samples_signal_sample_vec[i],"*25bp_CPM.bw", sep = ""))}
for (i in 1:length(pMN_samples_signal_sample_vec)){pMN_compute_mat_vec <- c(pMN_compute_mat_vec,paste("/data/cephfs-1/work/groups/henssen/users/xuro_c/sMNseq/pilot/samples/*",pMN_samples_signal_sample_vec[i],"*/coverage/*",pMN_samples_signal_sample_vec[i],"*25bp_CPM.bw", sep = ""))}


paste(as.character(PN_compute_mat_vec),collapse=" ")
paste(as.character(MN_compute_mat_vec),collapse=" ")
paste(as.character(pMN_compute_mat_vec),collapse=" ")

#####Suppl. Fig

#called regions MN




##########

samples_region_call_length <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename %in% MN_pass_vec$sample,]

sample_count_vec = c()
bp_covered_vec_vec = c()
n_short_fragments_vec = c()
n_long_fragments_vec = c()
window_size_vec = c()
window_size_vec_kb = c()
for(sample in unique(samples_region_call_length$samplename)){
  window_calls_counts_normalised_amplicon_anno_sample_df <- window_calls_counts_normalised_amplicon_anno[window_calls_counts_normalised_amplicon_anno$samplename == sample & window_calls_counts_normalised_amplicon_anno$window_call_final == 1 & window_calls_counts_normalised_amplicon_anno$amplicon_anno == "z",]
  fragment_i <- 1
  current_index <- 0
  prev_index <- 0
  fragment_vec <- c()
  if (dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1] != 0){
    for (i in seq(1:dim(window_calls_counts_normalised_amplicon_anno_sample_df)[1])){
      sample_row <- window_calls_counts_normalised_amplicon_anno_sample_df[i,]
      current_index <- as.numeric(row.names(sample_row))
      diff <- current_index - prev_index 
      if (i == 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff <= 4 & i != 1){
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
      if (diff > 4 & i != 1){
        fragment_i <- fragment_i + 1
        fragment_vec <- c(fragment_vec,fragment_i)
        prev_index <- current_index
      }
    }
    
    window_calls_counts_normalised_amplicon_anno_sample_df$fragment <- fragment_vec
    
    if(sample == "B5"){print(window_calls_counts_normalised_amplicon_anno_sample_df)}
    
    window_size_vec = c(window_size_vec,as.numeric(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)))
#    print(as.numeric(table(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)))
#    print(sample)
    for(fragment in unique(window_calls_counts_normalised_amplicon_anno_sample_df$fragment)){
      frag_df <- window_calls_counts_normalised_amplicon_anno_sample_df[window_calls_counts_normalised_amplicon_anno_sample_df$fragment == fragment,]
      fragment_size <- dim(window_calls_counts_normalised_amplicon_anno_sample_df[window_calls_counts_normalised_amplicon_anno_sample_df$fragment == fragment,])[1]
     # print(fragment_size)
      if(fragment_size != 1){
        size_frag_kb<- (as.numeric(rownames(frag_df[dim(frag_df)[1],])) - as.numeric(rownames(frag_df[1,]))) *25
        window_size_vec_kb = c(window_size_vec_kb,size_frag_kb)
       # print(window_size_vec_kb)
        }
      else{
        window_size_vec_kb = c(window_size_vec_kb,25)
       # print(window_size_vec_kb)
        }
      
    }
    
  }
  
}


pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/figures_version2/barlot_suppl_length_nonecDNA_calls.pdf", width =6 , height = 6)

ggplot(data.frame(table(window_size_vec_kb)), aes(x=window_size_vec_kb, y=Freq)) + 
  geom_bar(stat = "identity") + theme_classic()

dev.off()


#####################################
########## MNseq benchmark for reviewers
#####################################

#plot1: mean coverage over cir vs lin in empties
recon_all_amplicon<- makeGRangesFromDataFrame(reconstrucion_bed,keep.extra.columns = TRUE) #reconstruction
ov_vec <- countOverlaps(makeGRangesFromDataFrame(window_calls_counts_normalised[,1:3]),recon_all_amplicon)
window_calls_counts_normalised_benchmark <- window_calls_counts_normalised
window_calls_counts_normalised_benchmark$amplicon <- ov_vec
#empties:
meta[meta$class == "Empty",]$sample

#linear cpm - empties

mean_CPM_benchmark_circ_empty <- c()
mean_CPM_benchmark_lin_empty <- c()
for(sample in meta[meta$class == "Empty",]$sample){
  mean_CPM_benchmark_circ <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon != 0,]$CPM)
  mean_CPM_benchmark_lin <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon == 0,]$CPM)
  mean_CPM_benchmark_circ_empty <- c(mean_CPM_benchmark_circ_empty,mean_CPM_benchmark_circ)
  mean_CPM_benchmark_lin_empty <- c(mean_CPM_benchmark_lin_empty,mean_CPM_benchmark_lin)
  
  }
benchmark_circ_empty <-data.frame(CPM = mean_CPM_benchmark_circ_empty, type= "circ" ,class ="empty")
benchmark_lin_empty <-data.frame(CPM = mean_CPM_benchmark_lin_empty, type= "lin" ,class ="empty")



#linear cpm - PN

mean_CPM_benchmark_circ_PN <- c()
mean_CPM_benchmark_lin_PN <- c()
for(sample in meta[meta$class == "PN",]$sample){
  mean_CPM_benchmark_circ <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon != 0,]$CPM)
  mean_CPM_benchmark_lin <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon == 0,]$CPM)
  mean_CPM_benchmark_circ_PN <- c(mean_CPM_benchmark_circ_PN,mean_CPM_benchmark_circ)
  mean_CPM_benchmark_lin_PN <- c(mean_CPM_benchmark_lin_PN,mean_CPM_benchmark_lin)
  
}

benchmark_circ_PN <-data.frame(CPM = mean_CPM_benchmark_circ_PN, type= "circ" ,class ="PN")
benchmark_lin_PN <-data.frame(CPM = mean_CPM_benchmark_lin_PN, type= "lin" ,class ="PN")



#linear cpm - MN

mean_CPM_benchmark_circ_MN <- c()
mean_CPM_benchmark_lin_MN <- c()
for(sample in meta[meta$class == "MN",]$sample){
  mean_CPM_benchmark_circ <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon != 0,]$CPM)
  mean_CPM_benchmark_lin <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon == 0,]$CPM)
  mean_CPM_benchmark_circ_MN <- c(mean_CPM_benchmark_circ_MN,mean_CPM_benchmark_circ)
  mean_CPM_benchmark_lin_MN <- c(mean_CPM_benchmark_lin_MN,mean_CPM_benchmark_lin)
  
}

benchmark_circ_MN <- data.frame(CPM = mean_CPM_benchmark_circ_MN, type= "circ" ,class ="MN")
benchmark_lin_MN <- data.frame(CPM = mean_CPM_benchmark_lin_MN, type= "lin" ,class ="MN")

#linear cpm - mMN

mean_CPM_benchmark_circ_mMN <- c()
mean_CPM_benchmark_lin_mMN <- c()
for(sample in meta[meta$class == "mMN",]$sample){
  mean_CPM_benchmark_circ <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon != 0,]$CPM)
  mean_CPM_benchmark_lin <- mean(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample & window_calls_counts_normalised_benchmark$amplicon == 0,]$CPM)
  mean_CPM_benchmark_circ_mMN <- c(mean_CPM_benchmark_circ_mMN,mean_CPM_benchmark_circ)
  mean_CPM_benchmark_lin_mMN <- c(mean_CPM_benchmark_lin_mMN,mean_CPM_benchmark_lin)
  
}

benchmark_circ_mMN <-data.frame(CPM = mean_CPM_benchmark_circ_mMN, type= "circ" ,class ="mMN")
benchmark_lin_mMN <-data.frame(CPM = mean_CPM_benchmark_lin_mMN, type= "lin" ,class ="mMN")

benchmark_fig1 <- rbind(benchmark_circ_empty,benchmark_lin_empty,benchmark_circ_PN,benchmark_lin_PN,benchmark_circ_MN,benchmark_lin_MN,benchmark_circ_mMN,benchmark_lin_mMN)



pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/circ_lin_CPM_mean.pdf", width =8 , height = 6)

ggplot(data=benchmark_fig1, aes(x=factor(class, level=c("mMN","MN","PN","empty")), y=(CPM), fill=type)) +
  geom_boxplot(width = 1,outlier.shape = NA) +
  ggtitle("") +
  xlab("") + ylab("CPM")  + theme_classic()  + facet_wrap(~factor(class, levels=c("mMN","MN","PN","Empty")), scales = "free_x", nrow = 1)#+ scale_y_log10()
dev.off()


####
#####mappings to chr


idxstats_files <- list.files(path="/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings", pattern="*.idxstats", full.names=F, recursive=FALSE)
per_chromosome_chr_size <- read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V2

chr_map_df_list <- list()

for (file in idxstats_files){
  chr_map_df <- data.frame(chr_map = read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V3)
  
  per_chr_nreads <- read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V3
  library_size <- sum(per_chr_nreads)
  per_chr_reads_size_norm <- per_chr_nreads/per_chromosome_chr_size
  per_chr_perc <- per_chr_reads_size_norm/sum(per_chr_nreads/per_chromosome_chr_size)
  
  chr_map_df <- as.data.frame(per_chr_perc)
  colnames(chr_map_df) <-strsplit(file, split = "_")[[1]][3]
  chr_map_df_list <- c(chr_map_df_list,chr_map_df)
  }

per_chr_reads_df <- t(data.frame(chr_map_df_list))


table(meta$class)
per_chr_perc_order <- c(meta[meta$class=="mMN",][,2:3]$sample,meta[meta$class=="MN",][,2:3]$sample,meta[meta$class=="PN",][,2:3]$sample,meta[meta$class=="Empty",][,2:3]$sample)
per_chr_perc_order_anno <- c(meta[meta$class=="mMN",][,2:3]$class,meta[meta$class=="MN",][,2:3]$class,meta[meta$class=="PN",][,2:3]$class,meta[meta$class=="Empty",][,2:3]$class)




######

##Per_chr_perc  Heatmap



per_chr_perc_order_heatmap_anno <- data.frame(sample = per_chr_perc_order, class = per_chr_perc_order_anno)




pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/per_chr_heatmap.pdf", width =10 , height = 10)

heatmap_main_perchr <- pheatmap(as.matrix(per_chr_reads_df[per_chr_perc_order,]*100), color=colorRampPalette(c("white","red","darkred"))(100),cluster_cols = FALSE,cluster_rows=F, show_rownames = T, show_colnames = FALSE, annotation_row = per_chr_perc_order_heatmap_anno)
heatmap_main_perchr
dev.off()




####NON NORMALIZED
#####mappings to chr


idxstats_files <- list.files(path="/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings", pattern="*.idxstats", full.names=F, recursive=FALSE)
per_chromosome_chr_size <- read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V2

chr_map_df_list <- list()

for (file in idxstats_files){
  chr_map_df <- data.frame(chr_map = read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V3)
  
  per_chr_nreads <- read.table(paste("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/chromosome_mappings/",file,sep=""), header=F)[1:24,]$V3
  library_size <- sum(per_chr_nreads)
  per_chr_perc <- per_chr_nreads/library_size
  
  chr_map_df <- as.data.frame(per_chr_perc)
  colnames(chr_map_df) <-strsplit(file, split = "_")[[1]][3]
  chr_map_df_list <- c(chr_map_df_list,chr_map_df)
}

per_chr_reads_df <- t(data.frame(chr_map_df_list))
per_chr_reads_df*100

table(meta$class)
per_chr_perc_order <- c(meta[meta$class=="mMN",][,2:3]$sample,meta[meta$class=="MN",][,2:3]$sample,meta[meta$class=="PN",][,2:3]$sample,meta[meta$class=="Empty",][,2:3]$sample)
per_chr_perc_order_anno <- c(meta[meta$class=="mMN",][,2:3]$class,meta[meta$class=="MN",][,2:3]$class,meta[meta$class=="PN",][,2:3]$class,meta[meta$class=="Empty",][,2:3]$class)




######

##Per_chr_perc  Heatmap



per_chr_perc_order_heatmap_anno <- data.frame(sample = per_chr_perc_order, class = per_chr_perc_order_anno)




pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/per_chr_heatmap_non_norm.pdf", width =10 , height = 10)

heatmap_main_perchr <- pheatmap(as.matrix(per_chr_reads_df[per_chr_perc_order,]*100), color=colorRampPalette(c("white","red","darkred"))(100),cluster_cols = FALSE,cluster_rows=F, show_rownames = T, show_colnames = FALSE, annotation_row = per_chr_perc_order_heatmap_anno)
heatmap_main_perchr
dev.off()

#make same data as boxplots
perc_per_chr_per_class_df <- data.frame()


df_temp <- t(per_chr_reads_df)[,meta[meta$class == "mMN",]$sample]
rownames(df_temp) <- c(paste("chr",1:22, sep=""),"chrX","chrY")
for (i in 1:dim(df_temp)[1]){
  perc_per_chr_per_class_df <- rbind(perc_per_chr_per_class_df,data.frame(cov=as.numeric(df_temp[i,]), chr = rownames(df_temp)[i], class = "mMN"))
  }
df_temp <- t(per_chr_reads_df)[,meta[meta$class == "MN",]$sample]
rownames(df_temp) <- c(paste("chr",1:22, sep=""),"chrX","chrY")
for (i in 1:dim(df_temp)[1]){
  perc_per_chr_per_class_df <- rbind(perc_per_chr_per_class_df,data.frame(cov=as.numeric(df_temp[i,]), chr = rownames(df_temp)[i], class = "MN"))
}
df_temp <- t(per_chr_reads_df)[,meta[meta$class == "PN",]$sample]
rownames(df_temp) <- c(paste("chr",1:22, sep=""),"chrX","chrY")
for (i in 1:dim(df_temp)[1]){
  perc_per_chr_per_class_df <- rbind(perc_per_chr_per_class_df,data.frame(cov=as.numeric(df_temp[i,]), chr = rownames(df_temp)[i], class = "PN"))
}
df_temp <- t(per_chr_reads_df)[,meta[meta$class == "Empty",]$sample]
rownames(df_temp) <- c(paste("chr",1:22, sep=""),"chrX","chrY")
for (i in 1:dim(df_temp)[1]){
  perc_per_chr_per_class_df <- rbind(perc_per_chr_per_class_df,data.frame(cov=as.numeric(df_temp[i,]), chr = rownames(df_temp)[i], class = "Empty"))
}






###

#annotate amplicon
recon_all_amplicon<- makeGRangesFromDataFrame(reconstrucion_bed,keep.extra.columns = TRUE) #reconstruction
ov_vec <- countOverlaps(makeGRangesFromDataFrame(window_calls_counts_normalised[,1:3]),recon_all_amplicon)
window_calls_counts_normalised_benchmark <- window_calls_counts_normalised
window_calls_counts_normalised_benchmark$amplicon <- ov_vec


window_calls_counts_normalised_benchmark$amplicon_name <- "NA"

for (i in 1:dim(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$amplicon !=0 & window_calls_counts_normalised_benchmark$samplename == "A1",])[1]){
  gene_overlsp <- toString(gene_overlap <- data.frame(subsetByOverlaps(recon_all_amplicon, makeGRangesFromDataFrame(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$amplicon !=0 & window_calls_counts_normalised_benchmark$samplename == "A1",][i,1:3])))$circ_id[1])
  window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$amplicon !=0 & window_calls_counts_normalised_benchmark$samplename == "A1",][i,]$amplicon_name  <- gene_overlsp
  print(gene_overlsp)
}

window_calls_counts_normalised_benchmark$amplicon_name <- rep(window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == "A1",]$amplicon_name,48)

window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$amplicon != 0,]


#get top 20 bins for each 

top20bins_df <- data.frame()
for(sample in unique(window_calls_counts_normalised_benchmark$samplename)){
  print(sample)
  df_sample <- window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample,]
  df_sample_ordered <- df_sample[order(df_sample$nReads, decreasing = TRUE),][1:20,]
  
  top20bins_df <- rbind(top20bins_df,data.frame(sample = df_sample_ordered$samplename,nReads = df_sample_ordered$nReads,circlin = df_sample_ordered$amplicon, class = meta[meta$sample == sample,]$class))
  
}
top20bins_df$circlin_class <- "nope"
top20bins_df[top20bins_df$circlin ==1,]$circlin_class <- "circular"
top20bins_df[top20bins_df$circlin ==0,]$circlin_class <- "linear"


pdf("/Users/robinxu/Documents/Projects/ecDNA/Micronuclei/MNseq/region_calling/MNseq_for_reviewers/top20bins.pdf", width =8 , height = 6)

ggplot(data=top20bins_df, aes(x=factor(class, level=c("mMN","MN","PN","Empty")), y=(nReads))) +
  geom_jitter(aes(colour = circlin_class)) +
  ggtitle("") +
  xlab("") + ylab("n Reads")  + theme_classic()  #+ scale_y_log10()
dev.off()




#make readcount table
amplicon_read_count_table <- data.frame()
for(sample in unique(window_calls_counts_normalised_benchmark$samplename)){
  print(sample)
  df_sample <- window_calls_counts_normalised_benchmark[window_calls_counts_normalised_benchmark$samplename == sample,]
  
  amplicon_read_count_table_df_sample <- data.frame(sample=sample, SMC6 = sum(df_sample[df_sample$amplicon_name =="SMC6",]$nReads),MYCN = sum(df_sample[df_sample$amplicon_name =="MYCN",]$nReads),ODC1 = sum(df_sample[df_sample$amplicon_name =="ODC1",]$nReads),CDK4 = sum(df_sample[df_sample$amplicon_name =="CDK4",]$nReads),MDM2 = sum(df_sample[df_sample$amplicon_name =="MDM2",]$nReads) )
  amplicon_read_count_table <- rbind(amplicon_read_count_table,amplicon_read_count_table_df_sample)
  
}

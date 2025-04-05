# ecMN

run_window_calling.R identifies sequences in single micronucleus sequencing data.


# Installation

Requirements:
  - r-base=4.2.2
  - r-ggplot2=3.5.1
  - r-ggpubr=0.6.0
  - r-mass=7.3_60.0.1
  - r-stringr=1.5.1
  - r-tidyr=1.3.1
  - r-fitdistrplus=1.2.1

This code has been tested on aarch64-apple-darwin20 (64-bit) running under: macOS Monterey 12.3.

# How to run this code

To run this code: Rscript run_window_calling.R <meta_table> <reconstruction_bed> <output_directory> 

meta_table: must consist of 3 columns (Path_to_bed_file  Sample  Class). An example can be found under ./test_data/meta_pilot.txt

reconstruction_bed: this currently has no function, but will be used in future releases. Still has to be supplied. An example can be found under ./test_data/TR14_hg38_reconstruction_decoil_with_CDK4_shasta_assembly.bed.txt

output_directory: Directory, where to put QC plots and the window calls.

# Input format
Input files are coverage bed files generated using bedtools coverage. 

Must have following format: chr  start  end  number_of_reads  number_of_bases_covered  binsize  percentage_covered

# Output

Outputs a plain text files in the following format: chr	start	end	nReads	nBases_covered	binsize	percent_covered	samplename	window_class	window_call_final

window_call_final: 1 denotes an positive window



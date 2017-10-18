#!/usr/bin/env Rscript

#-----------------------------------------------------------
#
# The purpose of this program is to count the number of simulated
# insertions in each Snap-Velcro site category that lie within
# exons. The results are stored in a .rda table and also 
# .csv file in
# ../Data/gen-sim_exon_counts.rda 
# and
# ../Data/gen-sim_exon_counts.csv

#
# It assumes the file 
# ../Data/humangenome/Homo_sapiens.GRCh38.89.gff3.exonfilt 
# exists
#
#-----------------------------------------------------------

#--- Load necessary libraries
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
#--- Load simulation output file
load("../Data/gen-sim-out.rda")

#--- Read exonic regions data into a table 'exann'
exann <- read.table("../Data/humangenome/Homo_sapiens.GRCh38.89.gff3.exonfilt")

#--- Allocate memory for a 2D array containing results
exon_counts<-array(0,dim=c(24,9))

#--- Loop through chromosome names
j<-1
for (i in names(Hsapiens)[1:24]) {
	
	# Extract exonic regions for the current chromosome
	chrno<-strsplit(i,"chr")[[1]][2] # Get the current chromosome number (or letter)
        exann_i <- exann[exann$V1 == chrno,] # Extract a subset of the 'exann' table for regions in the current chromosome

	# For each five insertion categories (the four S-V types, plus endonuclease-independent insertions)
	# count the number which lie inside exonic regions provided in 'exann'
	for (k in 1:5) {
		tmp<-inrange(sites_loci[sites_chrm==i & sites_classes==k],exann_i$V2,exann_i$V3)
		exon_counts[j,k] <- length(which(tmp ==TRUE))
	}
	# Count the number of S-V sites of each category for the current chromosome and store in columns 6-9
	for (k in 6:9) {
		exon_counts[j,k]<-length(sites_loci[sites_chrm==i & sites_classes==k-5])
	}
	j<-j+1
}

#--- Store results as a table in a .rda and also in a .csv
save(exon_counts,file="../Data/gen-sim_exon_counts.rda")
write.csv(exon_counts,file="../Data/gen-sim_exon_counts.csv")



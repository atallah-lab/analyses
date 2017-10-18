#!/usr/bin/env Rscript

#-----------------------------------------------------------
#
# The purpose of this program is to count the number of sites
# from each Snap-Velcro category that lie inside exonic
# regions. The results are stored in a .rda table and also 
# .csv file in
# ../Data/exon_counts.rda 
# and
# ../Data/exon_counts.csv
#
# These counts can be used to develop probability models for
# the case of an L1 insertion in an exonic region.
#
# It assumes the file 
# ../Data/humangenome/Homo_sapiens.GRCh38.89.gff3 
# exists, and automatically extracts the exonic regions from 
# it into a table.
#
#-----------------------------------------------------------

#--- Load necessary libraries
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
# Print any warnings during the running of the script
options(warn=1)

#--- Use system (Unix) commands to extract exon locations from a hg38 .gff3 file
cat("\nFiltering hg38 gff3 file for exon locations... \n")
system(paste("grep exon ../Data/humangenome/Homo_sapiens.GRCh38.89.gff3 > tmp")) # Write lines containing "exon" to a tmp file
system(paste("cut -f1,4,5 tmp > ../Data/humangenome/Homo_sapiens.GRCh38.89.gff3.exonfilt")) # Extract columns 1, 4, 5 from tmp file, and write result to ../Data directory
system(paste("rm tmp"))
cat("Reading exon location table...\n")
exann <- read.table("../Data/humangenome/Homo_sapiens.GRCh38.89.gff3.exonfilt") # Read exonic regions into table (chromName	start	end)

exon_counts <- 	array(0,dim=c(24,9)) # Allocate memory for a 2D array containing counts of exonic insertions of each S-V type for each chromosome, and other information

#--- Loop through chromosome names
j<-1
for (i in names(Hsapiens)[1:24]){

	cat("\nProccessing ",i,"...")
	load(paste0("../Data/",i,"map.rda")) # Load map file for current chromosome

	# Data objects containing indices of S-V sites for each category are labeled with the chromosome name.
	# Here we copy the data objects to a set of new names which can be used consistently in the following
        # code.	
        ict<-get(paste0(i,"ict")) 
        icl<-get(paste0(i,"icl"))
        iot<-get(paste0(i,"iot"))
        iol<-get(paste0(i,"iol"))
        insites<-get(paste0(i,"insites"))

	# Extract exonic regions for the current chromosome
	chrno<-strsplit(i,"chr")[[1]][2] # Get the current chromosome number (or letter)
	exann_i <- exann[exann$V1 == chrno,] # Extract a subset of the 'exann' table for regions in the current chromosome

	# Fill in columns 1-4 of the row of the exon_counts array corresponding to current chromosome.
	# These columns contain the number of sites of each S-V category lying within exons 
	# Categories are in this order: closed-tight, closed-loose, open-tight, open-loose
	tmp<-inrange(insites[as.vector(ict[which(!is.na(ict))])],exann_i$V2,exann_i$V3) # Check if any Closed-Tight category sites are within the start-end range of exann_i
	exon_counts[j,1]<-length(which(tmp == TRUE)) 					# Fill an element of the exon_counts table with the count
        tmp<-inrange(insites[as.vector(icl[which(!is.na(icl))])],exann_i$V2,exann_i$V3)
        exon_counts[j,2]<-length(which(tmp == TRUE))
        tmp<-inrange(insites[as.vector(iot[which(!is.na(iot))])],exann_i$V2,exann_i$V3)
        exon_counts[j,3]<-length(which(tmp == TRUE))
        tmp<-inrange(insites[as.vector(iol[which(!is.na(iol))])],exann_i$V2,exann_i$V3)
        exon_counts[j,4]<-length(which(tmp == TRUE))

	# Fill in columns 1-4 of the row of the exon_counts array corresponding to current chromosome.
	# These columns contain the total number of sites of each category
	exon_counts[j,5]<-length(which(!is.na(ict)))
	exon_counts[j,6]<-length(which(!is.na(icl)))
	exon_counts[j,7]<-length(which(!is.na(iot)))
	exon_counts[j,8]<-length(which(!is.na(iol)))

	# The last column contains the exon fraction of the current chromosome
	tmp<-IRanges(exann_i$V2,exann_i$V3)
	exon_counts[j,9] <- sum(width(tmp))/length(Hsapiens[[i]])

	j<-j+1
}
#--- Store results as a table in a .rda and also in a .csv
cat("\n")
save(exon_counts,file="../Data/exon_counts.rda")
write.csv(exon_counts,file="../Data/exon_counts.csv")







#!/usr/bin/env Rscript

library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
for (i in names(Hsapiens)[1:24]){ # load all chromosome map files
        cat(paste0('Loading target site annotation: ',i,'\n'))
        load(paste0('../sim-develop/data/root_maps/',i,'.rda'))
}
load('../sim-develop/data/exann.rda')

gnames <- c()
numGenes <- 0
for (ii in names(Hsapiens)[1:24]){
    map<-get(paste0(ii,'Map')) # Get the site map data for the current chrom
    ict<-    map[[2]]
    icl<-    map[[3]]
    iot<-    map[[4]]
    iol<-    map[[5]]
    insites<-map[[1]]

    tmpann_1 <- exann[exann$chrom==ii,]

    numGenes <- numGenes + length(unique(unlist(tmpann_1$geneSym)))

}
counts <-  array(NA,dim=c(numGenes,4)) # Allocate matrix for counts

l<-1
for (ii in names(Hsapiens)[1:24]){

    ptm <- proc.time()
    cat(paste0('######## Chromosome ',ii,' #########\n'))

    map<-get(paste0(ii,'Map')) # Get the site map data for the current chrom
    ict<-    map[[2]]
    icl<-    map[[3]]
    iot<-    map[[4]]
    iol<-    map[[5]]
    insites<-map[[1]]

    tmpann_1 <- exann[exann$chrom==ii,]

    currGenes <- unique(unlist(tmpann_1$geneSym))
    numCurrGenes <- length(currGenes)
    gnc <- 0
    for (jj in currGenes) {

	gnc <- gnc+1
	cat(paste0('gene ',gnc,'/',numCurrGenes,' for ',ii,':\t',jj,'\n'))

        # count sites in exons of current gene
        tmpann_2 <- tmpann_1[grep(jj,tmpann_1$geneSym),]

        sen        <-inrange(insites[ict[!is.na(ict[,1]),1],1],tmpann_2$start,tmpann_2$end) # Check if any Closed-Tight category sites are within the start-end range of tmpann_2
        antisen    <-inrange(insites[ict[!is.na(ict[,2]),2],2],tmpann_2$start,tmpann_2$end)
        counts[l,1]<-length(which(sen)) + length(which(antisen))                            # Fill an element of the counts table

        sen        <-inrange(insites[icl[!is.na(icl[,1]),1],1],tmpann_2$start,tmpann_2$end)
        antisen    <-inrange(insites[icl[!is.na(icl[,2]),2],2],tmpann_2$start,tmpann_2$end)
        counts[l,2]<-length(which(sen)) + length(which(antisen))

        sen        <-inrange(insites[iot[!is.na(iot[,1]),1],1],tmpann_2$start,tmpann_2$end)
        antisen    <-inrange(insites[iot[!is.na(iot[,2]),2],2],tmpann_2$start,tmpann_2$end)
        counts[l,3]<-length(which(sen)) + length(which(antisen))

        sen        <-inrange(insites[iol[!is.na(iol[,1]),1],1],tmpann_2$start,tmpann_2$end)
        antisen    <-inrange(insites[iol[!is.na(iol[,2]),2],2],tmpann_2$start,tmpann_2$end)
        counts[l,4]<-length(which(sen)) + length(which(antisen))

        gnames <- append(gnames,jj)

        l <- l+1

    }
    print(proc.time() - ptm)
}

save(counts, gnames, file='./gene_exon_site_counts.rda')

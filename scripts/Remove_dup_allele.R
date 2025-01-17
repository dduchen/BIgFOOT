#!/usr/bin/env Rscript

# -- use this to find alleles for the gene of interest shared by the IGenotyper reference allele set
suppressMessages(suppressWarnings(library(DECIPHER)))
args<-commandArgs(TRUE)
alleles=args[1]
graph_alleles=gsub("order.","",alleles)
alleles<-readDNAStringSet(alleles)
alleles_aln<-AlignSeqs(alleles,verbose=F)
alleles_dist<-data.frame(DistanceMatrix(alleles_aln,verbose=F),check.names=F)
dup_alleles <- sapply(alleles_dist, function(column) sum(column == 0) > 1)
# Get the names of such columns
if(length(names(dup_alleles[dup_alleles]))>0){
    print(paste0("Duplicate alleles found - removing... ",names(dup_alleles[dup_alleles])[grep("/",names(dup_alleles[dup_alleles]),invert=T)]))
    remove_alleles<-names(dup_alleles[dup_alleles])[grep("/",names(dup_alleles[dup_alleles]),invert=T)]
    alleles_filt<-alleles[setdiff(names(alleles),remove_alleles)]
} else {
    alleles_filt<-alleles
}
#
writeXStringSet(alleles_filt,file=graph_alleles)
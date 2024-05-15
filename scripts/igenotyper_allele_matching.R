#!/usr/bin/env Rscript

# -- use this to find alleles for the gene of interest shared by the IGenotyper reference allele set
suppressMessages(library(DECIPHER))
args<-commandArgs(TRUE)
igeno_alleles=args[1]
#igeno_alleles="/home/dd392/pi_kleinstein/bigfoot/1kgenomes/crams/igl_samples/NA19761_wg_immunovar_genotyping/familywise_pe_haplotype_inference/NA19761.wg_immunovar.IGHV3-66.igeno_alleles.fasta"
graph_alleles=gsub("igeno_alleles.fasta","igeno_filtering.fasta",igeno_alleles)
igeno_alleles_fa<-readDNAStringSet(igeno_alleles)
graph_alleles_fa<-readDNAStringSet(graph_alleles)
# should be exact matches
graph_alleles_filt<-graph_alleles_fa[grep("IMGT",names(graph_alleles_fa))]
graph_alleles_filt<-graph_alleles_filt[which(graph_alleles_filt %in% igeno_alleles_fa)]
if(length(graph_alleles_filt)!=length(igeno_alleles_fa)){
    print("Mismatch btwn igenotyper alleles and reference alleles - taking closest matches of appropriate gene...")
    gene=gsub(".*wg_immunovar.","",gsub(".igeno_.*$","",igeno_alleles))
    graph_alleles_filt<-graph_alleles_fa[grep(gene,names(graph_alleles_fa))]
    graph_alleles_filt<-graph_alleles_filt[grep("IMGT",names(graph_alleles_filt))]
    matches<-"NA"
    for(i in 1:length(igeno_alleles_fa)){
        igeno_alleles_tmp<-igeno_alleles_fa[i]
        igeno_alleles_tmp_dist<-DistanceMatrix(AlignSeqs(c(igeno_alleles_tmp,graph_alleles_filt),verbose=F),verbose=F)[-1,1]
        matches<-c(matches,names(which(igeno_alleles_tmp_dist==min(igeno_alleles_tmp_dist))))
    }
    matches<-matches[matches!="NA"]
    writeXStringSet(graph_alleles_filt[matches],file=gsub("igeno_filtering.fasta","igeno_filtered.fasta",graph_alleles))
} else {
    print("Matches!")
    writeXStringSet(graph_alleles_filt,file=paste0(gsub("igeno_filtering.fasta","igeno_filtered.fasta",graph_alleles)))
}
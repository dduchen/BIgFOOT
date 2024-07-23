#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
# ensure relevant packages can be loaded - either uncomment + add specific libpaths or download
#.libPaths( c( .libPaths(), "/vast/palmer/home.mccleary/dd392/R/x86_64-pc-linux-gnu-library/4.2", "/vast/palmer/apps/avx2/software/R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0") )
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
list.of.packages <- c("DECIPHER", "data.table","dplyr","Biostrings")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) suppressMessages(suppressWarnings(suppressPackageStartupMessages(BiocManager::install(new.packages))))
#
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(DECIPHER))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(data.table))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(dplyr))))

contig_file=args[1]
sample<-gsub("\\..*","",contig_file)
haps_file<-gsub("contigs","haps.final",contig_file)
graph_file<-gsub("contigs.fasta","genome_graph.gfa",contig_file)
#
contigs<-readDNAStringSet(contig_file)
# remove contig with 0 abundance!
contigs<-contigs[grep("ab=0.0000",names(contigs),invert=T)]
# limit to alleles of interest - important if VG-flow used on haps + alleles:
# second round of inference - path names have ":path" embedded in them
contigs<-contigs[grep("IMGT|IGv2|OGRDB|HLA|KIR|:path",names(contigs))]
haps<-readDNAStringSet(haps_file)
for(i in 1:length(haps)){
    contig_match<-contigs[which(contigs==haps[i])]
    if(length(contig_match)==0){
        aln_dist<-DistanceMatrix(c(haps,contigs))
        aln_seqs<-AlignSeqs(c(haps,contigs))
        aln_dist<-as.matrix(stringDist(c(aln_seqs[names(haps[i])],aln_seqs[names(contigs)]),method="hamming"))
        aln_dist<-aln_dist[names(contigs),names(haps[i])]
        if(length(aln_dist[which(aln_dist==min(aln_dist))])>1){
            tmp_name_haps<-paste0(gsub("^.*#1#","",gsub("#0.*$","",names(aln_dist[which(aln_dist==min(aln_dist))]))),":",names(haps)[i],"_HammingDistBased")
            tmp_name_haps[2:length(tmp_name_haps)]<-gsub(":path.*$","",tmp_name_haps[2:length(tmp_name_haps)])
            tmp_name_haps<-paste(tmp_name_haps[2:length(tmp_name_haps)],tmp_name_haps[1],sep=" or ")
            names(haps)[i]<-tmp_name_haps
        } else {
            names(haps)[i]<-paste0(gsub("^.*#1#","",gsub("#0.*$","",names(aln_dist[which(aln_dist==min(aln_dist))]))),":",names(haps)[i],"_HammingDistBased")
        }
        #Biostrings::stringDist --> minimal levenshtein distance or align first --> minimal hamming distance. With this method could potentially use the either the haplotypes/IMGT alleles here as queries for inference
    } else {
    names(haps)[i]<-paste0(gsub("^.*#1#","",gsub("#0.*$","",names(contig_match))),":",names(haps)[i])
    }
}
writeXStringSet(haps,file=gsub("final","final.annot",haps_file))
#
gfa<-fread(graph_file,fill=T,header=F)
plines<-gfa[grep("^P",gfa$V1),]
haps_path_ids<-cbind(names(haps),gsub("^.*path","",gsub(" .*$","",names(haps))))
for(i in 1:nrow(plines)){
   graph_path<-plines$V2[i]
   graph_new_id<-haps_path_ids[haps_path_ids[,2]==graph_path,][1]
   gfa[gfa$V1=="P" & gfa$V2==graph_path,]$V2<-graph_new_id
}
outfile<-gsub("graph.gfa","graph.annot.gfa",graph_file)
fwrite(gfa,file=outfile,sep="\t",quote=F,col.names=F,row.names=F)
system(paste0("tr -s '\t' < ",outfile," > tmp.gfa && mv tmp.gfa ",outfile))
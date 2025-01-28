#!/usr/bin/env Rscript
workingdir=getwd()
# -- use this to find nodes not shared with alleles we care about  -- then vg find to extract reads that preferentially map to those non-target alleles and remove them from consideration
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
args<-commandArgs(TRUE)
graph_paths_file=args[1]
#graph_paths_file="/home/dduchen/Documents/bigfoot/1kgenomes/HG00438_wg_immunovar_genotyping/familywise_pe_haplotype_inference/HG00438.wg_immunovar.IGHV4-4.haplotypes.pathnodes"
graph_paths=fread(graph_paths_file,header=F)
dat<-graph_paths[,c(2,3)]
dat<-dat[unique(grep("^IM|^IGv|^OGR|^KIR|^HLA",dat$V2))]
dat$V3<-gsub("\\+","",dat$V3)
dat$V3<-gsub("\\-","",dat$V3)
#--
pathref<-as.list(unique(dat$V2))
names(pathref)<-pathref
#-- for each reference (names), entries are all nodes
for(i in 1:length(unique(dat$V2))){
#  nodes<-as.numeric(strsplit(dat[i,]$V3,split=",")[[1]])
  nodes<-as.numeric(unique(strsplit(paste0(dat[dat$V2==unique(dat$V2)[i],]$V3,collapse=","),split=",")[[1]]))
  pathref[[unique(dat$V2)[i]]]<-nodes
}
#
gene_id=gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file));
target_nodes<-unique(unname(unlist(pathref[grep(gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file)),names(pathref))])))
off_target_nodes<-unique(unname(unlist(pathref[grep(gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file)),names(pathref),invert = T)])))
nodes_for_filtering<-setdiff(off_target_nodes,target_nodes)
#
if(length(target_nodes)==0){
  nodes_for_filtering<-NULL
  gene_prefix<-gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file))
  #check for orphon alleles -- if orphon involvement, this is on another component of the graph - all nodes should be used for filtering
  if(length(grep("/",names(pathref)))>0){
    orphon_nodes<-unique(unname(unlist(pathref[grep("/",names(pathref))])))
    nodes_for_filtering<-orphon_nodes
    if(length(grep("/",names(pathref),invert=T))==0){
      print("Orphon-only ASC cluster")
      nodes_for_filtering<-NULL
    }
  }
}
###############################################################################################################
suppressMessages(suppressWarnings(fwrite(data.frame(nodes_for_filtering),file=gsub("pathnodes","filteringnodes",graph_paths_file),quote=F,sep="\t",col.names=F,row.names=F)))
###############################################################################################################
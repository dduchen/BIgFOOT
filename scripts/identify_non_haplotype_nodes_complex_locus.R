#!/usr/bin/env Rscript
workingdir=getwd()
# -- use this to find nodes not shared with alleles we care about  -- then vg find to extract reads that preferentially map to those non-target alleles and remove them from consideration
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
args<-commandArgs(TRUE)
hap_exact_file=args[1]
hap_offtarget_file=gsub("exact.txt","offtarget.txt",hap_exact_file)
graph_paths_file=args[2]
#
graph_paths=fread(graph_paths_file,header=F)
dat<-graph_paths[,c(2,3)]
#dat<-dat[unique(grep("^IM|^IGv|^OGR|^KIR|^HLA",dat$V2,invert=T))]
dat$V3<-gsub("\\+","",dat$V3)
dat$V3<-gsub("\\-","",dat$V3)
#--
pathref<-as.list(unique(dat$V2))
names(pathref)<-pathref
#-- for each reference (names), entries are all nodes
for(i in 1:length(unique(dat$V2))){
  nodes<-as.numeric(strsplit(dat[i,]$V3,split=",")[[1]])
  pathref[[i]]<-nodes
}
# gene-specific nodes
target_gene_nodes<-unique(unname(unlist(pathref[grep(gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file)),names(pathref))])))
off_target_gene_nodes<-pathref[grep("^IM|^IGv|^OGR|^KIR|^HLA",names(pathref))]
off_target_gene_nodes<-unique(unname(unlist(off_target_gene_nodes[grep(gsub("^.*wg_immunovar.","",gsub(".haplotypes.*$","",graph_paths_file)),names(off_target_gene_nodes),invert=T)])))
#
nodes_for_filtering<-setdiff(off_target_gene_nodes,target_gene_nodes)
# -- haps containing these nodes?
target_haps=fread(hap_exact_file,header=F)
offtarget_haps=fread(hap_offtarget_file,header=F)
#
target_haps$offtarget_hits<-0
for(hap in target_haps$V1){
    tmp_path<-pathref[grep(hap,names(pathref))]
    prop_offtarget<-prop.table(table(nodes_for_filtering %in% unname(unlist(tmp_path))))["TRUE"]
    if(is.na(prop_offtarget)){
        prop_offtarget<-0
    }
    target_haps[target_haps$V1==hap,]$offtarget_hits<-as.numeric(prop_offtarget)
}
offtarget_haps$offtarget_hits<-0
for(hap in offtarget_haps$V1){
    tmp_path<-pathref[grep(hap,names(pathref))]
    prop_offtarget<-prop.table(table(nodes_for_filtering %in% unname(unlist(tmp_path))))["TRUE"]
    if(is.na(prop_offtarget)){
        prop_offtarget<-0
    }
    offtarget_haps[offtarget_haps$V1==hap,]$offtarget_hits<-as.numeric(prop_offtarget)
}
# maybe try find cluster - use find_threshold
target_haps<-target_haps[target_haps$offtarget_hits<=0.5,]
offtarget_haps<-offtarget_haps[offtarget_haps$offtarget_hits>=0.8,]
#
#offtarget_haps<-setdiff(offtarget_haps$V1,target_haps$V1)
target_hap_nodes<-unique(unname(unlist(pathref[grep(paste(target_haps$V1,collapse="|"),names(pathref))])))
off_target_hap_nodes<-unique(unname(unlist(pathref[grep(paste(offtarget_haps$V1,collapse="|"),names(pathref))])))
#
if(length(target_hap_nodes)>0 & length(off_target_hap_nodes)>0){
    nodes_for_filtering_haps<-setdiff(off_target_hap_nodes,target_hap_nodes)
    nodes_for_filtering<-unique(c(nodes_for_filtering,nodes_for_filtering_haps))
}
###############################################################################################################
fwrite(data.frame(nodes_for_filtering),file=gsub("pathnodes","filteringnodes",graph_paths_file),quote=F,sep="\t",col.names=F,row.names=F)
###############################################################################################################
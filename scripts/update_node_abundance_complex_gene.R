#!/usr/bin/env Rscript
workingdir=getwd()
library(data.table);library(stringr)
args<-commandArgs(TRUE)
#arg_1=${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt
#arg_2=${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
#arg_3=${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt
# -- use this to find nodes not shared with alleles we care about  -- then vg find to extract reads that preferentially map to those non-target alleles and remove them from consideration

graph_paths_file=args[2]
#graph_paths_file="/home/dduchen/Documents/bigfoot/1kgenomes/HG00438_wg_immunovar_genotyping/familywise_pe_haplotype_inference/HG00438.wg_immunovar.IGHV4-4.vgflow.final.pathnodes"
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
  nodes<-as.numeric(strsplit(dat[i,]$V3,split=",")[[1]])
  pathref[[i]]<-nodes
}
path_summary_file=args[3]
#path_summary_file="/home/dduchen/Documents/bigfoot/1kgenomes/HG00438_wg_immunovar_genotyping/familywise_pe_haplotype_inference/HG00438.wg_immunovar.IGHV4-4.component.summary.txt"
path_summary=fread(path_summary_file,header=F)
component_indices <- which(str_detect(path_summary$V1, "^component_"))
components <- lapply(1:length(component_indices), function(i) {
  if (i < length(component_indices)) {
    path_summary$V1[(component_indices[i]+1):(component_indices[i+1]-1)]
  } else {
    path_summary$V1[(component_indices[i]+1):nrow(path_summary)]
  }
})
names(components)<-path_summary$V1[grep("component",path_summary$V1)]
#







node_abundance_file=args[1]
#node_abundance_file="/home/dduchen/Documents/bigfoot/1kgenomes/HG00438_wg_immunovar_genotyping/familywise_pe_haplotype_inference/HG00438.wg_immunovar.IGHV4-4.vgflow.node_abundance.txt"
node_abundance=fread(node_abundance_file,header=F)
node_abundance$node_id<-as.numeric(gsub(":.*$","",node_abundance$V1));
node_abundance$node_abundance<-as.numeric(gsub(".*:","",node_abundance$V1));

node_abundance<-data.frame(node_abundance)
node_abundance$alt_abundance<-node_abundance$node_abundance
for(i in names(components[grep("IMGT|OGRDB|IGv2",components)])){
  target_genes<-unname(unlist(components[i]))[grep("IMGT|OGRDB|IGv2",unname(unlist(components[i])))]
  other_genes<-length(unname(unlist(components[i]))[-grep("IMGT|OGRDB|IGv2",unname(unlist(components[i])))])
  comp_nodes<-unique(as.numeric(unname(unlist(pathref[target_genes]))))
  node_abundance[node_abundance$node_id %in% c(comp_nodes),]$alt_abundance<-node_abundance[node_abundance$node_id %in% c(comp_nodes),]$alt_abundance/(other_genes/(other_genes+1))
}
node_abundance$node_abundance_update<-paste0(node_abundance$node_id,":",node_abundance$alt_abundance)
###############################################################################################################
fwrite(file=node_abundance_file,data.frame(node_abundance$node_abundance_update),quote=F,col.names=F,sep="\t")
###############################################################################################################
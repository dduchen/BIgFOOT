#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
covdat_file=args[1]
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
covdat<-data.frame(fread(covdat_file))
colnames(covdat)<-gsub("node.","",colnames(covdat))
gfa_header<-fread(gsub(".coverage",".gfa.header",covdat_file),sep="\t",fill=T,header=F)
gfa_s<-fread(gsub(".coverage",".gfa.slines",covdat_file),sep="\t",fill=T,header=F)
gfa_l<-fread(gsub(".coverage",".gfa.llines",covdat_file),sep="\t",fill=T,header=F)
gfa_p<-fread(gsub(".coverage",".gfa.plines",covdat_file),sep="\t",fill=T,header=F)
gfa<-rbind(gfa_header,gfa_s,gfa_l,gfa_p,use.names=F,fill=T)
gfa_node_order<-gfa[gfa$V1=="S",]$V2
gfa[gfa$V1=="S",4]<-paste0("DP:f:",as.vector(covdat[gfa_node_order]))
header<-gfa[1,]
sline<-gfa[gfa$V1=="S",]
lline<-gfa[gfa$V1=="L",]
pline<-gfa[gfa$V1=="P",]
# alignment-implied nodes?
gene_interest<-gsub(".*\\.","",gsub(".genome_graph_ref.*","",covdat_file))
sample_identifier<-gsub(".*\\/","",gsub(paste0(".wg_immunovar.",gene_interest,".*"),"",covdat_file))
gene_nodes<-pline[grep(paste0(gene_interest),pline$V2),]$V3
ref_path_nodes<-pline[grep(paste0("grch38|chm13"),pline$V2),]$V3
aln_nodes<-pline[grep(paste0(gene_interest,"|grch38|chm13"),pline$V2,invert=T),]$V3
#
gene_nodes_uniq<-unique(unlist(strsplit(paste(gsub("\\+|\\-","",gene_nodes),collapse=","),split=",")))
ref_path_nodes_uniq<-unique(unlist(strsplit(paste(gsub("\\+|\\-","",ref_path_nodes),collapse=","),split=",")))
aln_nodes_uniq<-unique(unlist(strsplit(paste(gsub("\\+|\\-","",aln_nodes),collapse=","),split=",")))
# nodes with alignments vs. gene-related nodes
novel_nodes<-setdiff(aln_nodes_uniq,gene_nodes_uniq)
#exclude 'novel' node if its a tip (only 1 edge)
novel_node_edges<-rbind(lline[lline$V4 %in% novel_nodes,c(2,4)],lline[lline$V2 %in% novel_nodes,c(2,4)])
novel_nodes<-intersect(novel_nodes,unlist(novel_node_edges)[duplicated(unlist(novel_node_edges))])
# should also care if novel node has an edge to gene of interest - allow for degree=1 away from gene-associated node
gene_node_edges<-rbind(lline[lline$V4 %in% gene_nodes_uniq,c(2,4)],lline[lline$V2 %in% gene_nodes_uniq,c(2,4)])
#
# novel_nodes<-intersect(novel_nodes,unlist(gene_node_edges)) # also want edges in/out - otherwise will report node attached to last gene-related node
novel_nodes<-intersect(novel_nodes,unlist(gene_node_edges)[duplicated(unlist(gene_node_edges))])
#
# -- could just get alternative node and extract a local subgraph 'vg mod -x 1'
if(length(novel_nodes)>0){
    names(novel_nodes)<-paste0(sline[match(as.character(novel_nodes),sline$V2),]$V3,"_",sline[match(as.character(novel_nodes),sline$V2),]$V4)
    for(i in seq_along(novel_nodes)){
        newpline<-gfa[gfa$V1=="P",][1,]
        node_id<-gsub(".*_","",names(novel_nodes)[i])
        newpline[,2]<-paste0(sample_identifier,"#1#",gene_interest,"_variant_",node_id,"#",i)
        newpline[,3]<-paste0(novel_nodes[i],"+")
        newpline[,4]<-"*"
        pline<-rbind(pline,newpline)
        read_support<-pline[grep(gsub("\\+","",newpline$V3),pline$V3),]
        ref_convergence<-read_support[grep("grch|chm",read_support$V2),]
        read_support<-read_support[grep(paste0(gene_interest,"|grch|chm"),read_support$V2,invert=T),]
        # get local sequence + alt sequence
        node_connects<-lline[lline$V4==as.character(novel_nodes[i]) | lline$V2==as.character(novel_nodes[i]),]
        node_connects<-distinct(node_connects)
        node_other<-setdiff(unique(c(node_connects$V2,node_connects$V4)),as.character(novel_nodes[i]))
        if(length(node_other)==2){
            node_from<-node_other[1]
            node_to<-node_other[2]
            seq_from<-sline[sline$V2==as.character(node_other[1]),]$V3
            seq_to<-sline[sline$V2==as.character(node_other[2]),]$V3
            var_seq=paste0(seq_from,gsub("_.*","",names(novel_nodes[i])),seq_to)
        } else if(length(node_other[node_other %in% gene_nodes_uniq])==2){
            # limiting to nodes surrounded by gene of interest
            node_other<-node_other[node_other %in% gene_nodes_uniq]
            node_from<-node_other[1]
            node_to<-node_other[2]
            seq_from<-sline[sline$V2==as.character(node_other[1]),]$V3
            seq_to<-sline[sline$V2==as.character(node_other[2]),]$V3
            var_seq=paste0(seq_from,gsub("_.*","",names(novel_nodes[i])),seq_to)
        } else {
            print(paste0("Tricky local sequence for variant @ node:",novel_nodes[i],":",names(novel_nodes)[i]," - look at graph"))
        }
#        nonvar_node<-lline[lline$V4 %in% node_other | lline$V2 %in% node_other,]
#        nonvar_node<-nonvar_node[grep(as.character(novel_nodes[i]),nonvar_node$V2,invert=T),]
#        nonvar_node<-nonvar_node[grep(as.character(novel_nodes[i]),nonvar_node$V4,invert=T),]
#        nonvar_node<-intersect(unique(c(nonvar_node[nonvar_node$V2 %in% node_other,]$V4)),unique(c(nonvar_node[nonvar_node$V4 %in% node_other,]$V2)))
#        nonvar_seq<-""
#        for(alt_path in 1:length(nonvar_node)){
#            nonvar_node_seq<-sline[sline$V2==as.character(nonvar_node[alt_path]),]$V3
#            node_from<-lline[lline$V4==as.character(nonvar_node[alt_path]) | lline$V2==as.character(nonvar_node[alt_path]),]
#            node_from<-setdiff(c(node_from[node_from$V3=="-",]$V4,node_from[node_from$V3=="+",]$V2),as.character(nonvar_node[alt_path]))
#            seq_from<-sline[sline$V2==as.character(node_from),]$V3
#            node_to<-lline[lline$V4==as.character(nonvar_node[alt_path]) | lline$V2==as.character(nonvar_node[alt_path]),]
#            node_to<-setdiff(c(node_to[node_to$V3=="+",]$V4,node_to[node_to$V3=="-",]$V2),as.character(nonvar_node[alt_path]))
#            seq_to<-sline[sline$V2==as.character(node_to),]$V3
#            nonvar_seq=c(nonvar_seq,paste0(seq_from,nonvar_node_seq,seq_to))
#        }
 #       nonvar_seq<-nonvar_seq[-1]
 #       nonvar_seq<-paste0(nonvar_seq,collapse=":")
        if(nrow(ref_convergence>0)){
            print(paste0(newpline$V2,",",nrow(read_support),"_reads,",nchar(gsub("_.*","",names(novel_nodes[i]))),"_nt,",paste0("ref_overlaps:",paste0(ref_convergence$V2,collapse="_")),",",var_seq,",",paste0(gsub(".*\\/","",gsub(".coverage",".depth.gfa",covdat_file)))))
        } else {
            print(paste0(newpline$V2,",",nrow(read_support),"_reads,",nchar(novel_nodes[i]),"_nt,","novel_nonref,",var_seq,",",paste0(gsub(".*\\/","",gsub(".coverage",".depth.gfa",covdat_file)))))
        }
        pline<-pline[!duplicated(pline$V2),]
#        print(paste0(gene_interest,": Alignment-implied sample specific variation implies ",length(novel_nodes)," variants totaling ",sum(nchar(names(novel_nodes)))," nt deviation from known alleles"))
    }
}
gfa<-rbind(header,sline,lline,pline)
fwrite(gfa,file=paste0(gsub(".coverage",".depth.gfa",covdat_file)),col.names=F,quote=F,sep="\t")
# prune low-confidence variants:
paths_to_prune<-gfa[grep("DP:f:1#",gfa$V2),] 
nodes_to_prune<-gsub("\\+|\\-","",paths_to_prune$V3)
for(node in nodes_to_prune){
    # remove S line = the node
    # remove both/all L lines containing the node
    # edit read path associated with the node, simply excise it out -- edit read ID to indicate this
    # remove P line containing the node-assocaited variant
}





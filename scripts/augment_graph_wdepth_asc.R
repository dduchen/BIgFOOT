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
# account for asc_inference approach:
gene_interest_fullpaths<-unique(pline[[2]][grep("TR|IG",pline[[2]])])
gene_interest<-unique(gsub("\\*.*","",pline[[2]][grep("TR|IG",pline[[2]])]))
gene_interest_prefix<-unique(gsub(".-.*","",gene_interest))
if(length(grep("/",gene_interest_prefix))>0){
    print("Orphon involvement")
    gene_interest_prefix<-gsub(".$","",gene_interest_prefix[1])
}
if(all(gene_interest_prefix %in% c("IGKV","IGKV1","IGKV2","IGKV3","IGKV4","IGKV5","IGKV6","IGKV7"))){
    print("IGKV + IGKV2 ASC cluster")
    gene_interest_prefix<-gsub("2$","",gene_interest_prefix[1])
}
#gene_interest<-gsub(".*\\.","",gsub(".genome_graph_ref.*","",covdat_file))
sample_identifier<-gsub(".*\\/","",gsub(".wg_immunovar.*","",covdat_file))
gene_nodes<-pline[pline$V2 %in% gene_interest_fullpaths,]$V3
ref_path_nodes<-pline[grep(paste0("grch38|chm13"),pline$V2),]$V3
aln_nodes<-pline[grep("IG|TR|grch38|chm13",pline$V2,invert=T),]$V3
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
# -- deletion: edges not traversed by alleles - get reads that support the traveral to get read support and embed this as a novel variant - perhaps by including a node with no sequence?
aln_nodes_candidates<-vector()
for(iter in seq_along(gene_nodes)){
    aln_nodes_candidates_tmp<-aln_nodes[grep(gsub("\\+","\\\\+",gsub(",","|",gene_nodes[iter])),aln_nodes)]
    aln_nodes_candidates<-unique(c(aln_nodes_candidates,aln_nodes_candidates_tmp))
}
aln_nodes_candidates<-aln_nodes[aln_nodes %in% aln_nodes_candidates]
#
# look for alns which skip any nodes - require adequate support
if(length(gene_nodes_uniq)>1){
    gene_node_pairs<-combn(gene_nodes_uniq,2,simplify=F)
    } else {
    gene_node_pairs<-gene_nodes_uniq
}
names(gene_node_pairs)<-0
for(i in seq_along(gene_node_pairs)){
    combo<-paste0(gene_node_pairs[[i]],collapse=",")
    combo_rev<-paste0(rev(gene_node_pairs[[i]]),collapse=",")
    if(length(grep(paste0(combo,"|",combo_rev),paste0(gsub("\\+|-","",gene_nodes))))>0){
        gene_node_pairs[[i]]<-NA
    } else {
        # get number of reads supporting the combo - putative deletion - exact matches to node ids
        names(gene_node_pairs)[i]<-length(grep(paste0("^",combo,"$|^",combo,",|^",combo_rev,"$|^",combo_rev,","),paste0(gsub("\\+|-","",aln_nodes))))
    }
}
gene_node_pairs<-gene_node_pairs[!is.na(gene_node_pairs)]
gene_node_pairs_candidates<-gene_node_pairs[as.numeric(names(gene_node_pairs))>2]
#
# -- add the deletion path (2 nodes) to the graph
for(i in seq_along(gene_node_pairs_candidates)){
        newpline<-gfa[gfa$V1=="P",][1,];
        deletion_path<-gene_node_pairs_candidates[i];
        deletion_path_id<-paste0(sample_identifier,"#1#",gene_interest_prefix,"_deletion_DP:f:",names(gene_node_pairs_candidates[i]),"#",i);
        newpline[,2]<-deletion_path_id;
        node1<-deletion_path[[1]][1];
        node2<-deletion_path[[1]][2];
        node_connects<-lline[lline$V4==as.character(node1) & lline$V2==as.character(node2) | lline$V4==as.character(node2) & lline$V2==as.character(node1),];
        node_connects<-distinct(node_connects);
        if(nrow(node_connects)==1){
        newpline[,3]<-paste0(node_connects$V2,node_connects$V3,",",node_connects$V4,node_connects$V5)
        } else {
            print(paste0("Tricky local sequence for deletion @ node:",node1,"-",node2," - look at graph"))
        }
        pline<-rbind(pline,newpline)
}
# -- could just get alternative node and extract a local subgraph 'vg mod -x 1'
if(length(novel_nodes)>0){
    names(novel_nodes)<-paste0(sline[match(as.character(novel_nodes),sline$V2),]$V3,"_",sline[match(as.character(novel_nodes),sline$V2),]$V4)
    for(i in seq_along(novel_nodes)){
        newpline<-gfa[gfa$V1=="P",][1,]
        node_id<-gsub(".*_","",names(novel_nodes)[i])
        newpline[,2]<-paste0(sample_identifier,"#1#",gene_interest_prefix,"_variant_",node_id,"#",i)
        # we care about direction - try match with adjacent nodes if they're in allele sequences:
        node_support<-aln_nodes[grep(paste0(",",novel_nodes[[i]],"\\+|,",novel_nodes[[i]],"\\-"),aln_nodes)]
        suffix_pos<-length(grep(paste0(",",novel_nodes[[i]],"\\+"),node_support))
        local_read_pos<-node_support[grep(paste0(",",novel_nodes[[i]],"\\+"),node_support)]
        local_read_pos<-strsplit(local_read_pos,",")
        # checking allele paths for direction - if not then use the number of reads supporting the node as a proxy for the direction
        for(j in seq_along(local_read_pos)){
            local_read_pos[[j]]<-local_read_pos[[j]][c((grep(novel_nodes[[i]],local_read_pos[[j]])-1),(grep(novel_nodes[[i]],local_read_pos[[j]]))+1)]
            if(all(local_read_pos[[j]] %in% unlist(strsplit(gene_nodes,",")))){
                names(local_read_pos)[j]<-"1"
            }
        }
        suffix_neg<-length(grep(paste0(",",novel_nodes[[i]],"\\-"),node_support))
        local_read_neg<-node_support[grep(paste0(",",novel_nodes[[i]],"\\-"),node_support)]
        local_read_neg<-strsplit(local_read_neg,",")
        for(j in seq_along(local_read_neg)){
            local_read_neg[[j]]<-local_read_neg[[j]][c((grep(novel_nodes[[i]],local_read_neg[[j]])-1),(grep(novel_nodes[[i]],local_read_neg[[j]]))+1)]
            if(all(local_read_neg[[j]] %in% unlist(strsplit(gene_nodes,",")))){
                names(local_read_neg)[j]<-"1"
            }
        }
        if(any(as.numeric(c(names(local_read_neg),names(local_read_pos)))>0)){
            suffix_neg<-sum(as.numeric(names(local_read_neg)))
            suffix_pos<-sum(as.numeric(names(local_read_pos)))
        }
        if(suffix_pos>=suffix_neg){
            newpline[,3]<-paste0(novel_nodes[i],"+")
        } else {
            newpline[,3]<-paste0(novel_nodes[i],"-")
        }
        newpline[,4]<-"*"
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
        # append local region to graph
        sign_in_pos<-length(grep(paste0(node_to,"\\+,",newpline[,3]),node_support))
        sign_in_neg<-length(grep(paste0(node_to,"\\-,",newpline[,3]),node_support))
        if(sign_in_pos>=sign_in_neg){
            sign_in="+"
        } else {
            sign_in="-"
        }
        sign_out_pos<-length(grep(paste0(",",gsub("\\+","\\\\+",newpline[,3]),",",node_to,"\\+"),node_support))
        sign_out_neg<-length(grep(paste0(",",gsub("\\+","\\\\+",newpline[,3]),",",node_to,"-"),node_support))
        if(sign_out_pos>=sign_out_neg){
            sign_out="+"
        } else {
            sign_out="-"
        }
        newpline$V3<-paste0(node_from,sign_in,",",newpline$V3,",",node_to,sign_out)
        pline<-rbind(pline,newpline)
        read_support<-pline[grep(gsub("\\+|\\-","",newpline$V3),gsub("\\+|\\-","",pline$V3)),]
        ref_convergence<-read_support[grep("grch|chm",read_support$V2),]
        read_support<-read_support[grep(paste0(gene_interest_prefix,"|grch|chm"),read_support$V2,invert=T),]
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
# prune low-confidence variants:
gfa<-rbind(header,sline,lline,pline)
#
paths_to_prune<-gfa[grep("DP:f:0\\.#|DP:f:0#|DP:f:1#|DP:f:2#|DP:f:1\\.#|DP:f:2\\.#",gfa$V2),] 
if(nrow(paths_to_prune)>0){
    print(paste0("Pruning ",nrow(paths_to_prune)," low-confidence variants"))
    nodes_to_prune<-gsub("\\+|\\-","",paths_to_prune$V3)
    # remove single-read supported variants + prep the gfa for vg deconstruct - want variants and alleles to be listed in vcf file
    sline_orig<-sline
    pline_orig<-pline
    lline_orig<-lline
    for(node in nodes_to_prune){
        variant_node<-strsplit(node,split=",")[[1]][2]
        edge_in<-strsplit(node,split=",")[[1]][1:2]
        edge_out<-strsplit(node,split=",")[[1]][2:3]
        # remove S line = the node
        sline_filt<-sline[grep(paste0("^",variant_node,"$"),sline$V2,invert=T),]
        lline_filt<-lline[-intersect(which(lline$V2 %in% edge_in),which(lline$V4 %in% edge_in)),]
        lline_filt<-lline_filt[-intersect(which(lline_filt$V2 %in% edge_out),which(lline_filt$V4 %in% edge_out)),]
        # remove both/all L lines containing the node
        # edit read path associated with the node, simply excise it out -- edit read ID to indicate this
        paths_traversing_variant_path<-pline[grep(node,gsub("\\+|\\-","",pline$V3))]$V2
        if(length(grep("grch|chm|^IG|^TR",paths_traversing_variant_path))>0){
            print("Important path also traverses the variant chunk - removing read path ID only - not topology")
            specific_path_to_remove<-paths_traversing_variant_path[grep("grch|chm|^IG|^TR",paths_traversing_variant_path,invert=T)]
            if(length(specific_path_to_remove)>0){
                pline_filt<-pline[-which(pline$V2 %in% specific_path_to_remove),]
                pline<-pline_filt
            }
            #sline<-sline_filt
            #lline<-lline_filt
        } else {
            pline_filt<-pline[grep(node,gsub("\\+|\\-","",pline$V3),invert=T),]
            pline<-pline_filt
        }
        # remove paths containing filtered-out variant node provided they are just read paths
        other_var_paths<-pline_filt[grep(variant_node,gsub("\\+|\\-","",pline_filt$V3),invert=F),]
        if(all(nrow(other_var_paths)>0 & length(grep("grch|chm|^IG|^TR",other_var_paths$V2))>0)){
            print("Important path also traverses the variant node - removing read path ID only - not topology")
            specific_path_to_remove<-other_var_paths[grep("grch|chm|^IG|^TR",other_var_paths$V2,invert=T)]
            if(length(specific_path_to_remove$V2)>0){
                pline_filt<-pline[grep(specific_path_to_remove$V2,pline$V2,invert=T),]
                pline<-pline_filt
            }
            #sline<-sline_filt
            #lline<-lline_filt
            #next
        } else {
        # remove other P lines containing the putatively spurious node
            if(nrow(other_var_paths)>0){
                pline_filt<-pline_filt[-which(pline_filt$V2 %in% other_var_paths$V2),]
            }
            sline<-sline_filt
            lline<-lline_filt
            pline<-pline_filt
        }
    }
}
#
gfa<-rbind(header,sline,lline,pline)
fwrite(gfa,file=paste0(gsub(".coverage",".depth.gfa",covdat_file)),col.names=F,quote=F,sep="\t")
#
pline$V2<-gsub(":path.","#1#",pline$V2);pline$V2<-gsub("#1#_","#1#",pline$V2)
gfa_filt<-rbind(gfa_header,sline,lline,pline,fill=T)
fwrite(gfa_filt,file=paste0(gsub(".coverage",".depth.filt.gfa",covdat_file)),col.names=F,quote=F,sep="\t")
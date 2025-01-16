#!/usr/bin/env Rscript
workingdir=getwd()
# -- use this to find nodes not shared with alleles we care about  -- then vg find to extract reads that preferentially map to those non-target alleles and remove them from consideration
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
args<-commandArgs(TRUE)
path_file=args[1]
#path_file="NA18874.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines" # many variants - multi-haps
#path_file="NA18559.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines" # 2 variants
#path_file="NA19916.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines" # 2 variants, 1 below thresh
#path_file="HG00451.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG00119.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG02182.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines" #< -- mismatch found - figure this out
#path_file="HG02886.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG02620.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG02622.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG02886.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="HG03492.wg_immunovar.IGHV3-9.genome_graph_ref.augmented.depth.plines"
#path_file="NA20787.wg_immunovar.IGHV1-69.genome_graph_ref.augmented.depth.plines"
#path_file="HG02717.wg_immunovar.IGHV1-69.genome_graph_ref.augmented.depth.plines"
graph_paths=fread(path_file,header=F)
dat<-graph_paths[,c(2,3)]
#dat<-dat[unique(grep("^IM|^IGv|^OGR|^KIR|^HLA",dat$V2,invert=T))]
dat$V3<-gsub("\\+","",dat$V3)
dat$V3<-gsub("\\-","",dat$V3)
# paired end reads - ensure same read ID
dat$V2<-gsub("-1$","",dat$V2)
dat$V2<-gsub("-2$","",dat$V2)
dat$V2<-gsub("/1$","",dat$V2)
dat$V2<-gsub("/2$","",dat$V2)
#--
pathref<-as.list(unique(dat$V2))
names(pathref)<-pathref
#-- for each reference (names), entries are all nodes
for(i in 1:length(unique(dat$V2))){
  nodes<-as.numeric((strsplit(paste0(dat[dat$V2==unique(dat$V2)[i],]$V3,collapse=","),split=",")[[1]]))
  pathref[[unique(dat$V2)[i]]]<-nodes
}
#
variant_paths<-pathref[grep("variant",names(pathref))]
variant_paths_filt<-variant_paths[grep(":f:1#|:f:2#",names(variant_paths),invert=T)]
if(length(setdiff(names(variant_paths),names(variant_paths_filt)))>0){
    path_remove<-setdiff(names(variant_paths),names(variant_paths_filt))
    graph_paths<-graph_paths[grep(paste(path_remove,collapse="|"),graph_paths$V2,invert=T),]
}
variant_paths<-variant_paths_filt
#
allele_paths<-pathref[grep("^IG|^OG|^TR",names(pathref))]
names(allele_paths)<-gsub("\\*","::",names(allele_paths))
allele_paths_phasing_nodes<-allele_paths
pseudoploidy=length(allele_paths)
read_paths<-pathref[grep("variant|grch|chm|^IG|^OG|^TR",names(pathref),invert=T)]
if(pseudoploidy>1){
#    allele_paths<-allele_paths[grep("freq=0.0",names(allele_paths),invert=T)]
    for(allele in names(allele_paths)){
        alt_alleles<-allele_paths[grep(allele,names(allele_paths),invert=T)]
        allele_nodes_tmp<-setdiff(allele_paths[[allele]],unname(unlist(alt_alleles)))
        # if no singularly distinct nodes - use nodes that differentiate the allele from each other allele and rely on table...
        if(length(allele_nodes_tmp)==0){
            for(alt_allele_tmp in names(alt_alleles)){
                allele_nodes_tmp<-c(allele_nodes_tmp,setdiff(allele_paths[[allele]],alt_alleles[[alt_allele_tmp]]))
            }
        } else {
            allele_nodes_tmp<-setdiff(allele_paths[[allele]],unname(unlist(alt_alleles)))
        }
        allele_paths_phasing_nodes[[allele]]<-as.vector(allele_nodes_tmp)
    }
}
if(length(variant_paths)>0){
    if(length(variant_paths)>1){
        multivar=TRUE
        var_df <- do.call(rbind, lapply(names(variant_paths), function(name) {
            filtered_content <- variant_paths[[name]][-c(1, length(variant_paths[[name]]))]
            data.frame(variant = name,
                node = filtered_content,
                stringsAsFactors = FALSE)
        }))
        phasing_reads<-read_paths[sapply(read_paths, function(vec) any(var_df$node %in% vec))]
        # find co-occurring variants
        var_df$cooc_var<-"NA"
        for(varcomb in 1:ncol(data.frame(combn(var_df$node,2)))){
            cooc_nodes<-data.frame(combn(var_df$node,2))[,varcomb]
            cooc_reads<-read_paths[sapply(read_paths, function(vec) all(cooc_nodes %in% vec))]
            if(length(cooc_reads)>0){
                var_df$cooc_var[var_df$node %in% cooc_nodes]<-paste0(sort(var_df$variant[var_df$node %in% cooc_nodes]),collapse=" ")
            }
        }
        if(all(var_df$cooc_var=="NA")){
            print("No variants phased")
            multivar_indep=TRUE
            varsame_allele=FALSE
        }
        if(all(var_df$cooc_var!="NA")){
            print("All variants phased")
            varpairs<-unique(var_df$cooc_var[grep("^NA$",var_df$cooc_var,invert=T)])
            multivar_indep=FALSE
            varsame_allele=TRUE
        } else if(any(var_df$cooc_var!="NA")){
            print("Some variants phased")
            varpairs<-unique(var_df$cooc_var[grep("^NA$",var_df$cooc_var,invert=T)])
            multivar_indep=FALSE
            varsame_allele=FALSE
        }
    } else {
        multivar=FALSE
        multivar_indep=FALSE
        varsame_allele=TRUE
        var_df <- do.call(rbind, lapply(names(variant_paths), function(name) {
            filtered_content <- variant_paths[[name]][-c(1, length(variant_paths[[name]]))]
            data.frame(variant = name,
                node = filtered_content,
                stringsAsFactors = FALSE)
        }))
        var_df$cooc_var<-"NA"
    }
}
# if there are variants - haplotag them to an allele and update paths
#
linked_vars_total<-vector();
if(length(variant_paths)>0){
    phasing_df<-data.frame(matrix(0, ncol = length(variant_paths), nrow = length(allele_paths)))
    colnames(phasing_df)<-names(variant_paths)
    rownames(phasing_df)<-names(allele_paths)
    phasing_df_alleles<-data.frame(matrix(0, ncol = 1, nrow = length(allele_paths)))
    colnames(phasing_df_alleles)<-"phased_variants"
    rownames(phasing_df_alleles)<-names(allele_paths)
    phasing_df_alleles$allele_path<-names(allele_paths)
    phasing_df_alleles$indep<-0
    hap=0
    phasing_df_alleles$hap<-hap
    # iterate over variant paths
    variant_depths<-as.numeric(gsub(".*:f:","",gsub("#.$","",names(variant_paths))));
    allele_depth<-as.numeric(gsub(".*_","",gsub("x_freq.*","",names(allele_paths))));
    names(allele_depth)<-gsub(":path.*","",names(allele_paths));
    variant_processing_order<-order(variant_depths,decreasing=T);
    for(variant in variant_processing_order){
        print(names(variant_paths)[variant])
      #  if(names(variant_paths)[variant] %in% linked_vars_total){
      #      print("Variant already phased")
      #      next
      #  }
        variant_nodes<-variant_paths[[names(variant_paths)[variant]]];
        variant_nodes<-variant_nodes[-c(1,length(variant_nodes))]
        if(multivar==TRUE & multivar_indep==FALSE & var_df[var_df$node==variant_nodes,]$cooc_var!="NA"){
            # parse var_df table to indetify full set of linked reads
            linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")));
            variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
            while(length(setdiff(var_df[var_df$variant %in% c(unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))),]$variant,linked_vars)>0)){
                linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))
                variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
            }
            phasing_reads<-read_paths[sapply(read_paths, function(vec) any(variant_nodes %in% vec))]
            for(allele in names(allele_paths_phasing_nodes)){
                phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) any(unname(unlist(allele_paths_phasing_nodes[allele])) %in% vec))]
                if(length(phasing_reads_tmp)>0){
                    phasing_df[allele,names(variant_paths)[variant]]<-length(phasing_reads_tmp)
                }
                # update the phasing_df_alleles table
                var_ids_temp<-var_df[var_df$node %in% variant_nodes,]$variant
                if(length(intersect(var_ids_temp,phasing_df_alleles$phased_variant))==0){
                    colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                    if(nrow(phasing_df)>1 & any(phasing_df[,colindex]>0)){
#                        allele_w_max_vardepth <- rownames(phasing_df)[which.max(phasing_df[,colindex])]
                        subset_data<-as.matrix(phasing_df[, colindex])
                        max_row_col <- arrayInd(which.max(subset_data), dim(subset_data))
                        max_row_index <- max_row_col[1]
                        allele_w_max_vardepth <- rownames(phasing_df)[max_row_index]
                        hap=hap+1
                        var_ids_temp<-var_ids_temp[var_ids_temp %in% names(variant_paths)]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=hap)),check.names=F)
                    } else if(nrow(phasing_df)==1 & any(phasing_df[,colindex]>0)){
                        allele_w_max_vardepth <- rownames(phasing_df)
                        hap=hap+1
                        var_ids_temp<-var_ids_temp[var_ids_temp %in% names(variant_paths)]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=hap)),check.names=F)
                    } else {
                        print("Unable to phase variants with a specific allele")
                        hap=hap+1
                        var_ids_temp<-var_ids_temp[var_ids_temp %in% names(variant_paths)]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = "ambiguous", indep=0,hap=hap)),check.names=F)
                    }
                } else {
                    # check for conflicts and if not - use phasing_df_alleles
                    colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                    if(nrow(phasing_df)>1 & any(phasing_df[,colindex]>0)){
                        subset_data<-as.matrix(phasing_df[, colindex])
                        max_row_col <- arrayInd(which.max(subset_data), dim(subset_data))
                        max_row_index <- max_row_col[1]
                        allele_w_max_vardepth <- rownames(phasing_df)[max_row_index]
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    } else if(nrow(phasing_df)==1 & any(phasing_df[,colindex]>0)){
                        allele_w_max_vardepth <- rownames(phasing_df)
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    } else {
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                    }
                    if(any(phasing_df[,colindex]>0) & length(matching_hap)>1){
                        print("error - multiple haplotypes tagged for the same variant")
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = "ambiguous", indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    } else if(all(phasing_df[,colindex]==0) & length(matching_hap)==1){
                        print("Unable to phase variants with a specific allele")
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = "ambiguous", indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    }
                }
            }
            linked_vars_total<-c(linked_vars_total,linked_vars)
        } else {
            print("Link this independent variant to specific allele")
            phasing_reads<-read_paths[sapply(read_paths, function(vec) any(variant_nodes %in% vec))]
            for(allele in names(allele_paths_phasing_nodes)){
                phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) any(unname(unlist(allele_paths_phasing_nodes[allele])) %in% vec))]
                if(length(phasing_reads_tmp)>0){
                    phasing_df[allele,names(variant_paths)[variant]]<-length(phasing_reads_tmp)
                }
            }
            var_ids_temp<-var_df[var_df$node %in% variant_nodes,]$variant
            if(length(intersect(var_ids_temp,phasing_df_alleles$phased_variants))==0){
                colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                allele_w_max_vardepth <- rownames(phasing_df)
                if(nrow(phasing_df)>1){
                    allele_w_max_vardepth <- rownames(phasing_df)[which.max(phasing_df[,colindex])]
                }
                hap=hap+1
                phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=1,hap=hap)),check.names=F)
            }
        }
    }
    phasing_df_alleles<-phasing_df_alleles[grep("^0$",phasing_df_alleles$hap,invert=T),]
    phasing_df_alleles<-phasing_df_alleles[!duplicated(phasing_df_alleles$phased_variants),]
    # single allele inferred - but 1+ variants present - use depth to differentiate homozygous variant vs. 1 exact + 1 variant allele
#    if(length(allele_paths)==1){
#        if(sum(allele_depth)>5 & length(variant_depths)==1 & max(variant_depths)/sum(allele_depth)<0.667){ # HET
#            allele_graph_path<-graph_paths[grep("^IG|^OG|^TR",graph_paths$V2),];
#            allele_graph_path<-rbind(allele_graph_path,allele_graph_path);
#            allele_graph_path[1,]$V2<-sub("path.*?_","exact_",allele_graph_path[1,]$V2);
#            graph_paths<-rbind(graph_paths,allele_graph_path[1,])
#        } else if(sum(allele_depth)>5 & multivar==TRUE){
#            print("Multiple variants within the contex of a single inferred allele");
#            if(multivar_indep==TRUE){
#                print("Multi-variants are unable to be phased and may occur across all alleles");
#            } 
#            if(varsame_allele==TRUE & max(variant_depths)/sum(allele_depth)<0.667){
#                print("All variants are phased and thus on one of the 2+ alleles - other allele=exact match");
#                allele_graph_path<-graph_paths[grep("^IG|^OG|^TR",graph_paths$V2),];
#                allele_graph_path<-rbind(allele_graph_path,allele_graph_path);
#                allele_graph_path[1,]$V2<-sub("path.*?_","exact_",allele_graph_path[1,]$V2);
#                graph_paths<-rbind(graph_paths,allele_graph_path[1,])
#            }
#        } else {
#            print("Multiple independent variants or sufficient depth to be homozygous")
#            # note - for instsances of >2 haplotypes, unable to say exact match
#            # could use depth of allele-specific nodes within same positions as variants to determine presence of exact allele sequence
#        }
#    }
#
    print("Multiple alleles inferred")
    graph_paths_update<-graph_paths[-which(graph_paths$V2 %in% phasing_df_alleles$phased_variants),]
    # use - phasing_df_alleles -- phasing_df_alleles$haps
    for(allelic_backbone_tmp in unique(phasing_df_alleles$allele_path)){
        if(nrow(phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,])>1){
            phasing_df_alleles$phased_variants
            variant_depths_tmp<-as.numeric(gsub(".*:f:","",gsub("#.$","",phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,]$phased_variants)));
            allele_depth_tmp<-as.numeric(gsub(".*_","",gsub("x_freq.*","",unique(phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,]$allele_path))))
            if(all(variant_depths_tmp <=1.33*allele_depth_tmp) & all(variant_depths_tmp >=0.67*allele_depth_tmp)){
                print("Multiple unphased variants linked to the same allelic backbone within 33% of allele depth - merging")
                updating_hap<-min(phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,]$hap)
                phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,]$hap<-updating_hap
            } else {
                print("Multiple unphased variants linked to the same allelic backbone of various depths - resulting in allelic/variant ambiguity")
                phasing_df_alleles[phasing_df_alleles$allele_path==allelic_backbone_tmp & phasing_df_alleles$indep==1,]$allele_path<-"ambiguous"
                phasing_df_alleles[phasing_df_alleles$allele_path=="ambiguous" & phasing_df_alleles$indep==1,]$indep<-0
            }
        }
    }
    for(hap in unique(phasing_df_alleles$hap)){
        phased_df_tmp<-phasing_df_alleles[phasing_df_alleles$hap==hap,]
        allele_backbone<-unique(phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path)
        allele_backbone<-allele_backbone[!is.na(allele_backbone)]
        allele_backbone_pattern<-gsub(".*::","",allele_backbone)
        if(length(allele_backbone)>1){
            print("error - multiple alleles tagged for the same haplotype")
        } else {
            variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
            if(length(variant_replace)>1){
                var_length<-length(variant_replace)
                #avg_depth<-round(mean(as.numeric(gsub("variant_DP:f:","",gsub("#.*","",variant_replace)))),2)
                max_depth<-round(max(as.numeric(gsub("variant_DP:f:","",gsub("#.*","",variant_replace)))),2)
                variant_replace<-paste0("variant_DP:f:",max_depth,"#",hap,"_")
            }
            allele_update<-sub("path.*?_",variant_replace,allele_backbone)
            if(allele_update=="ambiguous"){
                print("Unable to link variant(s) with alleles")
                variant_graph_paths<-graph_paths[graph_paths$V2 %in% phased_df_tmp[phased_df_tmp$allele_path=="ambiguous",]$phased_variants,]
                variant_graph_paths$V2<-gsub("^.*#1#","",gsub("_variant_DP",":ambiguous_variant",variant_graph_paths$V2))
                # could potentially condense the variant space - if variants are phased and share immediate nodes with each other
                # cant say there are any exact matches to reference alleles for this individual - set to ambiguous
                ambiguous_allele_paths<-graph_paths_update[graph_paths_update$V2 %in% names(pathref[grep("^IG|^OG|^TR",names(pathref))]),]
                ambiguous_allele_paths$V2<-sub("path.*?_","ambiguous_exact_",ambiguous_allele_paths$V2)
                graph_paths_update<-rbind(graph_paths_update,variant_graph_paths)
                # remove allele pahts
                if(length(which(graph_paths_update$V2 %in% names(pathref[grep("^IG|^OG|^TR",names(pathref))])))>0){
                    graph_paths_update<-graph_paths_update[-which(graph_paths_update$V2 %in% names(pathref[grep("^IG|^OG|^TR",names(pathref))])),]
                    graph_paths_update<-rbind(graph_paths_update,ambiguous_allele_paths)
                }
            } else {
                allele_graph_path_replace<-graph_paths[grep(allele_backbone_pattern,graph_paths$V2),];
                variant_graph_path<-graph_paths[which(graph_paths$V2 %in% phased_df_tmp$phased_variants),];
                allele_graph_path_replace$V2<-sub("path.*?_",variant_replace,allele_graph_path_replace$V2)
                for(new_nodes_iterate in variant_graph_path$V3){
                    new_nodes<-gsub("\\+|\\-","",new_nodes_iterate)
                    new_nodes<-as.numeric((strsplit(paste0(new_nodes,collapse=","),split=",")[[1]]))
                    old_nodes<-allele_graph_path_replace$V3
                    old_nodes<-gsub("\\+|\\-","",allele_graph_path_replace$V3)
                    old_nodes<-as.numeric((strsplit(paste0(old_nodes,collapse=","),split=",")[[1]]))
                    old_nodes_replace<-old_nodes[which(old_nodes==new_nodes[1]):which(old_nodes==new_nodes[length(new_nodes)])]
                    old_nodes_replace<-setdiff(old_nodes_replace,new_nodes)
                    new_nodes_replace<-setdiff(new_nodes,old_nodes)
                    old_pattern_pos<-paste0(",",old_nodes_replace,"\\+|^",old_nodes_replace,"\\+")
                    old_pattern_neg<-paste0(",",old_nodes_replace,"\\-|^",old_nodes_replace,"\\-")
                    if(length(grep(old_pattern_pos,allele_graph_path_replace$V3))>0){
                        allele_graph_path_replace$V3<-gsub(old_pattern_pos,paste0(new_nodes_replace,"\\+"),allele_graph_path_replace$V3)
                        allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                        allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                    } else {
                        allele_graph_path_replace$V3<-gsub(old_pattern_neg,paste0(new_nodes_replace,"\\-"),allele_graph_path_replace$V3)
                        allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                        allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                    }
                }
                #remove appropriate backbone path/allele
                graph_paths_update<-rbind(graph_paths_update,allele_graph_path_replace)
            }
        }
    }
    for(orig_allele in unique(gsub(":variant_DP.*","",graph_paths_update$V2[grep(":variant_DP", graph_paths_update$V2)]))){
        print(orig_allele)
        # if sum of variant haplotype depths <0.667*sum of backbone allele depth - then ambiguous exact
        backbone_depth_tmp<-graph_paths_update[grep(gsub("\\*","\\\\*",orig_allele),graph_paths_update$V2)]$V2;
#        backbone_depth_tmp<-backbone_depth_tmp[grep("exact",backbone_depth_tmp,invert=T)];
        backbone_depth_tmp_allele<-as.numeric(gsub(".*_","",gsub("x_freq.*","",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=T)])));
        backbone_depth_tmp_variants<-sum(as.numeric(gsub("#.*","",gsub(".*variant_DP:f:","",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=F)]))));
        if(length(backbone_depth_tmp_allele)>0){
            if(backbone_depth_tmp_allele>5 & backbone_depth_tmp_variants/backbone_depth_tmp_allele<0.667){ # ambiguous_exact
                differential_depth<-backbone_depth_tmp_allele-backbone_depth_tmp_variants
                if(varsame_allele==TRUE){
                    ambiguous_path_old<-sub("path.*?_","exact_",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=T)])
                    ambiguous_path_new<-sub("_\\d+x_freq", paste0("_", differential_depth, "x_freq"), ambiguous_path_old)                    
                } else {
                    ambiguous_path_old<-sub("path.*?_","ambiguous_exact_",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=T)])
                    ambiguous_path_new<-sub("_\\d+x_freq", paste0("_", differential_depth, "x_freq"), ambiguous_path_old)
                }
                ambiguous_exact_newpath<-graph_paths_update[which(graph_paths_update$V2==backbone_depth_tmp[grep("path",gsub(".*\\*","",backbone_depth_tmp))]),]
                graph_paths_update<-graph_paths_update[-which(graph_paths_update$V2==backbone_depth_tmp[grep("path",gsub(".*\\*","",backbone_depth_tmp))]),]
                ambiguous_exact_newpath$V2<-ambiguous_path_new
                graph_paths_update<-rbind(graph_paths_update,ambiguous_exact_newpath)
            }
            if(backbone_depth_tmp_allele>5 & backbone_depth_tmp_variants/backbone_depth_tmp_allele>0.667){ # no ref allele - remove and use variant allele instead
                orig_allele_old<-backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=T)]
                variant_path_new<-backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=F)]
                graph_paths_update<-graph_paths_update[-which(graph_paths_update$V2==orig_allele_old),]
            }
        }
    }
    # for reference alleles with no matched variants - set to exact match
    for(orig_allele in unique(graph_paths_update$V2[grep(":path", graph_paths_update$V2)])){
        print("Variants found and phased to appropriate allele - remaining alleles are exact matches to reference db")
        print(orig_allele)
        replacement_graph_allele_path<-graph_paths_update[which(graph_paths_update$V2==orig_allele)]
        replacement_graph_allele_path$V2<-sub("path.*?_","exact_",orig_allele);
        graph_paths_update<-graph_paths_update[-which(graph_paths_update$V2==orig_allele)]
        graph_paths_update<-rbind(graph_paths_update,replacement_graph_allele_path)
    }
} else {
    print("No variants identified for this sample:gene (or variants found with insufficient read depth) - setting to exact reference allele matches!")
    graph_paths_update<-graph_paths
    for(orig_allele in unique(graph_paths_update$V2[grep(":path", graph_paths_update$V2)])){
        print(orig_allele)
        replacement_graph_allele_path<-graph_paths_update[which(graph_paths_update$V2==orig_allele)]
        replacement_graph_allele_path$V2<-sub("path.*?_","exact_",orig_allele);
        graph_paths_update<-graph_paths_update[-which(graph_paths_update$V2==orig_allele)]
        graph_paths_update<-rbind(graph_paths_update,replacement_graph_allele_path)
    }
}
###############################################################################################################
fwrite(graph_paths_update,file=gsub(".plines",".cleaned_paths",path_file),quote=F,col.names=F,row.names=F,sep="\t")  #
###############################################################################################################
#
#final_sweep<-graph_paths[grep("^IG|^OG|^TR",graph_paths$V2),]
#for(i in final_sweep$V2){
#    replacement_graph_allele_path<-sub("path.*?_","exact_",i);
#    graph_paths[graph_paths$V2==i,]$V2<-replacement_graph_allele_path
#}

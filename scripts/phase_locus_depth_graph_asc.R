#!/usr/bin/env Rscript
workingdir=getwd()
# -- use this to find nodes not shared with alleles we care about  -- then vg find to extract reads that preferentially map to those non-target alleles and remove them from consideration
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
args<-commandArgs(TRUE)
path_file=args[1]
# asc-based inference
#path_file="NA21130.wg_immunovar.IGHV_F12-G65.genome_graph_ref.augmented.depth.plines"
graph_paths=fread(path_file,header=F)
link_edges<-fread(gsub(".plines",".llines",path_file),header=F)
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
deletion_paths<-pathref[grep("deletion",names(pathref))]
#
aln_nodes<-graph_paths[grep("IG|TR|grch38|chm13",graph_paths$V2,invert=T),c(2,3)]
#
variant_paths_filt<-variant_paths[grep(":f:1#|:f:2#|:f:1\\.|:f:2\\.",names(variant_paths),invert=T)]
if(length(setdiff(names(variant_paths),names(variant_paths_filt)))>0){
    path_remove<-setdiff(names(variant_paths),names(variant_paths_filt))
    graph_paths<-graph_paths[grep(paste(path_remove,collapse="|"),graph_paths$V2,invert=T),]
}
variant_paths<-variant_paths_filt
#
deletion_paths_filt<-deletion_paths[grep(":f:1#|:f:2#",names(deletion_paths),invert=T)]
if(length(setdiff(names(deletion_paths),names(deletion_paths_filt)))>0){
    path_remove<-setdiff(names(deletion_paths),names(deletion_paths_filt))
    graph_paths<-graph_paths[grep(paste(path_remove,collapse="|"),graph_paths$V2,invert=T),]
}
deletion_paths<-deletion_paths_filt
if(length(deletion_paths)>0){
    variant_paths<-c(variant_paths,deletion_paths)
}
#
allele_paths<-pathref[grep("^IG|^OG|^TR",names(pathref))];
names(allele_paths)<-gsub("\\*","::",names(allele_paths));
allele_paths_phasing_nodes<-allele_paths;
pseudoploidy=length(allele_paths);
read_paths<-pathref[grep("variant|deletion|grch|chm|^IG|^OG|^TR",names(pathref),invert=T)];
generate_pairs <- function(vec) {
  combn(vec, 2, simplify = TRUE)  # Generate pairs of co-occurring values
}
convert_columns_to_list <- function(df) {
  return(lapply(df, function(col) as.vector(col)))
}
if(pseudoploidy>1){
#    allele_paths<-allele_paths[grep("freq=0.0",names(allele_paths),invert=T)]
    for(allele in names(allele_paths)){
        print(allele)
        allele_gene<-gsub("::.*","::",allele);
        alt_alleles<-allele_paths[grep(allele,names(allele_paths),invert=T)]
        allele_nodes_tmp<-setdiff(allele_paths[[allele]],unname(unlist(alt_alleles)))
        # if no singularly distinct nodes - use combination of nodes to differentiate the allele from each other allele and rely on table...
#        if(length(allele_nodes_tmp)==0){
        print("Also using unique co-occuring pairs of nodes to aid allelic phasing")
        allele_combns<-allele_paths[grep(allele,names(allele_paths),invert=F)]
        if(length(unlist(allele_paths[grep(allele,names(allele_paths),invert=F)]))>1){
            allele_combns <- lapply(allele_combns, generate_pairs)
            allele_combns <- convert_columns_to_list(data.frame(allele_combns))
            for(iter in seq_along(allele_combns)){
                allele_combns[[iter]]<-paste0(sort(allele_combns[[iter]]),collapse=",")
            }
        }
        alt_allele_combns_combined<-list()
        for(alt_allele_tmp in names(alt_alleles)){
            if(length(unlist(allele_paths[grep(alt_allele_tmp,names(allele_paths),invert=F)]))>1){
                alt_allele_combns <- lapply(allele_paths[grep(alt_allele_tmp,names(allele_paths),invert=F)], generate_pairs)
                alt_allele_combns <- convert_columns_to_list(data.frame(alt_allele_combns))
                alt_allele_combns_combined<-c(alt_allele_combns_combined,alt_allele_combns)
                for(iter in seq_along(alt_allele_combns_combined)){
                    alt_allele_combns_combined[[iter]]<-paste0(sort(alt_allele_combns_combined[[iter]]),collapse=",")
                }
            }
        }
        allele_nodes_tmp<-unique(c(allele_nodes_tmp,unname(unlist(setdiff(allele_combns,alt_allele_combns_combined)))))
#        } 
        allele_paths_phasing_nodes[[allele]]<-as.vector(allele_nodes_tmp)
    }
}
if(length(variant_paths)>0){
    if(length(variant_paths)>1){
        multivar=TRUE
        var_df <- do.call(rbind, lapply(names(variant_paths), function(name) {
            variant_path <- variant_paths[[name]]
            # If there are only two values, keep both in the same row
            if (length(variant_path) == 2) {
                filtered_content <- paste(variant_path, collapse = ",")  # Combine the two values into a single string
            } else {
                # Otherwise, exclude the first and last values
                filtered_content <- variant_path[-c(1, length(variant_path))]
            }
            data.frame(variant = name,
                    node = filtered_content,
                    stringsAsFactors = FALSE)
        }))
#        var_df <- do.call(rbind, lapply(names(variant_paths), function(name) {
#            filtered_content <- variant_paths[[name]][-c(1, length(variant_paths[[name]]))]
#            data.frame(variant = name,
#                node = filtered_content,
#                stringsAsFactors = FALSE)
#        }))
        phasing_reads<-read_paths[sapply(read_paths, function(vec) any(unlist(strsplit(as.character(var_df$node),",")) %in% vec))]
#        phasing_reads<-read_paths[sapply(read_paths, function(vec) any(var_df$node %in% vec))]
        # find co-occurring variants
        var_df$cooc_var<-"NA"
        for(varcomb in 1:ncol(data.frame(combn(var_df$node,2)))){
            cooc_nodes<-data.frame(combn(var_df$node,2))[,varcomb]
            cooc_nodes_wdel<-unlist(strsplit(as.character(cooc_nodes),","))
            cooc_reads<-read_paths[sapply(read_paths, function(vec) all(cooc_nodes_wdel %in% vec))]
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
            variant_path <- variant_paths[[name]]
            # If there are only two values, keep both in the same row
            if (length(variant_path) == 2) {
                filtered_content <- paste(variant_path, collapse = ",")  # Combine the two values into a single string
            } else {
                # Otherwise, exclude the first and last values
                filtered_content <- variant_path[-c(1, length(variant_path))]
            }
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
    variant_depths<-as.numeric(gsub(".*:f:","",gsub("#.$|#..$","",names(variant_paths))));
    allele_depth<-as.numeric(gsub(".*_","",gsub("x_freq.*","",names(allele_paths))));
    names(allele_depth)<-gsub(":path.*","",names(allele_paths));
    variant_processing_order<-order(variant_depths,decreasing=T);
    for(variant in variant_processing_order){
        print(names(variant_paths)[variant])
        variant_nodes<-variant_paths[[names(variant_paths)[variant]]];
        if(length(variant_nodes)==3){
            variant_nodes<-variant_nodes[-c(1,length(variant_nodes))]
        } else {
            variant_nodes<-variant_nodes
        }
        linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% c(paste0(variant_nodes,collapse=","),variant_nodes),]$cooc_var,split=" ")));
#        if(multivar==TRUE & multivar_indep==FALSE & var_df[var_df$node==paste0(variant_nodes,collapse=","),]$cooc_var!="NA"){
        variant_id_tmp<-names(variant_paths)[variant]
        if(all(multivar==TRUE & multivar_indep==FALSE & var_df[which(var_df$variant==variant_id_tmp),]$cooc_var!="NA")){
            variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
            variant_nodes_strict<-variant_nodes
            while(length(setdiff(var_df[var_df$variant %in% c(unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))),]$variant,linked_vars)>0)){
                linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))
                variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
            }
#            variant_nodes<-variant_nodes_strict
            # parse var_df table to indetify full set of linked reads
            #linked_vars<-unique(unlist(strsplit(var_df[var_df$node==paste0(variant_nodes,collapse=","),]$cooc_var,split=" ")));
            in_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[1])))
            out_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[length(x)])))
            reqd_nodes<-unique(c(in_node,out_node))
            phasing_reads<-read_paths[sapply(read_paths, function(vec) any(unlist(strsplit(as.character(variant_nodes),split=",")) %in% vec))]
            candidate_alleles<-vector()
            for(allele in names(allele_paths_phasing_nodes)){
                print(allele)
#                allele_paths[[allele]]
                if(all(reqd_nodes %in% unlist(allele_paths[allele]))){
                    candidate_alleles<-unique(c(candidate_alleles,allele))
                }
            }
            if(length(candidate_alleles)>1){
                linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% c(paste0(variant_nodes_strict,collapse=","),variant_nodes_strict),]$cooc_var,split=" ")));
                variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
                while(length(setdiff(var_df[var_df$variant %in% c(unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))),]$variant,linked_vars)>0)){
                    linked_vars<-unique(unlist(strsplit(var_df[var_df$node %in% variant_nodes,]$cooc_var,split=" ")))
                    variant_nodes<-var_df[var_df$variant %in% linked_vars,]$node;
                }
                candidate_alleles_refined<-vector() # try use co-occuring nodes to find single backbone allele
                for(allele in candidate_alleles){
                    phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) all(unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=T)]) %in% vec))]
                    for(node_pairs in unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=F)])){
                        node_pairs<-strsplit(node_pairs,split=",")[[1]]
#                        print(phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                        phasing_reads_tmp<-unique(c(phasing_reads_tmp,phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))]))
                    }
                    if(length(phasing_reads_tmp)>0){
                        in_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[1])))
                        out_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[length(x)])))
                        reqd_nodes<-unique(c(in_node,out_node))
                        if(sum(reqd_nodes %in% unlist(allele_paths[allele]))/(length(reqd_nodes))>0.9){
                            candidate_alleles_refined<-unique(c(candidate_alleles_refined,allele))
                        }
                    }
                }
                candidate_alleles<-candidate_alleles_refined
                print(paste0(candidate_alleles," contains topology consistent with variant call"));
            }
            if(length(candidate_alleles)==0){
                # allele should have flanking nodes in common with variant nodes
                candidate_alleles<-vector() # try use co-occuring nodes to find single backbone allele
                for(allele in names(allele_paths_phasing_nodes)){
                    phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) all(unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=T)]) %in% vec))]
                    for(node_pairs in unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=F)])){
                        node_pairs<-strsplit(node_pairs,split=",")[[1]]
#                        print(phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                        phasing_reads_tmp<-c(phasing_reads_tmp,phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                    }
                    if(length(phasing_reads_tmp)>0){
                        candidate_alleles<-unique(c(candidate_alleles,allele))
                    }
                }
                if(length(candidate_alleles)>0){
                    print(paste0(candidate_alleles," contains topology consistent with variant call"));
                }
            }
            for(allele in candidate_alleles){
                #print(paste0("Candidate allele: ", allele));
                phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) any(unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=T)]) %in% vec))]
                for(node_pairs in unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=F)])){
                    node_pairs<-strsplit(node_pairs,split=",")[[1]]
#                    print(phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                    phasing_reads_tmp<-unique(c(phasing_reads_tmp,phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))]))
                }
                var_ids_temp<-var_df[var_df$node %in% variant_nodes,]$variant
#                phasing_df[allele,names(variant_paths)[variant]]<-length(phasing_reads_tmp)
                # update the phasing_df_alleles table
                if(length(intersect(var_ids_temp,phasing_df_alleles$phased_variant))==0){
                    phasing_df[allele,unique(var_ids_temp)]<-length(phasing_reads_tmp);
                    colindex<-which(colnames(phasing_df) %in% var_ids_temp);
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
                        print("Unable to phase variants with a specific allele 1")
                        hap=hap+1
                        var_ids_temp<-var_ids_temp[var_ids_temp %in% names(variant_paths)]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = "ambiguous", indep=0,hap=hap)),check.names=F)
                    }
                } else {
                    # check for conflicts and if not - use phasing_df_alleles
                    if(phasing_df[allele,names(variant_paths)[variant]]<length(phasing_reads_tmp)){
                        phasing_df[allele,names(variant_paths)[variant]]<-length(phasing_reads_tmp)
                    }
#                    colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                    colindex<-which(colnames(phasing_df) %in% names(variant_paths)[variant])
                    if(all(nrow(phasing_df)>1 & any(phasing_df[,colindex]>0))){
                        subset_data<-as.matrix(phasing_df[, colindex])
                        max_row_col <- arrayInd(which.max(subset_data), dim(subset_data))
                        max_row_index <- max_row_col[1]
                        allele_w_max_vardepth <- rownames(phasing_df)[max_row_index]
#                        if(max(phasing_df_alleles$hap)==0){
#                            # initialize hap
#                            hap=hap+1
#                            matching_hap=hap
#                        } else {
                           matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
#                        }
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    } else if(all(nrow(phasing_df)==1 & any(phasing_df[,colindex]>0))){
                        allele_w_max_vardepth <- rownames(phasing_df)
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                        vartemp_phased<-data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=matching_hap)
                        vartemp_phased<-vartemp_phased[vartemp_phased$phased_variants %in% names(variant_paths)[variant],]
                        phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, vartemp_phased),check.names=F)
                    } else {
                        matching_hap<-unique(phasing_df_alleles[phasing_df_alleles$phased_variants %in% var_ids_temp,]$hap)
                    }
                    if(all(any(phasing_df[,colindex]>0) & length(matching_hap)>1)){
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
            if(length(candidate_alleles)==0){
                var_ids_temp<-var_df[var_df$node %in% variant_nodes,]$variant
            #    if(length(intersect(var_ids_temp,phasing_df_alleles$phased_variants))==0){
                    colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                    allele_w_max_vardepth <- rownames(phasing_df)
                    if(nrow(phasing_df)>1){
                        if(sum(phasing_df[,colindex],na.rm=T)>1){
                            allele_w_max_vardepth <- rownames(phasing_df)[which.max(phasing_df[,colindex])]
                        } else {
                            allele_w_max_vardepth<-"ambiguous";
                        }
                    }
                    hap=hap+1
                    phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=0,hap=hap)),check.names=F)
            #    }
            }
            linked_vars_total<-unique(c(linked_vars_total,linked_vars))
        } else {
            print("Link this independent variant to specific allele")
            in_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[1])))
            out_node<-unlist(unname(sapply(variant_paths[var_df[var_df$node %in% variant_nodes,]$variant], function(x) x[3])))
            reqd_nodes<-unique(c(in_node,out_node))
            candidate_alleles<-vector()
            for(allele in names(allele_paths_phasing_nodes)){
                print(allele)
#                allele_paths[allele]
                if(all(reqd_nodes %in% unlist(allele_paths[allele]))){
                    candidate_alleles<-unique(c(candidate_alleles,allele))
                    print(paste0("Candidate allele: ", allele," contains topology consistent with variant call"));
                }
            }
            phasing_reads<-read_paths[sapply(read_paths, function(vec) any(variant_nodes %in% vec))]
            if(length(candidate_alleles)>0){
                for(allele in candidate_alleles){
                    phasing_reads_tmp<-phasing_reads[sapply(phasing_reads, function(vec) any(unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=T)]) %in% vec))]
                    for(node_pairs in unname(unlist(allele_paths_phasing_nodes[allele])[grep(",",unlist(allele_paths_phasing_nodes[allele]),invert=F)])){
                        node_pairs<-strsplit(node_pairs,split=",")[[1]]
#                       print(phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                        phasing_reads_tmp<-c(phasing_reads_tmp,phasing_reads[sapply(phasing_reads, function(vec) all(node_pairs %in% vec))])
                    }
#                    if(length(phasing_reads_tmp)>0){
                        phasing_df[allele,names(variant_paths)[variant]]<-length(phasing_reads_tmp)
 #                   }
                }
                var_ids_temp<-var_df[var_df$node %in% c(paste0(variant_nodes,collapse=","),variant_nodes),]$variant
                if(length(intersect(var_ids_temp,phasing_df_alleles$phased_variants))==0){
                    colindex<-which(colnames(phasing_df) %in% var_ids_temp)
                    allele_w_max_vardepth <- rownames(phasing_df)
                    if(nrow(phasing_df)>1){
                        if(sum(phasing_df[,colindex],na.rm=T)>1){
                            if(length(grep(max(phasing_df[,colindex]),phasing_df[,colindex]))==1){
                                allele_w_max_vardepth <- rownames(phasing_df)[which.max(phasing_df[,colindex])]
                            } else {
                                allele_w_max_vardepth<-"ambiguous"; # ties
                                print("Multiple reference alleles with equivalent read support")
                            }
                        } else {
                            allele_w_max_vardepth<-"ambiguous";
                        }
                    }
                    hap=hap+1
                    phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=1,hap=hap)),check.names=F)
                }
            } else {
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
                        if(sum(phasing_df[,colindex],na.rm=T)>1){
                            allele_w_max_vardepth <- rownames(phasing_df)[which.max(phasing_df[,colindex])]
                        } else {
                            allele_w_max_vardepth<-"ambiguous";
                        }
                    }
                    hap=hap+1
                    phasing_df_alleles <- data.frame(rbind(phasing_df_alleles, data.frame(phased_variants = var_ids_temp, allele_path = allele_w_max_vardepth, indep=1,hap=hap)),check.names=F)
                }
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
#    print("Multiple alleles inferred")
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
    #function to check for duplication in allele sequence
    trim_duplicates <- function(string) {
        elements <- unlist(strsplit(string, ","))
        # Loop through the elements to check for the first duplication
        seen <- character(0)  # Vector to track seen elements
        for (i in 1:length(elements)) {
        if (elements[i] %in% seen) {
            # Identify the sequence before and after the duplicate
            initial_sequence <- elements[1:(i-1)]
            initial_sequence_collapsed<-paste0(initial_sequence,collapse=",")
            remaining_sequence <- elements[(i):length(elements)]
            remaining_sequence_collapsed<-paste0(remaining_sequence,collapse=",")
            if(startsWith(initial_sequence_collapsed, remaining_sequence_collapsed)){
            #print("Duplicate sequence found - excluding duplicate region")
            return(initial_sequence_collapsed)
            break
            }
        }
        seen <- c(seen, elements[i])
        }
        # Return the original string if no duplicates found
        return(string)
    }
    # split haplotypes if necessary
    for(hap in unique(phasing_df_alleles$hap)){
        phased_df_tmp<-phasing_df_alleles[phasing_df_alleles$hap==hap,]
        allele_backbone<-unique(phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path)
        allele_backbone<-allele_backbone[!is.na(allele_backbone)]
        allele_backbone_pattern<-gsub(".*::","",allele_backbone)
        if(length(allele_backbone)>1){
            print("Multiple alleles tagged for the same haplotype - Splitting or attempting to use the backbone with consistently high read support")
            best_allele<-vector()
            for(iter in 1:ncol(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants])){
                best_allele<-c(best_allele,which.max(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants][,iter]))
            }
            if(length(unique(best_allele))==1){
                allele_backbone<-rownames(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants])[unique(best_allele)]
                allele_backbone<-gsub("\\.1$|\\.2$","",allele_backbone)
                allele_backbone_pattern<-gsub(".*::","",allele_backbone)
                variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
            } else {
                for(splitting_hap in allele_backbone[2:length(allele_backbone)]){
                    phasing_df_alleles[phasing_df_alleles$allele_path==splitting_hap,]$hap<-max(phasing_df_alleles$hap)+1
                }
                allele_backbone<-allele_backbone[1]
                allele_backbone_pattern<-allele_backbone_pattern[1]
                phased_df_tmp<-phasing_df_alleles[phasing_df_alleles$hap==hap,]
                variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
            }
        }
    }
    for(hap in unique(phasing_df_alleles$hap)){
        phased_df_tmp<-phasing_df_alleles[phasing_df_alleles$hap==hap,]
        allele_backbone<-unique(phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path)
        allele_backbone<-allele_backbone[!is.na(allele_backbone)]
        allele_backbone_pattern<-gsub(".*::","",allele_backbone)
        if(length(allele_backbone)>1){
            print("Multiple alleles tagged for the same haplotype - Splitting or attempting to use the backbone with consistently high read support")
            best_allele<-vector()
            for(iter in 1:ncol(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants])){
                best_allele<-c(best_allele,which.max(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants][,iter]))
            }
            if(length(unique(best_allele))==1){
                allele_backbone<-rownames(phasing_df[phasing_df_alleles[phasing_df_alleles$hap==hap,]$allele_path,phasing_df_alleles[phasing_df_alleles$hap==hap,]$phased_variants])[unique(best_allele)]
                allele_backbone<-gsub("\\.1$|\\.2$","",allele_backbone)
                allele_backbone_pattern<-gsub(".*::","",allele_backbone)
                variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
            } else {
                for(splitting_hap in allele_backbone[2:length(allele_backbone)]){
                    phasing_df_alleles[phasing_df_alleles$allele_path==splitting_hap,]$hap<-max(phasing_df_alleles$hap)+1
                }
                allele_backbone<-allele_backbone[1]
                allele_backbone_pattern<-allele_backbone_pattern[1]
                phased_df_tmp<-phasing_df_alleles[phasing_df_alleles$hap==hap,]
                variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
#                phasing_df_alleles[phasing_df_alleles$allele_path %in% allele_backbone & phasing_df_alleles$hap==hap,]$allele_path<-"ambiguous"
#                phasing_df_alleles[phasing_df_alleles$allele_path=="ambiguous",]$indep<-0
            }
        } else {
            variant_replace<-paste0(gsub(".*_variant","variant",gsub("#1","",phased_df_tmp$phased_variant)),"#1_")
            variant_replace<-paste0(gsub(".*_deletion","deletion",gsub("#1","",variant_replace)),"#1_")
        }
        if(length(variant_replace)>1){
            var_length<-length(variant_replace)
            avg_depth<-round(mean(as.numeric(gsub("#.*","",gsub("variant_DP:f:|deletion_DP:f:","",gsub("_#.*","",variant_replace))))),2)
            variant_replace<-paste0("variant_DP:f:",avg_depth,"#",hap,"_")
            #max_depth<-round(max(as.numeric(gsub("#.*","",gsub("variant_DP:f:|deletion_DP:f:","",gsub("_#.*","",variant_replace))))),2)
            #variant_replace<-paste0("variant_DP:f:",max_depth,"#",hap,"_")
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
            check_duplication<-trim_duplicates(allele_graph_path_replace$V3)
            if(check_duplication!=allele_graph_path_replace$V3){
                print("Duplicate sequence found - excluding duplicate region")
                allele_graph_path_replace$V3<-check_duplication
            }
            variant_graph_path<-graph_paths[which(graph_paths$V2 %in% phased_df_tmp$phased_variants),];
            allele_graph_path_replace$V2<-sub("path.*?_",variant_replace,allele_graph_path_replace$V2)
#            for(new_nodes_iterate in variant_graph_path[order(as.numeric(gsub("#.*$","",gsub("^.*_DP:f:","",variant_graph_path$V2))),decreasing=T),]$V3){ # depth order
            for(new_nodes_iterate in variant_graph_path[order(as.numeric(gsub("^.*#","",variant_graph_path$V2)),decreasing=F),]$V3){ # variant order
                #print(new_nodes_iterate)
                new_nodes<-gsub("\\+|\\-","",new_nodes_iterate)
                new_nodes<-as.numeric((strsplit(paste0(new_nodes,collapse=","),split=",")[[1]]))
                old_nodes<-allele_graph_path_replace$V3
                old_nodes<-gsub("\\+|\\-","",allele_graph_path_replace$V3)
                old_nodes<-as.numeric((strsplit(paste0(old_nodes,collapse=","),split=",")[[1]]))
                if(all(c(new_nodes[1],new_nodes[length(new_nodes)]) %in% old_nodes)){
                    old_nodes_replace<-old_nodes[which(old_nodes==new_nodes[1]):which(old_nodes==new_nodes[length(new_nodes)])]
                } else if(all(c(new_nodes[1],new_nodes[length(new_nodes)]) %in% as.numeric(unlist(strsplit(gsub("\\+|\\-","",graph_paths[grep(allele_backbone_pattern,graph_paths$V2),]$V3),split=","))))){
                    orig_old_nodes<-as.numeric(unlist(strsplit(gsub("\\+|\\-","",graph_paths[grep(allele_backbone_pattern,graph_paths$V2),]$V3),split=",")))
                    old_nodes_replace<-orig_old_nodes[which(orig_old_nodes==new_nodes[1]):which(orig_old_nodes==new_nodes[length(new_nodes)])]
                } else if(any(c(new_nodes[1],new_nodes[length(new_nodes)]) %in% c(old_nodes[1],old_nodes[length(old_nodes)]))){
                    print("Variant at tip of allele sequence - not currently supported - later versions will add nodes upstream/downstream of locus");
                    #reset graph path
                    allele_graph_path_replace<-graph_paths[grep(allele_backbone_pattern,graph_paths$V2),]
                    next
                }
                redundant_var=F
                if(length(old_nodes_replace)==length(new_nodes)){
                    if(all(sort(old_nodes_replace)==sort(new_nodes))){
                        redundant_var=T
                    }
                }
                old_nodes_replace<-setdiff(old_nodes_replace,new_nodes)
                new_nodes_replace<-setdiff(new_nodes,old_nodes)
                old_pattern_pos<-paste0(",",old_nodes_replace,"\\+|^",old_nodes_replace,"\\+")
                old_pattern_neg<-paste0(",",old_nodes_replace,"\\-|^",old_nodes_replace,"\\-")
                #print(paste0(new_nodes_iterate,": part 2"))
                if(all(length(old_nodes_replace)==1 & length(new_nodes_replace)==1)){
                    # check if signs match up with link_edges - if not invert the order/signs
                    old_nodes_wsigns<-strsplit(allele_graph_path_replace$V3,split=",")[[1]][which(old_nodes %in% new_nodes[new_nodes %in% old_nodes])]
                    varlocus_wsigns<-strsplit(new_nodes_iterate,split=",")[[1]]
                    varlocus_wsigns_flip<-rev(varlocus_wsigns)
                    varlocus_wsigns_flip[grep("\\+",varlocus_wsigns)]<-gsub("\\+","-",varlocus_wsigns[grep("\\+",varlocus_wsigns)]);
                    varlocus_wsigns_flip[grep("-",varlocus_wsigns)]<-gsub("-","\\+",varlocus_wsigns_flip[grep("-",varlocus_wsigns)]);
                    if(length(grep(old_pattern_pos,allele_graph_path_replace$V3))>0){
                        if(sum(old_nodes_wsigns %in% varlocus_wsigns)>=sum(old_nodes_wsigns %in% varlocus_wsigns_flip)){
                            allele_graph_path_replace$V3<-gsub(old_pattern_pos,paste0(new_nodes_replace,"\\+"),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub("\\+\\+","+",allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                        } else {
                            allele_graph_path_replace$V3<-gsub(old_pattern_pos,paste0(new_nodes_replace,"\\-"),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub("\\+\\+","+",allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                        }
                    } else if(length(grep(old_pattern_neg,allele_graph_path_replace$V3))>0){
                        if(sum(old_nodes_wsigns %in% varlocus_wsigns)>=sum(old_nodes_wsigns %in% varlocus_wsigns_flip)){
                            allele_graph_path_replace$V3<-gsub(old_pattern_neg,paste0(new_nodes_replace,"\\-"),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub("\\-\\-","-",allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                        } else {
                            allele_graph_path_replace$V3<-gsub(old_pattern_neg,paste0(new_nodes_replace,"\\+"),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub("\\-\\-","-",allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\+",new_nodes_replace),paste0("+,",new_nodes_replace),allele_graph_path_replace$V3)
                            allele_graph_path_replace$V3<-gsub(paste0("\\-",new_nodes_replace),paste0("-,",new_nodes_replace),allele_graph_path_replace$V3)
                        }
                    } else {
                        print(paste0("Investigate me: ",new_nodes_iterate))
                    }
                } else {
                    if(redundant_var==T){
                        print("No difference from partially updated graph")
                        redundant_var=F
                        next
                    } else if(length(new_nodes_replace)==0){
                        print("Putative deletion")
                        for(node_del in old_nodes_replace){
                            candidate_path_del<-gsub(paste0(old_nodes_replace,"."),"",allele_graph_path_replace$V3)
                            candidate_path_del<-gsub("^,","",candidate_path_del);candidate_path_del<-gsub(",$","",candidate_path_del);candidate_path_del<-gsub(",,",",",candidate_path_del);
                            allele_graph_path_replace$V3<-candidate_path_del
                        #    print(allele_graph_path_replace$V3)
                        }
                    } else {
                        print("Putative insertion")
                        insert_me<-strsplit(new_nodes_iterate,split=",")[[1]]
                        #link_edges
                        if(length(grep(gsub("\\+","\\\\+",gsub(paste0(new_nodes_replace,".,",collapse=""),"",new_nodes_iterate)),allele_graph_path_replace$V3))==1){
                            allele_graph_path_replace$V3<-gsub(gsub("\\+","\\\\+",gsub(paste0(new_nodes_replace,".,"),"",new_nodes_iterate)),new_nodes_iterate,allele_graph_path_replace$V3)
                        } else {
                            print("Multi-node insertion - trying to resolve")
                            flank_node1=insert_me[1];flank_node2=insert_me[length(insert_me)];
                            allele_seq_nodes_tmp<-strsplit(allele_graph_path_replace$V3,split=",")[[1]];
                            flank_node1_pos<-grep(paste0("^",gsub("\\+|-",".",flank_node1),"$"),allele_seq_nodes_tmp);
                            flank_node2_pos<-grep(paste0("^",gsub("\\+|-",".",flank_node2),"$"),allele_seq_nodes_tmp);
                            if(all(length(flank_node1_pos)==0 & flank_node2_pos==1)){
                                print("upstream of sequence")
                            } else if(all(length(flank_node2_pos)==0 & flank_node1_pos==length(allele_seq_nodes_tmp))){
                                print("upstream of sequence")
                            } else if(flank_node1_pos<flank_node2_pos){
                                if(flank_node1_pos==1){
                                    allele_seq_nodes_tmp_p1<-allele_seq_nodes_tmp[1];
                                } else {
                                    allele_seq_nodes_tmp_p1<-allele_seq_nodes_tmp[1:(flank_node1_pos)];
                                }
                                #allele_seq_nodes_tmp_p2<-allele_seq_nodes_tmp[(flank_node1_pos+1):(flank_node2_pos-1)];
                                if(flank_node2_pos==length(allele_seq_nodes_tmp)){
                                    allele_seq_nodes_tmp_p3<-allele_seq_nodes_tmp[length(allele_seq_nodes_tmp)];
                                } else {
                                    allele_seq_nodes_tmp_p3<-allele_seq_nodes_tmp[(flank_node2_pos):length(allele_seq_nodes_tmp)];
                                }
                                insert_me_clipped<-insert_me[-c(1,length(insert_me))];
                                seq_replace<-paste0(paste0(allele_seq_nodes_tmp_p1,collapse=","),",",paste0(insert_me_clipped,collapse=","),",",paste0(allele_seq_nodes_tmp_p3,collapse=","),collapse=",")
                                seq_replace<-gsub("^,","",seq_replace);seq_replace<-gsub(",$","",seq_replace);seq_replace<-gsub(",,",",",seq_replace);
                                allele_graph_path_replace$V3<-seq_replace
                            } else if(flank_node1_pos>flank_node2_pos){
                                if(flank_node2_pos==1){
                                    allele_seq_nodes_tmp_p1<-allele_seq_nodes_tmp[1];
                                } else {
                                    allele_seq_nodes_tmp_p1<-allele_seq_nodes_tmp[1:(flank_node2_pos)];
                                }
                                #allele_seq_nodes_tmp_p2<-allele_seq_nodes_tmp[flank_node1_pos:flank_node2_pos];
                                if(flank_node1_pos==length(allele_seq_nodes_tmp)){
                                    allele_seq_nodes_tmp_p3<-allele_seq_nodes_tmp[length(allele_seq_nodes_tmp)];
                                } else {
                                    allele_seq_nodes_tmp_p3<-allele_seq_nodes_tmp[(flank_node1_pos):length(allele_seq_nodes_tmp)];
                                }
                                insert_me_clipped<-insert_me[-c(1,length(insert_me))];
                                seq_replace<-paste0(paste0(allele_seq_nodes_tmp_p1,collapse=","),",",paste0(rev(insert_me_clipped),collapse=","),",",paste0(allele_seq_nodes_tmp_p3,collapse=","),collapse=",")
                                seq_replace<-gsub("^,","",seq_replace);seq_replace<-gsub(",$","",seq_replace);seq_replace<-gsub(",,",",",seq_replace);
                                allele_graph_path_replace$V3<-seq_replace
                            } else {
                                print("Investigate me: ",new_nodes_iterate)
                            }
                        }
                    }
                }
            }
            #remove appropriate backbone path/allele
            aln_nodes_candidates<-aln_nodes[grep(gsub("\\+|-",".",allele_graph_path_replace$V3),aln_nodes$V3),]$V3
            # -- use links from reads to re-orient signs/directions of things if those edges arent in the graph already...
            # if node = 1 bp, then it doesnt matter which direction!

            # --- confirm whether this step is necessary ---    
            #if(length(aln_nodes_candidates)>0){
            #    first_node<-substr(allele_graph_path_replace$V3, 1, 2)
            #    first_node<-gsub("\\+|-","",first_node)
            #    final_node <- gsub(".*,","",new_nodes_iterate)
            #    final_node <- gsub("\\+|-","",final_node)
            #    aln_nodes_candidates<-gsub(paste0("^.*",first_node),first_node,aln_nodes_candidates)
            #    aln_nodes_candidates<-gsub(paste0(final_node,"+.*"),paste0(final_node,"+"),aln_nodes_candidates)
            #    aln_nodes_candidates<-gsub(paste0(final_node,"-.*"),paste0(final_node,"-"),aln_nodes_candidates)
            #    new_varpath_update<-names(which.max(table(aln_nodes_candidates)))
            #    if(length(new_varpath_update)>0){
            #        allele_graph_path_replace$V3<-new_varpath_update
            #    }
            #}
            graph_paths_update<-rbind(graph_paths_update,allele_graph_path_replace)
        }
    }
    for(orig_allele in unique(gsub(":variant_DP.*","",graph_paths_update$V2[grep(":variant_DP", graph_paths_update$V2)]))){
        print(orig_allele)
        # if sum of variant haplotype depths <0.667*sum of backbone allele depth - then ambiguous exact
        backbone_depth_tmp<-graph_paths_update[grep(paste0(gsub("\\*","\\\\*",orig_allele),":"),graph_paths_update$V2)]$V2;
#        backbone_depth_tmp<-backbone_depth_tmp[grep("exact",backbone_depth_tmp,invert=T)];
        backbone_depth_tmp_allele<-as.numeric(gsub(".*_","",gsub("x_freq.*","",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=T)])));
        backbone_depth_tmp_variants<-sum(as.numeric(gsub("#.*","",gsub("_#.*","",gsub(".*variant_DP:f:","",backbone_depth_tmp[grep("variant",backbone_depth_tmp,invert=F)])))),na.rm=T);
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
        graph_paths_update<-graph_paths_update[!duplicated(graph_paths_update$V2),]
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

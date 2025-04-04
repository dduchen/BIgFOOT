suppressMessages(suppressWarnings(library(data.table)))
currdir=getwd()
subdir=list.files(pattern="familywise_")
sample_id=gsub("^.*/","",currdir);
if(length(list.files(pattern=paste0(gsub("_.*$","",sample_id),".results_raw.txt")))>0){
    sample_id=gsub("_.*$","",sample_id)
} else {
    sample_id=gsub("_wg_immunovar_genotyping.*$","",sample_id)
}
if(as.numeric(grepl("franken",currdir))==1){
    graph="franken"
    } else {
    graph="wg_immunovar"
}
if(as.numeric(grepl("_pe",subdir))==1){
    aln_type="pe"
    relevant_depth_file<-list.files(path="../",pattern=paste0("^",sample_id,"\\."))
    relevant_depth_file<-relevant_depth_file[grep("pe.depth",relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep(graph,relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep("immune",relevant_depth_file)]
    relevant_depth_file<-paste0("../",relevant_depth_file)
} else {
    aln_type="ngmerged"
    relevant_depth_file<-list.files(path="../",pattern=paste0("^",sample_id))
    relevant_depth_file<-relevant_depth_file[grep("ngmerged.depth",relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep(graph,relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep("immune",relevant_depth_file)]
    relevant_depth_file<-paste0("../",relevant_depth_file)
}
igh_depth<-read.table(relevant_depth_file[grep("immune",relevant_depth_file)])
#ttn_depth<-read.table(relevant_depth_file[grep("TTN",relevant_depth_file)])
#
depth_file=list.files(pattern="depth_raw.txt")
graph_results_file<-list.files(pattern="results_raw.txt")
depth<-fread(paste0(depth_file),header=T)
depth<-depth[,1:3]
graph_results<-fread(paste0(graph_results_file),header=F)
#
if(sample_id==colnames(depth)[1]){
    #
    ttn_depth=list.files(path="../",pattern="TTN.pe.depth")
    ttn_depth<-ttn_depth[grep(sample_id,ttn_depth)]
    if(length(ttn_depth)==1){
        ttn_depth<-fread(paste0("../",ttn_depth),header=F)
        ttn_depth$gene<-"TTN";ttn_depth<-ttn_depth[,c("gene","V1","V2")]
        colnames(ttn_depth)<-c(sample_id, "mean","sd")
        depth<-rbind(ttn_depth,depth)
    }
    #
    graph_results$gene<-gsub("^.*\\.","",gsub("\\*.*","",graph_results$V1))
    graph_results$allele<-gsub(":.*$","",gsub("^.*\\*","",graph_results$V1))
    graph_results$coverage<-gsub("^.* ","",gsub(" freq.*$","",graph_results$V2))
    graph_results$coverage<-as.numeric(gsub("x","",graph_results$coverage))
    #graph_results$relative_freq<-as.numeric(gsub("^.*freq=","",graph_results$V2))
    graph_results$relative_freq<-as.numeric(gsub("^.*freq=","",gsub("_Hamming.*$","",graph_results$V3)))
    #
    graph_results$passed<-"NA"
    graph_results$gene<-gsub("_.*","",graph_results$gene)
    for(i in seq_along(unique(graph_results$gene))){
        gene<-unique(graph_results$gene)[i]
        # accounting for non-standard gene nomenclature
        print(gene)
        gene_full<-gsub("__","/",gene)
        res<-data.frame(graph_results)[graph_results$gene==gene,]
    #    dep<-data.frame(depth)[grep(paste0(gene_full,".genotyping"),depth[[1]]),]
        dep<-data.frame(depth)[grep(paste0(gene_full,"$"),depth[[1]]),]
        res2<-cbind(res,dep)
        # If allele has an estimated frequency >=10% then should include
        res2$passed<-ifelse(res2$relative_freq>=0.1,TRUE,FALSE)
        # further filtering using relative frequency + locus-specific coverage
        # if coverage is low for the allele = 2*SD below the mean, set to False - or just use hard cutoff of 3
    #    res2$passed<-ifelse(res2$coverage<(res2$mean-(2*res2$sd)),FALSE,res2$passed)
        res2$passed<-ifelse(res2$coverage<2,FALSE,res2$passed)
        # However, if the allele has 4x lower coverage than the top allele, set the lower allele to FALSE 
        res2$passed<-ifelse((res2$coverage*4)<=max(res2$coverage),FALSE,res2$passed)
        # if coverage for all alleles falls below the locus depth mean - 1 sd - use relative frequency provided depth >10
        if(all(res2$passed==F)){
            res2$passed<-ifelse(res2$coverage>=2 & res2$relative_freq>=0.2,TRUE,FALSE)
        }
        if(all(graph_results[graph_results$gene %in% res2$gene,]$gene==res2$gene)){
            graph_results[graph_results$gene %in% res2$gene,]$passed<-res2$passed
        }
    }
    other_info<-gsub("^.*/","",getwd())
    graph_results$aln_type<-ifelse(grepl("_pe_",subdir)==1,"PE","NGmerge")
    graph_results$augmented_graph<-""
    #
    deletion_candidates<-list.files(path="familywise_pe_haplotype_inference",pattern="_files.txt")
    deletion_candidates_gene<-gsub(paste0(sample_id,"_"),"",gsub("_files.txt","",deletion_candidates))
    # for asc-based inference - ignore these 'gene ids'
    deletion_candidates_gene<-deletion_candidates_gene[grep("IGHV_F",deletion_candidates_gene,invert=T)]
    deletion_candidates_gene<-deletion_candidates_gene[grep("IGLV_F",deletion_candidates_gene,invert=T)]
    deletion_candidates_gene<-deletion_candidates_gene[grep("IGKV_F",deletion_candidates_gene,invert=T)]
    depth_del<-depth;colnames(depth_del)<-c("gene","mean","sd")
    # update v low coverage
    for(gene_low_depth in intersect(depth_del$gene,graph_results[graph_results$V2=="0x",]$gene)){
        graph_results[graph_results$V2=="0x" & graph_results$gene==gene_low_depth,]$V2<-paste0(round(depth_del[depth_del$gene==gene_low_depth,]$mean,2),"x")
    }
    # not necessarily raw depth - reads can be filtered out --> off-target alignment
#    deletion_candidates_gene<-deletion_candidates_gene[-which(deletion_candidates_gene %in% depth_del[depth_del$mean>0,]$gene)]
    deletion_candidates_gene<-deletion_candidates_gene[-which(deletion_candidates_gene %in% intersect(deletion_candidates_gene,unique(graph_results$gene)))]
    graph_results_wDels<-graph_results
    for(i in seq_along(deletion_candidates_gene)){
        temp_row<-graph_results[1,];
        temp_row[,1]<-paste0(deletion_candidates_gene[i],":deletion");
        temp_row[,2]<-"0x";
        temp_row[,3]<-"freq=0.000";
        temp_row[,4]<-deletion_candidates_gene[i];
        temp_row[,5]<-"00";
        temp_row[,6:7]<-0;
        temp_row[,8]<-"FALSE";
        graph_results_wDels<-rbind(graph_results_wDels,temp_row)
    }
    if(length(list.files(pattern="alleles.fasta"))>0){
        alleles<-Biostrings::readDNAStringSet(list.files(pattern="alleles.fasta"))
        alleles_novel<-alleles[grep("variant",names(alleles))]
        alleles_novel<-alleles_novel[grep("ambiguous",names(alleles_novel),invert=T)]
        alleles_exact<-alleles[grep("exact",names(alleles))]
        alleles_exact<-alleles_exact[grep("ambiguous",names(alleles_exact),invert=T)]
        alleles_ambig<-alleles[grep("ambiguous",names(alleles))]
        alleles_ambig_var<-alleles_ambig[grep("variant",names(alleles_ambig),invert=F)]
        alleles_ambig_ref<-alleles_ambig[grep("variant",names(alleles_ambig),invert=T)]
        graph_results_wDels$novel<-""
        graph_results_wDels$fasta<-""
        graph_results_wDels<-dplyr::distinct(graph_results_wDels)
        for(seqid in unique(names(alleles_novel))){
            seqid_match<-gsub(":.*","",seqid)
            graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$novel<-seqid
            graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$fasta<-as.character(alleles_novel[seqid])[[1]]
        }
        for(seqid in unique(names(alleles_exact))){
            seqid_match<-gsub(":.*","",seqid)
            if(nchar(graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$novel)==0){
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$novel<-seqid
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$fasta<-as.character(alleles_exact[seqid])[[1]]
            } else {
                graph_results_wDels<-rbind(graph_results_wDels,graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,])
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$novel[2]<-seqid
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$fasta[2]<-as.character(alleles_exact[seqid])[[1]]
            }
        }
        # for missing due to ambiguity - use ambiguous exact alleles
        for(seqid in unique(names(alleles_ambig_ref))){
            seqid_match<-gsub(":.*","",seqid)
            if(graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$fasta==""){
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$novel<-seqid
                graph_results_wDels[graph_results_wDels$gene==gsub("\\*.*","",seqid_match) & graph_results_wDels$allele==gsub("^.*\\*","",seqid_match) ,]$fasta<-as.character(alleles_ambig_ref[seqid])[[1]]
            }
        }
    # Use coverage.tsv to obtain confidence result - based on minimum node depth
        for(gene_id in unique(graph_results_wDels$gene)){
            #print(gene_id)
            file_tmp<-list.files(path="familywise_pe_haplotype_inference",pattern=paste0(gene_id,".variant.coverage.tsv"))
            if(length(file_tmp)>0 && file.size(paste0("familywise_pe_haplotype_inference/",file_tmp))>0){
                covdat<-fread(paste0("familywise_pe_haplotype_inference/",file_tmp),header=F,sep="\t")
                for(allele in unique(covdat$V1)){
                    #print(allele)
                    if(length(grep("variant|:exact",allele))>0){
                        allele_mincov<-min(covdat[covdat$V1==allele,]$V3)
                        allele_bplen<-max(covdat[covdat$V1==allele,]$V2)
                        mincov_prop<-nrow(covdat[covdat$V1==allele & covdat$V3==allele_mincov,])/allele_bplen
                        covindex_match<-grep(gsub("#.*$","",gsub("\\*","\\\\*",allele)),graph_results_wDels$novel)
                        if(length(covindex_match)>0){
                            graph_results_wDels[covindex_match,]$augmented_graph<-paste0(allele_mincov,":",round(mincov_prop*100,2))
                        }
                    }
                }
            } else {
                print(paste0("No coverage file for: ",gene_id))
            }
        }
    }
    fwrite(graph_results_wDels,file=paste0(sample_id,"_",graph_results$aln_type[1],"_",graph,"_results_cleaned.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    #################################################################
    #could derive copy number variation/alleles using the relative abundance of observed strains + strain-specific abundance/locus depth
    # or from TTN (immunotyper-SR) or the nodes immediately adjacent to the area?
} else {
    print(paste0("Sample ID mismatch for: ",sample_id))
}
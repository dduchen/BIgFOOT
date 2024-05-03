library(data.table)
currdir=getwd()
subdir=list.files(pattern="familywise_")
sample_id=gsub("^.*/","",currdir);
sample_id=gsub("_.*$","",sample_id)
if(as.numeric(grepl("franken",currdir))==1){
    graph="franken"
    } else {
    graph="wg_immunovar"
}
if(as.numeric(grepl("_pe",subdir))==1){
    aln_type="pe"
    relevant_depth_file<-list.files(path="../",pattern=paste0(sample_id))
    relevant_depth_file<-relevant_depth_file[grep("pe.depth",relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep(graph,relevant_depth_file)]
    relevant_depth_file<-relevant_depth_file[grep("immune",relevant_depth_file)]
    relevant_depth_file<-paste0("../",relevant_depth_file)
} else {
    aln_type="ngmerged"
    relevant_depth_file<-list.files(path="../",pattern=paste0(sample_id))
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
        # If allele has an estimated frequency >=20% then should include
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
    graph_results$augmented_graph<-ifelse(grepl("_augmented_",other_info)==1,T,F)
    #
    fwrite(graph_results,file=paste0(sample_id,"_",graph_results$aln_type[1],"_",graph,"_results_cleaned.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    #################################################################
    #could derive copy number variation/alleles using the relative abundance of observed strains + strain-specific abundance/locus depth
    # or from TTN (immunotyper-SR) or the nodes immediately adjacent to the area?
} else {
    print(paste0("Sample ID mismatch for: ",sample_id))
}
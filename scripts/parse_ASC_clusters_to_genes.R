#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
args<-commandArgs(TRUE)
#dat_file="/media/dduchen/Data/1kgenomes/NA21130_wg_immunovar_genotyping/familywise_pe_haplotype_inference/../NA21130.depth_raw.asc.txt"
#asc_table_file="~/Documents/github/BIgFOOT/custom_beds/ASC_metadata.matching.tsv"
dat_file=args[1]
asc_table_file=args[2]
#
dat<-fread(dat_file,header=T)
dat_recode<-dat
dat_recode<-data.frame(dat_recode)
asc_table<-fread(asc_table_file,header=T)
for(gene in unique(unname(unlist(dat[,1])))){
    print(gene)
    if(gene %in% unique(gsub("\\*.*","",asc_table$new_allele))){
        asc_tmp<-asc_table[gsub("\\*.*","",asc_table$new_allele)==gene,]
        asc_tmp_genes<-unique(gsub("\\*.*","",gsub(".*#1#","",asc_tmp$imgt_allele)))
        asc_tmp_genes<-asc_tmp_genes[grep("_P$|_F$",asc_tmp_genes,invert=T)]

        asc_depth<-dat_recode[grep(gene,dat_recode[,1]),]
        while(length(asc_tmp_genes)>nrow(asc_depth)){
            asc_depth<-rbind(asc_depth,asc_depth[1,])
        }
        asc_depth[,1]<-asc_tmp_genes
        for(gene_id in asc_tmp_genes){
            if(gene_id %in% dat_recode[,1]){
                old_mean=dat_recode[dat_recode[,1]==gene_id,]$mean
                old_sd=dat_recode[dat_recode[,1]==gene_id,]$sd
                new_mean<-mean(c(old_mean,asc_depth$mean[asc_depth[,1]==gene_id]))
                new_sd<-mean(c(old_mean,asc_depth$sd[asc_depth[,1]==gene_id]))
                rowdat_tmp<-asc_depth[asc_depth[,1]==gene_id,]
                rowdat_tmp$mean<-new_mean
                rowdat_tmp$sd<-new_sd
                dat_recode<-rbind(dat_recode,rowdat_tmp)
            } else {
                dat_recode<-rbind(dat_recode,asc_depth[asc_depth[,1]==gene_id,])
            }
        }
    dat_recode<-dat_recode[grep(gene,dat_recode[,1],invert=T),]
    }
}
dat_recode<-dat_recode[order(dat_recode[,1]),]
col_index<-as.character(colnames(dat_recode)[1])
df_summary <- dat_recode %>%
  group_by(dat_recode[[1]]) %>%
    mutate(
    mean = ifelse(n() > 1, mean(mean, na.rm = TRUE), mean),
    sd = ifelse(n() > 1, mean(sd, na.rm = TRUE), sd)
  ) %>%
  ungroup()
dat_recode<-df_summary[!duplicated(df_summary[,1]),] 
fwrite(dat_recode,file=gsub(".txt",".recode.txt",dat_file),col.names=T,quote=F,row.names=F,sep="\t")
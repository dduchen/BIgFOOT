#!/usr/bin/env Rscript
workingdir=getwd()
pathfiles<-list.files(pattern=".pathnodes")
.libPaths( c( .libPaths(), "/vast/palmer/home.mccleary/dd392/R/x86_64-pc-linux-gnu-library/4.2", "/vast/palmer/apps/avx2/software/R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0") )
library(data.table);library(ggplot2);library(viridis);library(MetBrewer)
library(ngram) #ngram
library(tm)
library(quanteda);library(quanteda.textmodels) #
for(k in 1:length(pathfiles)){
  outprefix<-gsub(".pathnodes","",pathfiles[k])
  dat<-fread(pathfiles[k],header=F)
  dat<-dat[,c(2,3)]
  #dat<-dat[unique(grep("^I|patch",dat$V2))]
  dat<-dat[unique(grep("^IM|^IG|^OGR|^KIR|^HLA",dat$V2))]
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
  str(pathref)
  #
  save(pathref,file=paste0(outprefix,"_Path_Specific_Nodes.R"))
  #
  # split by locus if more than IG present -- make sure combined in correct order at end of script (node_weights_combined)
  genes<-unique(gsub("#.*$","",names(pathref)))
  if(length(setdiff(genes,c("IMGT","IGv2","OGRDB")))==length(genes)){
    print("Custom graph - just process nodes and weights")
    total_nodes<-as.vector(NULL)
    for(i in 1:length(pathref)){
       # Want unique nodes traversed by each path - can traverse same node multiple times
      total_nodes<-c(total_nodes,as.character(unique(pathref[[i]])))
    }
    total_nodes<-table(total_nodes) # repeat for more intuitive plot
    node_dat<-data.frame(table(total_nodes))
    node_dat$nodes<-as.numeric(node_dat$total_nodes)
    path_specific_nodes_hist_immunovar<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=nodes)) + geom_bar(stat="identity") + xlab("Number of Paths/Node") + ylab("N") +
      scale_fill_gradient(name = '',low="darkblue",high="red") + theme(legend.direction = 'horizontal', legend.position = 'bottom') +
      theme_bw() + theme(legend.position="none")
    ggsave(paste0(outprefix,"_IMGT_PathSpecificNodes_hist.svg"), path_specific_nodes_hist_immunovar,device = "svg",width=4,height=5,units='in',dpi=1000)
  ######################################################################################################
  # - make lower values weighted more + normalize
    total_nodes_wt<-max(total_nodes)/total_nodes # max num. of shared nodes=17 --> downweight to 1
    total_nodes_wt = (total_nodes_wt-min(total_nodes_wt))/(max(total_nodes_wt)-min(total_nodes_wt))
    node_weights<-data.frame(total_nodes_wt)
    colnames(node_weights)<-c("node","weight")
    all(names(total_nodes)==node_weights$node) # TRUE
    node_weights$paths<-total_nodes
    node_weights$unscaled_weight<-max(total_nodes)/node_weights$paths
    # ensure numeric values treated appropriately
    node_weights$node2<-as.numeric(as.character(node_weights$node))
    node_weights<-node_weights[order(node_weights$node2),]
    node_weights$node<-node_weights$node2
    node_weights$node2<-NULL
    #fwrite(node_weights,file=paste0(outprefix,"_NodeWeightsPaths_Reference.csv"),quote=F,sep=",",col.names=T,row.names=F)
    ##############
  # remake distribution but use scaled weights as color/fill
    node_weights_forplot<-node_weights[order(node_weights$paths),]
    node_weights_forplot<-node_weights_forplot[!duplicated(node_weights_forplot$paths),]
    all(node_weights_forplot$paths==node_dat$total_nodes)
  # TRUE
    node_dat$weight<-node_weights_forplot$weight
    path_specific_nodes_hist_immunovar2<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=weight)) + geom_bar(stat="identity") +
      xlab("Path Depth") + ylab("Nodes (N)") +
      theme(legend.direction = 'horizontal', legend.position = 'bottom') +
      scale_fill_gradientn(name = 'Scaled Weight',trans="pseudo_log",colors = met.brewer("Egypt",n=15,type="continuous")) +
      theme_bw() #+ theme(legend.position="none")
    ggsave(paste0(outprefix,"_PathSpecificNodes_hist_weighted.png"), path_specific_nodes_hist_immunovar2,device = "png",width=7,height=4,units='in',dpi=1000)
    dat_corpus<-dat
    dat_corpus$V3<-gsub(","," ",dat_corpus$V3);
    dat_corpus$gene<-gsub("^.*#1#","",dat_corpus$V2);
    dat_corpus$gene_family<-substr(gsub("^.*#1#","",dat_corpus$V2),1,4);
    dat_corpus$gene_family<-gsub("\\*.*","",dat_corpus$gene_family);
    dat_corpus$genotype<-gsub("\\*.*$","",dat_corpus$gene)
    dat_corpus$genotype<-gsub("_.*$","",dat_corpus$genotype) # account for IgV2/OGRDB format
    my_corpus<-corpus(dat_corpus,docid_field = "V2",text_field = "V3",meta=dat_corpus$V2)
    # Create docvar with ID
    docvars(my_corpus, "id_numeric") <- 1:ndoc(my_corpus)
    #
    print(table(docvars(my_corpus)$gene_family))
    pathref_complete<-as.list(names(pathref))
    names(pathref_complete)<-names(pathref)
    for(i in 1:length(pathref)){
      ref<-names(pathref)[i]
      tmp.weights<-node_weights[node_weights$node %in% pathref[[i]],]$weight
      if(all(node_weights[which(node_weights$node %in% pathref[[i]]),]$node %in% pathref[[i]])){ #must be TRUE - mapping correct
        pathref_complete[[i]]<-list("nodes"=node_weights[which(node_weights$node %in% pathref[[i]]),]$node,"weights"=node_weights[which(node_weights$node %in% pathref[[i]]),]$weight)
      } else {
      print("Error!")
      }
    }
  #--
    fwrite(data.frame(node_weights$node),file=paste0(outprefix,"_Nodes.txt"),sep="\t",col.names=F,row.names=F,quote=F)
    save(my_corpus,node_weights,pathref_complete,file=paste0(outprefix,"_PathSpecificNodes_wWeights.R"))
  } else {
    if(length(setdiff(genes,c("IMGT","IGv2","OGRDB")))>0){
      additional_loci<-setdiff(genes,c("IMGT","IGv2","OGRDB"))
      print(paste0(additional_loci," in addition to IG alleles in this graph"))
      additional_loci<-c(additional_loci,"IMGT|IGv2|OGRDB")
      node_weights_combined<-as.data.frame(NULL)
  # these plots are misleading for IG genes - multiple possible paths for each gene, should instead do this at the genotype level rather than the path level
      for(g in seq_along(additional_loci)){
        locus<-paste0(additional_loci[g],"#1#")
        print(paste0("Parsing graph for ",gsub("#.*$","",locus)," alleles/paths"))
        total_nodes<-as.vector(NULL)
        for(i in 1:length(pathref[grep(locus,names(pathref))])){
          # Want unique nodes traversed by each path - can traverse same node multiple times
          total_nodes<-c(total_nodes,as.character(unique(pathref[grep(locus,names(pathref))][[i]])))
          }
            total_nodes<-table(total_nodes) # repeat for more intuitive plot
            node_dat<-data.frame(table(total_nodes))
            node_dat$nodes<-as.numeric(levels(node_dat$total_nodes)) # factor, so use the levels
            path_specific_nodes_hist_immunovar<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=nodes)) + geom_bar(stat="identity") + xlab("Number of Paths/Node") + ylab("N") +
              scale_fill_gradient(name = '',low="darkblue",high="red") + theme(legend.direction = 'horizontal', legend.position = 'bottom') +
              theme_bw() + theme(legend.position="none") + ggtitle(paste0(additional_loci[g]," Alleles: ",outprefix," graph"))
            ggsave(paste0(outprefix,"_",gsub("\\|","-",additional_loci[g]),"_PathSpecificNodes_hist.svg"), path_specific_nodes_hist_immunovar,device = "svg",width=4,height=5,units='in',dpi=1000) 
            # - make lower values weighted more + normalize
            total_nodes_wt<-max(total_nodes)/total_nodes # max num. of shared nodes=17 --> downweight to 1
            total_nodes_wt = (total_nodes_wt-min(total_nodes_wt))/(max(total_nodes_wt)-min(total_nodes_wt))
            node_weights<-data.frame(total_nodes_wt)
            colnames(node_weights)<-c("node","weight")
            all(names(total_nodes)==node_weights$node) # TRUE
            node_weights$paths<-total_nodes
            node_weights$unscaled_weight<-max(total_nodes)/node_weights$paths
            # ensure numeric values treated appropriately
            node_weights$node2<-as.numeric(as.character(node_weights$node))
            node_weights<-node_weights[order(node_weights$node2),]
            node_weights$node<-node_weights$node2
            node_weights$node2<-NULL
            # append to combined dataframe
            node_weights_combined<-rbind(node_weights_combined,node_weights)
            #fwrite(node_weights,file=paste0(outprefix,"_NodeWeightsPaths_Reference.csv"),quote=F,sep=",",col.names=T,row.names=F)
            ##############
            # remake distribution but use scaled weights as color/fill
            node_weights_forplot<-node_weights[order(node_weights$paths),]
            node_weights_forplot<-node_weights_forplot[!duplicated(node_weights_forplot$paths),]
            all(node_weights_forplot$paths==node_dat$nodes)
            # TRUE
            node_dat$weight<-node_weights_forplot$weight
            #
            path_specific_nodes_hist_immunovar2<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=weight)) + geom_bar(stat="identity") +
              xlab("Path Depth") + ylab("Nodes (N)") + theme(legend.direction = 'horizontal', legend.position = 'bottom') +
              scale_fill_gradientn(name = 'Scaled Weight',trans="pseudo_log",colors = met.brewer("Egypt",n=15,type="continuous")) +
              theme_bw() + ggtitle(paste0(additional_loci[g]," Alleles: ",outprefix," graph")) #+ theme(legend.position="none")
            ggsave(paste0(outprefix,"_",gsub("\\|","-",additional_loci[g]),"_PathSpecificNodes_hist_weighted.svg"), path_specific_nodes_hist_immunovar2, device = "svg",width=7,height=4,units='in',dpi=1000) 
  #          png(paste0(outprefix,"_",additional_loci[g],"_PathSpecificNodes_hist_weighted.png"),type="cairo",width = 7, height = 4, units= 'in', res=1000)
  #          path_specific_nodes_hist_immunovar2
  #          dev.off()
      }
    node_weights<-node_weights_combined
    } else {
      print("Only IG alleles in this graph")
      total_nodes<-as.vector(NULL)
      for(i in 1:length(pathref)){
        # Want unique nodes traversed by each path - can traverse same node multiple times
        total_nodes<-c(total_nodes,as.character(unique(pathref[[i]])))
      }
      total_nodes<-table(total_nodes) # repeat for more intuitive plot
      node_dat<-data.frame(table(total_nodes))
      node_dat$nodes<-as.numeric(node_dat$total_nodes)
      path_specific_nodes_hist_immunovar<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=nodes)) + geom_bar(stat="identity") + xlab("Number of Paths/Node") + ylab("N") +
        scale_fill_gradient(name = '',low="darkblue",high="red") + theme(legend.direction = 'horizontal', legend.position = 'bottom') +
        theme_bw() + theme(legend.position="none")
      ggsave(paste0(outprefix,"_IMGT_PathSpecificNodes_hist.svg"), path_specific_nodes_hist_immunovar,device = "svg",width=4,height=5,units='in',dpi=1000)
    ######################################################################################################
    # - make lower values weighted more + normalize
      total_nodes_wt<-max(total_nodes)/total_nodes # max num. of shared nodes=17 --> downweight to 1
      total_nodes_wt = (total_nodes_wt-min(total_nodes_wt))/(max(total_nodes_wt)-min(total_nodes_wt))
      node_weights<-data.frame(total_nodes_wt)
      colnames(node_weights)<-c("node","weight")
      all(names(total_nodes)==node_weights$node) # TRUE
      node_weights$paths<-total_nodes
      node_weights$unscaled_weight<-max(total_nodes)/node_weights$paths
      # ensure numeric values treated appropriately
      node_weights$node2<-as.numeric(as.character(node_weights$node))
      node_weights<-node_weights[order(node_weights$node2),]
      node_weights$node<-node_weights$node2
      node_weights$node2<-NULL
      #fwrite(node_weights,file=paste0(outprefix,"_NodeWeightsPaths_Reference.csv"),quote=F,sep=",",col.names=T,row.names=F)
      ##############
    # remake distribution but use scaled weights as color/fill
      node_weights_forplot<-node_weights[order(node_weights$paths),]
      node_weights_forplot<-node_weights_forplot[!duplicated(node_weights_forplot$paths),]
      all(node_weights_forplot$paths==node_dat$total_nodes)
    # TRUE
      node_dat$weight<-node_weights_forplot$weight
    #
      path_specific_nodes_hist_immunovar2<-ggplot(data=node_dat,aes(x=nodes,y=Freq,fill=weight)) + geom_bar(stat="identity") +
        xlab("Path Depth") + ylab("Nodes (N)") +
        theme(legend.direction = 'horizontal', legend.position = 'bottom') +
        scale_fill_gradientn(name = 'Scaled Weight',trans="pseudo_log",colors = met.brewer("Egypt",n=15,type="continuous")) +
        theme_bw() #+ theme(legend.position="none")
    #
      png(paste0(outprefix,"_PathSpecificNodes_hist_weighted.png"),type="cairo",width = 7, height = 4, units= 'in', res=1000)
      path_specific_nodes_hist_immunovar2
      dev.off()
    }
    #--
    dat_corpus<-dat
    dat_corpus$V3<-gsub(","," ",dat_corpus$V3);
    dat_corpus$gene<-gsub("^.*#1#","",dat_corpus$V2);
    dat_corpus$gene_family<-substr(gsub("^.*#1#","",dat_corpus$V2),1,4);
    dat_corpus$gene_family<-gsub("\\*.*","",dat_corpus$gene_family);
    dat_corpus$genotype<-gsub("\\*.*$","",dat_corpus$gene)
    dat_corpus$genotype<-gsub("_.*$","",dat_corpus$genotype) # account for IgV2/OGRDB format
    my_corpus<-corpus(dat_corpus,docid_field = "V2",text_field = "V3",meta=dat_corpus$V2)
    # Create docvar with ID
    docvars(my_corpus, "id_numeric") <- 1:ndoc(my_corpus)
    #
    table(docvars(my_corpus)$gene_family)
    #IGH IGK IGL TRA TRB TRD TRG ...
    #515 129 137 185 167  14  31 ...
    #
  #  IGH_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="IGH"))))
  #  IGK_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="IGK"))))
  #  IGL_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="IGL"))))
  #  TRA_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="TRA"))))
  #  TRB_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="TRB"))))
  #  TRD_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="TRD"))))
  #  TRG_nodes<-unique(unname(unlist(tokens_subset(tokens(my_corpus),gene_family=="TRG"))))
    #
  ########################################################################
  # create new list containing reference, nodes, and weights
    pathref_complete<-as.list(names(pathref))
    names(pathref_complete)<-names(pathref)
    for(i in 1:length(pathref)){
      ref<-names(pathref)[i]
      tmp.weights<-node_weights[node_weights$node %in% pathref[[i]],]$weight
      if(all(node_weights[which(node_weights$node %in% pathref[[i]]),]$node %in% pathref[[i]])){ #must be TRUE - mapping correct
        pathref_complete[[i]]<-list("nodes"=node_weights[which(node_weights$node %in% pathref[[i]]),]$node,"weights"=node_weights[which(node_weights$node %in% pathref[[i]]),]$weight)
    } else {
      print("Error!")
      }
    }
  #--
    fwrite(data.frame(node_weights$node),file=paste0(outprefix,"_Nodes.txt"),sep="\t",col.names=F,row.names=F,quote=F)
    #save(my_corpus,node_weights,pathref_complete,coocCounts,IGH_nodes,IGK_nodes,IGL_nodes,TRA_nodes,TRB_nodes,TRD_nodes,TRG_nodes,file=paste0(outprefix,"_PathSpecificNodes_wWeights.R"))
    #save(my_corpus,node_weights,pathref_complete,IGH_nodes,IGK_nodes,IGL_nodes,TRA_nodes,TRB_nodes,TRD_nodes,TRG_nodes,file=paste0(outprefix,"_PathSpecificNodes_wWeights.R"))
    save(my_corpus,node_weights,pathref_complete,file=paste0(outprefix,"_PathSpecificNodes_wWeights.R"))
  }
}


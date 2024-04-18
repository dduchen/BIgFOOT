#!/bin/bash

cd ${datadir}
gene=${each%.nodes.txt}
gene=${gene%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
#
loci=$(vg paths -Lv ${graph_base}.xg | grep "#1#${gene}\*" | cut -f1 -d"#" | sort | uniq | grep "IMGT\|HLA\|KIR")
echo "$gene --> $loci"
#
if [[ "$loci" =~ ^(HLA)$ ]]; then
    echo "HLA inference"
    outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference/HLA
    mkdir -p ${outdir}/seqwish_${sample_id}.${graph}
elif [[ "$loci" =~ ^(KIR)$ ]]; then
    echo "KIR inference"
    outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference/KIR
    mkdir -p ${outdir}/seqwish_${sample_id}.${graph}
elif [[ "$loci" =~ ^(IMGT)$ ]]; then
    echo "IG/TR inference"
    outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference
    mkdir -p ${outdir}/seqwish_${sample_id}.${graph}
    # ASC cluster - if more than 1 gene included in ASC cluster - use clustering-based gene ID for graph construction:
    grep "${gene}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv > ${outdir}/potential_asc_for_${gene}
    if [ -s ${outdir}/potential_asc_for_${gene} ]; then
    # if n_genes > 1 then use asc cluster
        asc_cluster=$(cut -f4 ${outdir}/potential_asc_for_${gene} | sed s/'\*.*'//g | sort | uniq)
        if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
            echo "Multiple V genes in this cluster";
            for asc_multi in ${asc_cluster[@]}; do echo $asc_multi;
                echo "Concatentating nodes across cluster members";
                cat ${genotyping_nodes_dir}/ig_asc/${asc_multi}.immune_subset.nodes.txt >> ${outdir}/${gene}_asc_nodes.txt
            done
            sort ${outdir}/${gene}_asc_nodes.txt | uniq > ${genotyping_nodes_dir}/ig_asc/${gene}_asc_nodes.txt
            each=ig_asc/${gene}_asc_nodes.txt
        else
            grep ${asc_cluster} ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f1 > ${outdir}/potential_asc_for_${gene}
            grep -v ${gene}"\*" ${outdir}/potential_asc_for_${gene} | grep -v "${gene}_"  > ${outdir}/asc_additional_${gene}_fams.txt
            if [ -s ${outdir}/asc_additional_${gene}_fams.txt ]; then
                echo "Complex locus - using ASC-based allele";
                each=ig_asc/${asc_cluster}.immune_subset.nodes.txt
            else
                echo "Using core gene+haplotype graph for ${gene}";
            fi
        fi
    fi
else 
    echo "Unknown locus for ${gene}"
fi
# check if graph exists in subdir - if it does use it rather than reconstruct things
# check if we've completed analysis for this gene already
if [ -s ${outdir}/${sample_id}_${gene}_files.txt ]; then
    echo "Analysis completed for ${gene} - did you restart or fail to cleanup? (File=${outdir}/${sample_id}_${gene}_files.txt)"
else
    if [[ "$loci" =~ ^(IMGT)$ ]]; then
        echo "IMGT: Immunoglobulin/T-cell receptor gene inference"
        if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gcsa ]; then
            echo "Using stored locus specific graph/indexes";
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.genotyping.immune_subset.vg ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gbwt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
        else
            echo "subsetting graph and alignment to region surrounding locus of interest: $gene"
            vg find -x ${graph_base}.xg -c 150 -L -N ${genotyping_nodes_dir}/${each} > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
            echo "Extracting local haplotypes reads --> ${sample_id}.${graph}.${gene}.haplotypes.gfa"
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "grch\|chm" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt | 
            vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg -F > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
            seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "IGL\|hirh_\|IGK\|MHC\|IMG\|IGv\|OGR" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt
            # vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "IGL\|hirh_\|IGK\|MHC\|${gene_actual}" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt | 
            vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
            # haplotype graph:
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -F -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            # removing duplicates - ensuring 'gene' matches are at top of fasta
            seqkit grep -r -p ${gene}"\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.match.fasta
            seqkit grep -r -v -p ${gene}"\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.unmatch.fasta
            cat ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.match.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.unmatch.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            seqkit grep -r -p "IMG|IGv|OGR" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >  ${outdir}/${gene}.alleles.fasta
            # set min length of haplotypes = 100, max length should scale with allele length --> then downsample haplotypes
            seqkit stats ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.stats
            gene_min_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f6 -d' ' | sed s/","//g)
            #gene_max_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f8 -d' ')
            ##############################
            # use local haplotypes #
            ##############################
            seqkit grep -v -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.haps.fasta
            if [ -s ${outdir}/${gene}.haps.fasta ]; then
                echo "Local haplotypes found for: ${gene_actual}"
                seqkit seq --min-len $(bc -l <<< "scale=2;${gene_min_len}*0.9"| awk '{printf("%d\n",$1 + 0.5)}') --max-len 15000  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                # retain haplotypes with perfect match to one of our alleles
                mkdir -p $outdir/${gene}_haps
                if [[ ${gene} == *["VJ"]* ]]; then
                    echo "Exact allele:haplotype matching"
                    minimap2 -x sr -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                else
                    echo "Fuzzy allele:haplotype matching"
                    minimap2 -x sr -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | egrep "NM:i:0|NM:i:[0-9]$" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                fi
                echo "$(cut -f1 ${outdir}/${gene}_haps/haps.matching.txt | wc -l) out of $(grep ">" ${outdir}/${gene}.haps.fasta | cut -f1 | wc -l) haplotypes with a $gene (or ASC cluster member) allele mapping"
                if [ -s ${outdir}/${gene}_haps/haps.matching.txt ]; then
                    seqkit grep -r -f $outdir/${gene}_haps/haps.matching.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    seqkit split --quiet -i $outdir/${gene}_haps/${gene}.haps.fasta -f
                    mkdir -p $outdir/${gene}_haps_final
                    ${tools_dir}/Assembly-dereplicator/dereplicator.py --distance 0.000001 --count 25 $outdir/${gene}_haps/${gene}.haps.fasta.split/ $outdir/${gene}_haps_final
                    cat $outdir/${gene}_haps_final/*fasta > ${outdir}/${gene}.haps.fasta
                    seqkit stats ${outdir}/${gene}.haps.fasta
                    rm -rf $outdir/${gene}_haps_final ; rm -rf $outdir/${gene}_haps
                    cat ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                    seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                else
                    echo "No haplotypes containing the alleles"
                    rm -rf $outdir/${gene}_haps
                fi
            else
                echo "No local haplotypes found for ${gene_actual}"
                rm -f ${outdir}/${gene}.haps.fasta
            fi
            # remove duplicate sequences
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            # more than one chromosomal locus represented?
            if [ $(grep "chm13#chr" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt | cut -f1 -d':' | wc -l) -gt 1 ]; then
                echo "more than 1 chromosome - forcing multi-component graph construction"                
                seqkit grep -r -p "chm13" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | minimap2 -x asm20 -c - ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt
                for chr_loc in $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                    chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                    grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                    seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf -g ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                done
                echo "combine subgraphs"
                ls ${outdir}/${sample_id}.${graph}.${gene}.chr*gfa > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og 
                odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
                #delete temp files
                rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt; rm ${outdir}/${sample_id}.${graph}.${gene}.chr*
            else
                echo "single chromosomal locus - standard graph construction"
                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf
                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            fi 
            #odgi explode -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -b 1 -s P -O -P -p "${sample_id}.${gene}.exp" -g
            #mv ${sample_id}.${gene}.exp*gfa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg;
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -P --pass-paths
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz --gbz-format -P --pass-paths;
            vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
            vg index -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
            vg index -j ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz
            # save in genotyping directory so this doesnt have to be repeated!
            cp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.genotyping.immune_subset.vg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.xg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gcsa
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gcsa.lcp
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.gbwt
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes.pg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haplotypes_ref.fasta
        fi
        # parse alignments
        vg find -x ${graph_base}.xg -l ${gam_file%.gam}.sorted.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam
        vg view -a ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam | sed s/'.*"fragment_next": {"name": '//g | sed s/'.*"fragment_prev": {"name": '//g | sed s/', "name": ".*'//g > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt
        if [[ $(wc -l <${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt) -ge 1 ]]; then
            echo "Re-aligning reads to locus-specific haplotype graph"
            vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            # copy mapped reads -- use this for augmentation?
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
            vg filter -r 0 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            echo "1) Performing inference on haplotype graph labled only with alleles of interest";
            # if variable gene -- limit to OGRDB validated set of alleles?
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "IMGT\|IGv2\|OGRDB" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
            if [ "${valid_alleles}" = true ]; then
                echo "Where possible - retaining path labels only for OGRDB validated set of ${gene_actual} alleles";
                if [[ ${gene} == *"IGHV"* ]]; then
                    echo "IGHV";
                    ogrdb_refs="${bigfoot_dir}/../custom_beds/ogrdb_IGH_*fasta"
                    grep ${gene}\\* $ogrdb_refs | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    fi
                elif [[ ${gene} == *"IGKV"* ]]; then
                    echo "IGKV";
                    ogrdb_refs=${bigfoot_dir}/../custom_beds/ogrdb_IGK_*fasta
                    grep ${gene}\\* $ogrdb_refs | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    fi
                elif [[ ${gene} == *"IGLV"* ]]; then
                    echo "IGLV";
                    ogrdb_refs=${bigfoot_dir}/../custom_beds/ogrdb_IGL_*fasta
                    grep ${gene}\\* $ogrdb_refs | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    fi
                else
                    echo "Retaining path labels only for ${gene_actual}";
                    grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                fi
            else
                echo "Retaining path labels only for ${gene_actual} If you want to limit variable gene inference to OGRDB set - set 'valid_alleles=true'"
                grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
            fi
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
            vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
            # gene-specific read depth for minimum strain-level coverage
            depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
            # min_strain_depth=$(bc -l <<< "scale=3;${depth_locus}/$gene_min_len")
            min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.05"| awk '{printf("%d\n",$1 + 0.5)}')
            if (( $(awk 'BEGIN {print ("'"$min_strain_depth"'" < 0.1) ? "1" : "0"}') )); then
                min_strain_depth=0.1
            fi
            echo "Minimum strain depth required: $min_strain_depth"
            cd ${outdir}
            # allele inference:
            mkdir -p ${outdir}/${gene}_allele_inference
            python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
            for opt in {2..2}; do #relative difference performs best on average
                echo "basing inference on haplotypes + alleles embedded in graph"
                cd ${outdir}/${gene}_allele_inference
                python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c ${min_strain_depth} --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
                # python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --careful --min_depth 0 --trim 0 -m 0 -c ${min_strain_depth} --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
                if [[ "$opt" =~ 1 ]]; then
                    echo "optimization_approach = absolute difference"
                    mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.abs.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.abs.haps.final.fasta ; 
                    mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.abs.genome_graph.gfa
                    rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
                    Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.abs.contigs.fasta
                    echo "Allele-level abundance estimation completed for ${gene} ::"
                    grep ">" ${outdir}/${sample_id}.${graph}.${gene}.abs.haps.final.annot.fasta
                    cd ${outdir}
                elif [[ "$opt" =~ 2 ]]; then
                    echo "optimization_approach = relative difference";
                    if [ -s haps.final.fasta ]; then
                        echo "$(grep ">" haps.final.fasta | wc -l) alleles > threshold";
                        echo "Matching inferred alleles to allele set + constructing allele-specific ${gene} graph"
                    else
                        echo "No strains with estimated allele coverage >= ${min_strain_depth} - using basal threshold (0.01) to rescue inference";
                        python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c 0.1 --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
                    fi
                    mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.fasta ;
                    mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.rel.genome_graph.gfa;
                    rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
                    if [ -s ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.fasta ]; then
                        Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.rel.contigs.fasta
                        echo "Allele-level abundance estimation completed for ${gene} ::"
                        grep ">" ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta
                    else
                        echo "No final haplotypes > depth threshold"
                    fi
                    cd ${outdir}
                elif [[ "$opt" =~ 3 ]]; then
                    echo "optimization_approach = relative absolute differences sqrt approach"
                    mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.relabs.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.relabs.haps.final.fasta ; 
                    mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.relabs.genome_graph.gfa
                    rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
                    Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.relabs.contigs.fasta
                    echo "Allele-level abundance estimation completed for ${gene} ::"
                    grep ">" ${outdir}/${sample_id}.${graph}.${gene}.relabs.haps.final.annot.fasta
                    cd ${outdir}
                else
                    echo "Select an acceptable optimization approach (1: absolute diff, 2: relative diff, 3: sqrt(relative abs diff)"
                    cd ${outdir}
                fi
            done
            rm -rf ${outdir}/${gene}_allele_inference
            echo "2) Augmenting annotated post-flow inference graph with reads for association testing";
            sed s/' '/_/g ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta > ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta
            if [ -s "${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta" ]; then    
                seqkit grep -r -p "chm" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
                if [ $(grep "chm13#chr" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | cut -f1 -d':' | wc -l) -gt 1 ]; then
                    echo "Likely orphon involvement - ensuring we use best matching reference locus"
                    seqkit grep -r -p "chm13" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | minimap2 -x asm20 -c - ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt
                    for chr_loc in $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                        chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                        grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                        seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta
                        seqkit grep -p $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | cat - ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa
                    done
                    if [ $(ls ${outdir}/${sample_id}.${graph}.${gene}.chr*.genome_graph_ref.gfa | wc -l) -gt 1 ]; then
                        echo "combining relevant subgraphs"
                        ls ${outdir}/${sample_id}.${graph}.${gene}.chr*.genome_graph_ref.gfa > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                        odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.og
                        odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.og -g > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                        #delete temp files
                        rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.og;rm ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta; rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt; rm ${outdir}/${sample_id}.${graph}.${gene}.chr*
                    else
                        mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                        rm ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta;
                    fi
                else
                    cat ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                    echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local CHM13 reference sequence"
                    sed s/' '/_/g -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                fi
                vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg
                vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg;
                vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -P --pass-paths
                vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbz --gbz-format -P --pass-paths;
                vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
                vg index -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
            #   map locus-associated reads to augmented/annotated graph
                vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
                vg filter -r 0.95 -P -s 1 -q 60 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
                vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth;
                depth_aug=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth)
                aug_depth=$(bc -l <<< "scale=2;${depth_aug}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
                if [ "${aug_depth}" -gt 3 ]; then
                    augment_cov=${aug_depth}
                else
                    augment_cov=3
                fi
                echo "Minimum coverage to add breakpoint: ${augment_cov} (3 <--> 10% of strain depth)"
            #
                vg augment -m ${augment_cov} -q 5 -Q 60 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg;
                vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
                vg mod -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
            else
                echo "Insuffient coverage to infer alleles for $gene - not trying to embed variation"
            fi
            # remove files we dont need anymore
            ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v "${gene}.genome_graph_ref.augmented.gfa\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|node_abundance\|depth" > ${outdir}/${sample_id}_${gene}_files.txt
            #ls ${outdir}/${gene}\.*  | grep "haps.fasta\|alleles" >> ${outdir}/${sample_id}_${gene}_files.txt
            ls ${outdir}/*asc_*${gene}* >> ${outdir}/${sample_id}_${gene}_files.txt
            xargs rm < ${outdir}/${sample_id}_${gene}_files.txt
        else 
            echo "No reads aligning for ${gene}"
            ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v  "${gene}.genome_graph_ref.augmented.gfa\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|node_abundance\|depth" > ${outdir}/${sample_id}_${gene}_files.txt
            xargs rm < ${outdir}/${sample_id}_${gene}_files.txt
        fi
    elif [[ "$loci" =~ ^(HLA)$ ]]; then 
        downsampled=FALSE
        if [[ "$downsampled" =~ ^(TRUE)$ ]]; then
            echo "basing inference on haplotypes + alleles embedded in graph, for inference based on full set of alleles set downsampled=FALSE (this might be much slower for certain HLA genes)"
            # haplotype graph:
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -F -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            # downsample to 2 typing codes --> perform inference (helps with DRB1 / highly polymorphic genes):
            vg paths -Lx ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${sample_id}.${graph}.paths_to_keep.txt 
            while IFS= read -r line; do
                modified_line=$(echo "$line" | awk -F':' 'BEGIN{OFS=":"} {split($2, a, ":"); $2=a[1] ":" a[2];}1')
                modified_line=$(echo $modified_line | sed s/'::.*'// | sed s/'#0'// | sed s/':$'//)
                echo -e "${line}\t${modified_line}"
            done < ${sample_id}.${graph}.paths_to_keep.txt | sort -k2,2 | awk -F'\t' '!seen[$2]++' > ${sample_id}.${graph}.paths_to_keep_reduced.txt
            cut -f 2 ${sample_id}.${graph}.paths_to_keep_reduced.txt | awk -F':' '{ if (length($NF) >= 3) $NF = substr($NF, 1, 2) } 1' OFS=':' | paste -d'\t' ${sample_id}.${graph}.paths_to_keep_reduced.txt - | sort -k3,3 | awk -F'\t' '!seen[$3]++'> ${sample_id}.${graph}.paths_to_keep_reduced.txt.tmp && mv ${sample_id}.${graph}.paths_to_keep_reduced.txt.tmp ${sample_id}.${graph}.paths_to_keep_reduced.txt
            cut -f1 ${sample_id}.${graph}.paths_to_keep_reduced.txt > ${sample_id}.${graph}.paths_to_keep_reduced.txt.tmp && mv ${sample_id}.${graph}.paths_to_keep_reduced.txt.tmp ${sample_id}.${graph}.paths_to_keep_reduced.txt
            vg paths -v ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -F -p ${sample_id}.${graph}.paths_to_keep_reduced.txt > ${outdir}/${gene}.alleles.fasta
            # set min length of haplotypes = 100, max length should scale with allele length --> then downsample haplotypes
                seqkit stats ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.stats
                gene_min_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f6 -d' ' | sed s/,//g)
                wfmash_param=$(bc -l <<< "scale=2;${gene_min_len}/10" | awk '{printf("%d\n",$1 + 0.5)}')
                if [ "${wfmash_param}" -lt 100 ]; then
                    wfmash_param=100
                fi
            ####################################
            # use local downsampled haplotypes #
            ####################################
            seqkit grep -r -p "MHC" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.haps.fasta
            if [ -s ${outdir}/${gene}.haps.fasta ]; then
                echo "Local haplotypes found for: ${gene_actual}"
                seqkit seq --min-len ${gene_min_len} --max-len 1000000  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                # retain haplotypes with close match to one of our alleles
                mkdir -p $outdir/${gene}_haps
                # allow for an edit distance less than 10bp
                # time minimap2 -x asm20 -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | cut -f1,13 | egrep "NM:i:0|NM:i:[0-9]$" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                time blend -x asm20 -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | cut -f1,13 | egrep "NM:i:0|NM:i:[0-9]$" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                echo "$(cut -f1 ${outdir}/${gene}_haps/haps.matching.txt | wc -l) out of $(grep ">" ${outdir}/${gene}.haps.fasta | cut -f1 | wc -l) haplotypes with a close match (edit distance<=10) to a $gene allele mapping"
                seqkit grep -r -f $outdir/${gene}_haps/haps.matching.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                seqkit split --quiet -i $outdir/${gene}_haps/${gene}.haps.fasta -f
                mkdir -p $outdir/${gene}_haps_final
                ${tools_dir}/Assembly-dereplicator/dereplicator.py --distance 0.0000001 --count 25 $outdir/${gene}_haps/${gene}.haps.fasta.split/ $outdir/${gene}_haps_final
                cat $outdir/${gene}_haps_final/*fasta > ${outdir}/${gene}.haps.fasta
                seqkit stats -a -G 'N' ${outdir}/${gene}.haps.fasta
                rm -rf $outdir/${gene}_haps_final ; rm -rf $outdir/${gene}_haps
                cat ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            else
                echo "No local haplotypes found for ${gene_actual}"
                rm -f ${outdir}/${gene}.haps.fasta
            fi
            # remove duplicate sequences
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            samtools faidx ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            blend -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf
            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg;
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -P --pass-paths
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz --gbz-format -P --pass-paths;
            vg index -j ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz
            vg minimizer -p -d ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.min ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz
            rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_paths.gbwt
            vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -k 45 -m ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
            vg index -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
            # use VG map if lots of reads
            reads_aln=$(vg stats -a ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam | grep "aligned:" | cut -f2 -d":")
            echo "We're working with ${reads_aln} reads for ${gene}"
            echo "Re-aligning reads to locus-specific haplotype graph"
            if [ "${reads_aln}" -gt 250000 ]; then
                time vg giraffe -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -Z ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz -m ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.min -H ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                vg filter -r 0.9 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam
                time vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.filt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            else
                echo "VG map alignment"
                time vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            fi
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
            vg filter -r 0.8 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            # 1) should perform inference on haplotype graph with non-allele paths removed?
            echo "Limiting inference to alleles of interest"
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "HLA" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
            vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
            depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
            min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.05"| awk '{printf("%d\n",$1 + 0.5)}')
            if (( $(awk 'BEGIN {print ("'"$min_strain_depth"'" < 0.1) ? "1" : "0"}') )); then
                min_strain_depth=0.1
            fi
            cd ${outdir}
            python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
            time python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --min_depth 0 --trim 0 --greedy_mode all --ilp --max_strains 2 -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 --wfmash_param ${wfmash_param} ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
            mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.haps.final.fasta ; mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph.gfa
            rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
            Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
            echo "Allele-level abundance estimation completed for ${gene} ::"
            grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta
        else
            echo "basing inference on full set of alleles for ${gene}, for inference based on embedded haplotypes+alleles set downsampled=TRUE"
            # 2) - extract reference backbone to embed variation along with called alleles
            echo "MHC/HLA locus :: Using pre-constructed IPD-based HLA graphs for inference"
            # re-align to haplotype graph + embed variation with adequate support
            time vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -g ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gcsa -1 ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            vg filter -r 0 -P -s 1 -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.tmp.gam && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.tmp.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
            min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
            if (( $(awk 'BEGIN {print ("'"$min_strain_depth"'" < 0.1) ? "1" : "0"}') )); then
                min_strain_depth=0.1
            fi        
            echo "Minimum strain depth required: ${min_strain_depth}"
            vg convert -fW ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
            vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
            python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
            gene_min_len=$(vg paths -Ex ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg | grep "HLA" | cut -f2 | sort | uniq | head -1)
            wfmash_param=$(bc -l <<< "scale=2;${gene_min_len}/10" | awk '{printf("%d\n",$1 + 0.5)}')
            if [ "${wfmash_param}" -lt 100 ]; then
                wfmash_param=100
            fi
            time python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --min_depth 0 --trim 0 --greedy_mode all --ilp --max_strains 2 -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 --wfmash_param ${wfmash_param} ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
            mv trimmed_contigs_long.fasta ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
            mv haps_long.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.haps.final.fasta
            mv genome_graph_long.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph.gfa
            rm genome_graph_long.gt; rm haps_long.fasta; rm overlaps_long.minimap2.paf; rm trimmed_contigs_long.paths ; rm trimmed_contigs_long.gfa; rm trimmed_contigs_long.fasta.fai
            Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
            echo "Allele-level abundance estimation completed for ${gene} ::"
            grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta
            # 2) augment annotated post-flow inference graph with reads for association testing
            # embed along with reference for variant calling
            cat ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
            echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local CHM13 reference sequence"
            sed s/' '/_/g -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
            blend -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
            vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg
            vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
            vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg;
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -P --pass-paths
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbz --gbz-format -P --pass-paths;
            vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
            vg index -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
            #   map locus-associated reads to augmented/annotated graph
            vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
            vg filter -r 0.95 -P -s 1 -q 60 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth;
            depth_aug=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth)
            aug_depth=$(bc -l <<< "scale=2;${depth_aug}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
            if [ "${aug_depth}" -gt 3 ]; then
                augment_cov=${aug_depth}
            else
                augment_cov=3
            fi
            echo "Minimum coverage to add breakpoint: ${augment_cov} (3 <--> strain depth @ gene)"
            vg augment -m ${augment_cov} -q 5 -Q 60 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg;
            vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
            vg mod -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
            # remove files we dont need anymore
            ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v  "${gene}.genome_graph_ref.augmented.gfa\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|node_abundance\|depth" > ${outdir}/${sample_id}_${gene}_files.txt
            ls ${outdir}/${gene}\.*  | grep "haps.fasta\|alleles" >> ${outdir}/${sample_id}_${gene}_files.txt
            xargs rm < ${outdir}/${sample_id}_${gene}_files.txt
            rm -rf ${outdir}/seqwish_${sample_id}.${graph}/*
        fi
    else
        echo "Unknown locus or incorrect script for ${gene}"
    fi
fi

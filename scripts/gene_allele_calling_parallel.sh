#!/bin/bash

cd ${datadir}
gene=${each%.nodes.txt}
gene=${gene%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
#
loci=$(vg paths -Lv ${graph_base}.xg | grep "#1#${gene}\*" | cut -f1 -d"#" | sort | uniq | grep "IMGT\|HLA\|KIR")
echo "$gene --> $loci"
#
unset asc_cluster #ensure no carryover of variables from previous gene analysis
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
            echo "Complex locus - using ASC-based allele -- ALSO - looks like there are multiple V genes in this cluster";
            complex_gene=true
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
                complex_gene=true
            else
                echo "Using core gene+haplotype graph for ${gene}";
            fi
        fi
    else
        echo "Using core gene+haplotype graph for ${gene}" #; <--- might need to do more here!
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
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gbwt ];then
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp
            fi
            # ensure P lines are at bottom of gfa file:
            grep "^P" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.plines
            grep "^P" -v ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.noplines
            cat ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.noplines ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.plines > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.noplines
            rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa.plines
        else
            echo "subsetting graph and alignment to region surrounding locus of interest: $gene"
            context_length=150
           # if [[ $gene == *"J-|D-"* ]]; then 
           #     echo "Using synthetic graph with haplotype flanks"
           #    #use similar approach for D genes - from expressed repertoire-related graph + relevant haplotype flanks
           #     context_length=300
           # fi
            vg find -x ${graph_base}.xg -c ${context_length} -L -N ${genotyping_nodes_dir}/${each} > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
            echo "Extracting local haplotypes reads --> ${sample_id}.${graph}.${gene}.haplotypes.gfa"
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "grch\|chm" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt
            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt ]; then
                echo "Local reference sequences obtained for ${gene}"
            else 
                echo "Expanding search space to find local reference sequence ${gene}"
                vg find -x ${graph_base}.xg -c 100 -N ${genotyping_nodes_dir}/${each} > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
                echo "Extracting local haplotypes reads --> ${sample_id}.${graph}.${gene}.haplotypes.gfa"
                vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "grch\|chm" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt
            fi
            vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt | \
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
            seqkit grep -r -p "IMG|IGv|OGR" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.alleles.fasta
            seqkit grep -r -p ${gene}"\*" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.exact.fasta
            seqkit grep -r -v -p ${gene}"\*" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.offtarget.fasta
            # set min length of haplotypes = 100, max length should scale with allele length --> then downsample haplotypes
            seqkit stats ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.stats
            gene_min_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f6 -d' ' | sed s/","//g)
            gene_max_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f8 -d' ')
            ##############################
            # use local haplotypes #
            ##############################
            seqkit grep -v -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.haps.fasta
            if [ -s ${outdir}/${gene}.haps.fasta ]; then
                echo "Local haplotypes found for: ${gene_actual}"
#                seqkit seq --min-len $(bc -l <<< "scale=2;${gene_min_len}*0.9"| awk '{printf("%d\n",$1 + 0.5)}') --max-len 15000  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                seqkit seq --min-len $(bc -l <<< "scale=2;${gene_min_len}*0.9"| awk '{printf("%d\n",$1 + 0.5)}') --max-len $(bc -l <<< "scale=2;${gene_max_len}*3"| awk '{printf("%d\n",$1 + 0.5)}')  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                # retain haplotypes with perfect match to one of our alleles
                mkdir -p $outdir/${gene}_haps
                if [[ ${gene} == *["VJ"]* ]]; then
                    echo "Exact allele:haplotype matching"
                    minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                    minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.exact.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.exact.txt
                    minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.offtarget.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.offtarget.txt
                    cp ${outdir}/${gene}_haps/haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt
                    cp ${outdir}/${gene}_haps/haps.matching.offtarget.txt ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget.txt
                else
                    echo "Fuzzy allele:haplotype matching (< 1% gap-compressed sequence divergence)"
#                    minimap2 -x sr -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | egrep "NM:i:([0-9]|10)\b" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                    minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | egrep "de:f:0\.00[0-9][0-9]?" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                fi
                echo "$(cut -f1 ${outdir}/${gene}_haps/haps.matching.txt | wc -l) out of $(grep ">" ${outdir}/${gene}.haps.fasta | cut -f1 | wc -l) haplotypes with a $gene (or ASC cluster member) allele mapping"
                echo "$(cut -f1 ${outdir}/${gene}_haps/haps.matching.exact.txt | wc -l) out of $(grep ">" ${outdir}/${gene}.haps.fasta | cut -f1 | wc -l) haplotypes with a $gene-specific allele mapping"
            else
                echo "No local haplotypes found for ${gene_actual}"
                rm -f ${outdir}/${gene}.haps.fasta
            fi
            # remove duplicate sequences
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            # checking number of chromosomes + primary allele<->haplotype alignments
            seqkit grep -r -p "chm13" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | minimap2 -x asm20 --secondary=no -c - ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt
            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt ]; then
                echo "Using CHM13 reference backbone"
            else 
                echo "Using GRCh38 reference backbone"
                seqkit grep -r -p "grch" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | minimap2 -x asm20 --secondary=no -c - ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt
            fi
            # account for some genes with no local haplotypes!
            if [ -s ${outdir}/${gene}.haps.fasta ]; then
                seqkit grep -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | minimap2 -x asm20 --secondary=no -c ${outdir}/${gene}.haps.fasta - | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt
            fi
            if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
                echo "Accounting for multi-ASC clustering"
                if [ -s ${outdir}/${gene}_haps/haps.matching.txt ]; then
                    if [ -s ${outdir}/${gene}_haps/haps.matching.exact.txt ]; then
                        echo "Haplotypes exist with exact matches to alleles"
                        seqkit grep -r -f ${outdir}/${gene}_haps/haps.matching.exact.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    else
                        echo "Haplotypes exist with exact matches to ASC-related genes only"
                        seqkit grep -r -f ${outdir}/${gene}_haps/haps.matching.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    fi
                    seqkit split --quiet -i ${outdir}/${gene}_haps/${gene}.haps.fasta -f
                    mkdir -p ${outdir}/${gene}_haps_final
                    ${tools_dir}/Assembly-dereplicator/dereplicator.py --distance 0.000001 --count $(grep ">" ${outdir}/${gene}.alleles.exact.fasta | wc -l) ${outdir}/${gene}_haps/${gene}.haps.fasta.split/ $outdir/${gene}_haps_final
                    cat ${outdir}/${gene}_haps_final/*fasta > ${outdir}/${gene}.haps.fasta
                    seqkit stats ${outdir}/${gene}.haps.fasta
                    rm -rf $outdir/${gene}_haps_final ; rm -rf $outdir/${gene}_haps
                    cat ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                    seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                else
                    echo "No haplotypes containing the alleles"
                    rm -rf $outdir/${gene}_haps
                fi
                #
                cut -f4 ${outdir}/potential_asc_for_${gene} | sed s/'\*.*'//g | sort | uniq > ${outdir}/potential_asc_for_${gene}_ascs.txt
                for clust in $(cat ${outdir}/potential_asc_for_${gene}_ascs.txt);do echo ${clust};
                    grep "${clust}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f1 | sed s/'*'/'\\*'/g > ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt
                    # check if >1 chromsome -- (cluster specific) is so, split again
                    if [ $(grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f2 | sort | uniq | wc -l) -gt 1 ]; then
                        echo "More than 1 chromosome represented in $clust -- additionally splitting by chr"
                        for chr_loc in $(grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f2 | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                            chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                            grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt - | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                            grep -f <(sed s/'\*'/'\\*'/g ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt) ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt | cut -f2 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                            # Alleles might align best to a speciifc haplotype - but need ensure the haplotype belongs in the chromosome-specific subgraph
                            if [ $(grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep ${chr_loc} | wc -l) -ge 1 ]; then
                                seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                            fi
                            minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf
                            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa
                        done
                        ls ${outdir}/${sample_id}.${graph}.${gene}.${clust}.*.gfa | grep -v "haplotypes" - > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt
                        odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og 
                        odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa
                        #delete temp files
                        rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt; rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.chr*
                    else
                        echo "Single locus for this cluster of ${gene}"
                        seqkit grep -r -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta
                        if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt ]; then
                            grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt | cut -f2 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotype_assignments.patterns.txt
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta
                        fi
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa
                    fi
                done
                echo "Constructing pan-cluster graph for ${gene}"
                ls ${outdir}/${sample_id}.${graph}.${gene}.*.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -s "#" -O -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og 
                odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
                # final inference step after filtering out spurious reads - consolidated gene-specific graph when alleles spread across clusters:
                # --> add something during inference: if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
                cat ${outdir}/${gene}.alleles.exact.fasta ${outdir}/${gene}.haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta
                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf
                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp;
                mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa
                # indexes for vg map
                vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg
                vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.vg
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.vg > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa
                vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg
                vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg;
                vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt -P --pass-paths
                vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pruned.vg
                vg index -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pruned.vg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gfa
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gbwt
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.xg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.pg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gcsa
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_locus.gcsa.lcp
                #delete temp files
                xargs rm < ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
            else
                echo "Single locus/cluster for ${gene}"
                if [ -s ${outdir}/${gene}_haps/haps.matching.txt ]; then
                    if [ -s ${outdir}/${gene}_haps/haps.matching.exact.txt ]; then
                        echo "Haplotypes exist with exact matches to alleles"
                        seqkit grep -r -f ${outdir}/${gene}_haps/haps.matching.exact.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    else
                        echo "Haplotypes exist with exact matches to ASC-related genes only"
                        seqkit grep -r -f ${outdir}/${gene}_haps/haps.matching.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    fi
                    seqkit split --quiet -i ${outdir}/${gene}_haps/${gene}.haps.fasta -f
                    mkdir -p ${outdir}/${gene}_haps_final
                    ${tools_dir}/Assembly-dereplicator/dereplicator.py --distance 0.000001 --count $(grep ">" ${outdir}/${gene}.alleles.exact.fasta | wc -l) $outdir/${gene}_haps/${gene}.haps.fasta.split/ $outdir/${gene}_haps_final
                    cat ${outdir}/${gene}_haps_final/*fasta > ${outdir}/${gene}.haps.fasta
                    seqkit stats ${outdir}/${gene}.haps.fasta
                    rm -rf $outdir/${gene}_haps_final ; rm -rf $outdir/${gene}_haps
                    cat ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                    seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                else
                    echo "No haplotypes containing the alleles"
                    rm -rf $outdir/${gene}_haps
                fi
               # this still checks for multiple chromosomes - if single chromosome, then this shouldnt change things
                for chr_loc in $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                    chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                    grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                    seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    if [ $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq | wc -l) -le 1 ]; then
                        echo "adding all alleles to subgraph"
                        seqkit grep -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    fi
                    seqkit rmdup -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    #grep -r ${chr_graph_tmp} ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep -v "${gene}\|IMGT\|OGRDB\|IGv2" - | cut -f1 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.haplotype_assignments.patterns.txt
                    #seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf -g ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                done
                echo "combining chr specific subgraphs for ${gene}"
                # ensuring we only combine relevant subgraphs - avoid patten matching from other gfa files potentially created during an incomplete run
                ls ${outdir}/${sample_id}.${graph}.${gene}.*gfa | grep -v "haplotypes\|final\|vgflow\|genome_graph" > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -s "#" -O -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og 
                odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
                #delete temp files
                xargs rm < ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt;
            fi
            # alternative haplotype filtering:
#            sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
#            if [ -s ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt ]; then
#                Rscript ${bigfoot_dir}/prune_locus_graph_haplotypes.R ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
#                # could use this to create cleaner haplotype graph - then augment to subgraphs specific to relevant haplotypes
#                vg paths -F ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p
#                seqkit grep -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta |  seqkit grep -r -p "${gene}" > ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta;
#                seqkit grep -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta |  seqkit grep -r -v -p "${gene}" > ${outdir}/${sample_id}.${graph}.${gene}.off_target_alleles.fasta;
#                seqkit grep -r -f <(cat <(cut -f1 ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact_ids.txt) <(cut -f1 ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget_ids.txt) | sort | uniq) ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.fasta;
#                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.fasta ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.fasta > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.paf;
#                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.paf -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa -b ${outdir}/seqwish_${sample_id}.${graph};
#                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.tmp ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa;
#                    vg paths -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa -p <(cut -f1 ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact_ids.txt) -r > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.vg;
#                    vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.xg ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.vg;
#                    vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.xg -o ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gbwt -P --pass-paths
#                    vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.vg > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.pruned.vg
#                    vg index -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.pruned.vg
#                    #
#                    vg paths -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa -p <(cut -f1 ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget_ids.txt) -r > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.vg;
#                    vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.xg ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.vg;
#                    vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.xg -o ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.gbwt -P --pass-paths
#                    vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.vg > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.pruned.vg
#                    vg index -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.pruned.vg
#                    # map to relevant subgraphs
#                    vg map -F ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.xg -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.alleles.gam
#                    vg map -F ${outdir}/${sample_id}.${graph}.${gene}.off_target_alleles.fasta -x ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.xg -g ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget_alleles.gam
#                    cat ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.alleles.gam ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.offtarget_alleles.gam > ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.gam
#                    vg augment -iS ${outdir}/${sample_id}.${graph}.${gene}.strict_haps.gfa ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.gam > ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.vg
#                    vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.vg > ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.gfa
#                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.tmp ${outdir}/${sample_id}.${graph}.${gene}.hap_specific_aln.gfa;
#            fi
            #mv ${sample_id}.${gene}.exp*gfa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fix.gfa && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fix.gfa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
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
            cp ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt
            cp ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget.txt ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.offtarget.txt
            rm ${outdir}/${gene}.alleles.stats
            rm ${outdir}/${gene}.haps.fasta
        fi
        # parse alignments - probably too permissive currently with ..genotyping.immune_subset.vg for read extraction for complex genes e.g. 3-30
        vg find -x ${graph_base}.xg -l ${workdir}/${gam_file%.gam}.sorted.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam
        vg view -a ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -X | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt
        if [[ $(wc -l <${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt) -ge 1 ]]; then
            # filtering weights in final call set - downweight alleles in a component containing other genes
            # vg chunk -- dont think we need component-level info anymore
            #vg chunk -C -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -b ${outdir}/${sample_id}.${graph}.${gene}.component
            #> ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt;
            #for chunk in $(ls ${outdir}/${sample_id}.${graph}.${gene}.component*vg); do echo ${chunk};
            #    component=$(echo "${chunk##*$gene.}" | sed s/".vg"//g)
            #    alleles_comp=$(vg paths -Lv ${chunk} | grep "${gene}\*" | sort | uniq)
            #    other_genes_comp=$(vg paths -Lv ${chunk} | grep -v "${gene}\*" | grep "IMGT\|OGRDB\|IGv2" | sed s/".*#1#"//g | sort | uniq)
            #    echo $component >> ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt
            #    echo "${alleles_comp[@]}" >> ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt
            #    echo "${other_genes_comp[@]}" >> ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt
            #done
            #sed -i '/^$/d' ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt
            echo "Re-aligning reads to locus-specific haplotype graph"
            vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            # copy mapped reads -- use this for augmentation?
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
            # extra read filtering for complex genes
            if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
                if [ $(vg paths -Lv  ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep -v "${gene}\*" | grep "IMGT\|IGv2\|OGRDB" | wc -l) -gt 0 ]; then
                    echo "Complex locus detected for ${gene} - with multiple other genes in overlapping ASC clusters - filtering out reads aligning to non-gene nodes"
                    sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                    if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt ]; then
#                        Rscript ${bigfoot_dir}/identify_non_haplotype_nodes_complex_locus.R ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                        Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                    else
                        Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                    fi
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg
                    vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
#                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
#                    vg filter -r 0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                    vg filter -r 0.95 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
                    cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
                else
                    echo "Complex locus detected for ${gene} - however the multiple ASC clusters contain only taget gene alleles - retaining all reads"
                    vg filter -r 0 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                fi
            elif [ $(vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep "IMGT\|IGv2\|OGRDB" | grep -v "${gene}\*\|${gene}_\|${gene}" | wc -l) -gt 0 ]; then
                echo "Complex locus detected for ${gene} - single ASC cluster but additional genes present - filtering out reads aligning to non-gene nodes"
                sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt ]; then
#                    Rscript ${bigfoot_dir}/identify_non_haplotype_nodes_complex_locus.R ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                    Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                else
                    Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                fi
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes ]; then
                    echo "Filtering out reads aligning to non-gene nodes"
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg
                    vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                    vg filter -r 0.0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#                else
#                    Rscript ${bigfoot_dir}/identify_non_haplotype_nodes_complex_locus.R ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                fi
            else
                echo "Non-ASC-based analysis for ${gene}"
                vg filter -r 0 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            fi
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            echo "1) Performing inference on haplotype graph labled only with alleles of interest";
            # if variable gene -- limit to OGRDB validated set of alleles?
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "IMGT\|IGv2\|OGRDB" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
            if [ "${valid_alleles}" = true ]; then
                echo "Where possible - retaining path labels only for OGRDB validated set of ${gene_actual} alleles";
                if [[ ${gene} == *"IGHV"* ]]; then
                    echo "IGHV";
                    # joint filtering
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
            elif [ "${valid_alleles}" = "igenotyper" ]; then
                echo "Retaining path labels only for IGenotyper alleles for ${gene_actual}. If you want to limit variable gene inference to OGRDB set - set 'valid_alleles=true'"
                if [ -s ${bigfoot_source}/igenotyper_alleles/alleles.fasta ]; then
                    echo "using saved Igenotyper allele set (${bigfoot_source}/igenotyper_alleles/alleles.fasta)"
                else
                    wget -P ${bigfoot_source}/igenotyper_alleles https://raw.githubusercontent.com/oscarlr/IGenotyper/master/IGenotyper/data/alleles.fasta
                fi
                echo "Retaining path labels only for alleles matching IGenotyper alleles for ${gene} - sequence based matching";
                seqkit grep -n -r -p "${gene}_" ${bigfoot_source}/igenotyper_alleles/alleles.fasta > ${outdir}/${sample_id}.${graph}.${gene}.igeno_alleles.fasta
                vg paths -Fv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.igeno_filtering.fasta
                Rscript ${bigfoot_dir}/igenotyper_allele_matching.R ${outdir}/${sample_id}.${graph}.${gene}.igeno_alleles.fasta
                grep ">" ${outdir}/${sample_id}.${graph}.${gene}.igeno_filtered.fasta | sed s/">"//g | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.alleles
                rm ${outdir}/${sample_id}.${graph}.${gene}.igeno*fasta
            else
                echo "Retaining path labels only for IMGT-derived ${gene_actual} If you want to limit variable gene inference to OGRDB set - set 'valid_alleles=true'"
                grep "${gene}\*" ${outdir}/${sample_id}.${graph}.${gene}.alleles | grep "IMGT" > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
            fi
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
            # filter out reads mapping to nodes specific for non-ogrdb alleles?
#            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep -v -w -f <(sed s/\*/'\\*'/g ${outdir}/${sample_id}.${graph}.${gene}.alleles) | grep "^I\|^OG" > ${outdir}/${sample_id}.${graph}.${gene}.additional_filt_alleles
#            vg paths -Ax ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -p ${outdir}/${sample_id}.${graph}.${gene}.additional_filt_alleles | cut -f6 | sed s/'>'/' '/g | tr ' ' '\n' | sed "/^ *$/d" | sed "/ /d" | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.additional_filt_alleles_nodes
#            vg paths -Ax ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -p ${outdir}/${sample_id}.${graph}.${gene}.alleles | cut -f6 | sed s/'>'/' '/g | tr ' ' '\n' | sed "/^ *$/d" | sed "/ /d" | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.matching_filt_alleles_nodes
#            comm -13 ${outdir}/${sample_id}.${graph}.${gene}.matching_filt_alleles_nodes ${outdir}/${sample_id}.${graph}.${gene}.additional_filt_alleles_nodes > ${outdir}/${sample_id}.${graph}.${gene}.valid.filteringnodes
#            vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -L -N ${outdir}/${sample_id}.${graph}.${gene}.valid.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.vg
#            vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#            vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.txt
#            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.txt ];then
#                vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
#                vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
#                vg filter -r 0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#            fi
# Require a read mapping to nodes specific to valie alleles?
#            vg paths -Ax ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -p ${outdir}/${sample_id}.${graph}.${gene}.alleles | cut -f6 | sed s/'>'/' '/g | tr ' ' '\n' | sed "/^ *$/d" | sed "/ /d" | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.matching_filt_alleles_nodes
#            vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -L -N ${outdir}/${sample_id}.${graph}.${gene}.matching_filt_alleles_nodes > ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.vg
#            vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#            vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.txt
#            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.txt ];then
#                vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -n -f ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.fastq
#                vg map -f ${outdir}/${sample_id}.${graph}.${gene}.valid.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.valid.prefilt.gam
#                vg filter -r 0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.valid.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#            fi
            vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
            # gene-specific read depth for minimum strain-level coverage
            depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
            # min_strain_depth=$(bc -l <<< "scale=3;${depth_locus}/$gene_min_len")
            min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.05"| awk '{printf("%d\n",$1 + 0.5)}')
#            min_strain_depth=$(bc -l <<< "scale=2;(${depth_locus}*0.05)*$(echo "${asc_cluster[@]}" | wc -l)"| awk '{printf("%d\n",$1 + 0.5)}')
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
                echo "Round 2: ILP-based inference + allele abundance refinement"
                if [ $(grep ">" trimmed_contigs.fasta | grep -v "ab=0.000" | wc -l) -gt 1 ]; then
                    max_strains=$(grep ">" trimmed_contigs.fasta | grep -v "ab=0.000" | wc -l)
                else
                    max_strains=2
                fi
                python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --ilp --max_strains ${max_strains} --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c ${min_strain_depth} --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
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
                        echo "Round 3: No strains with estimated allele coverage >= ${min_strain_depth} - using basal threshold (0.01) + ILP";
                        python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --ilp --max_strains ${max_strains} --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c 0.01 --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
#                          vg depth -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.vg
                        # inferred coverage seems to be using the lowest abundance node... should be estimating the average - perhaps weighted path makes sense here vs. just averaging.
                        # use below chunks - select rows/nodes associated with alleles of interest - set lowest coverage node to avg and rerun
                        # mean_node_depth=$(perl -F: -lane '$total += $F[1]; END{print $total/$.}' ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt)
                        # mean_node_depth=$(bc -l <<< "scale=3;${mean_node_depth}/$(echo "${asc_cluster[@]}" | wc -l)")
#                        python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c 0.01 --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa

                    fi
                    mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.fasta ;
                    mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.rel.genome_graph.gfa;
                    rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
                    if [ -s ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.fasta ]; then
                    # account for ${outdir}/${sample_id}.${graph}.${gene}.component.summary.txt to up/downweight alleles influenced by non-gene of interest
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
            if [ "${assoc_testing}" = true ]; then
                echo "2) Augmenting annotated post-flow inference graph with reads for association testing";
                sed s/' '/_/g ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta > ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta ]; then
                    #seqkit grep -r -p "chm" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
                    seqkit grep -r -p "chm13" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.ref.fasta
                    if [ -s ${outdir}/${sample_id}.${graph}.${gene}.ref.fasta ]; then
                        echo "Using CHM13 reference backbone"
                        ref_backbone="chm"
                    else 
                        echo "Using GRCh38 reference backbone"
    #                    seqkit grep -r -p "grch" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.ref.fasta
                        ref_backbone="grch"
                    fi
                    if [ $(grep ${ref_backbone} ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | cut -f1 -d':' | wc -l) -gt 1 ]; then
                        echo "Likely orphon involvement - ensuring we use best matching reference locus"
                        seqkit grep -r -p ${ref_backbone} ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | minimap2 -x asm20 --secondary=no -c - ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt
                        for chr_loc in $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                            chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                            grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta
                            seqkit grep -p $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | cat - ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta
                            minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf
                            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa
                        done
                        if [ $(ls ${outdir}/${sample_id}.${graph}.${gene}.*.genome_graph_ref.gfa | wc -l) -gt 1 ]; then
                            echo "combining relevant subgraphs"
                            ls ${outdir}/${sample_id}.${graph}.${gene}.*.genome_graph_ref.gfa > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                            odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.og
                            odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.og -g > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                            #delete temp files
                            xargs rm < ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
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
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
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
                    vg augment -m ${augment_cov} -q 5 -Q 60 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg;
                    vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
                    # gaf of reads aligned to augmented graph - gafpack to get coverage/depth
                    vg convert -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf
                    # paired-end, append identified to duplicate read ids
                    awk -F'\t' -v OFS='\t' 'NR==FNR{f[$1]++} NR>FNR{if (f[$1]>1) {$1=$1 "-" (++count[$1])}; print}' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.dedup.gaf
                    # add read ids for testing/validation:
                    vg augment -B ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg -F ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.dedup.gaf > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg
                    vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                    # break up gfa file because larger textfiles are problematic:
                    head -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.header
                    grep "^P" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.plines
                    grep "^S" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.slines
                    grep "^L" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.llines
#                    cat ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.noplines ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.plines > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                    gafpack -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -a ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf -l > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage
                    Rscript ${bigfoot_dir}/augment_graph_wdepth.R ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage >> ${outdir}/${sample_id}.putative_variants.csv
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths
                    grep "chm\|grch\|${gene}" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths
                    vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.pg
                    vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa
                    grep "^P" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.plines
                    grep "^P" -v ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.noplines
                    cat ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.noplines ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.plines > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa
                    # should really re-align reads and then pack/call variants so we can incorporate depth info
                    vg deconstruct -a -n -L 1 -P "chm" -P "grch" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.coding.vcf
                    rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.header
                    rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.slines
                    rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.llines
                    rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.plines
                    # also output something indicating novel variant identified - read support for node with no reference path
                    # chop into single nucleotide nodes - then only merge nodes with equal depth/coverage in the R script after embedding reads
#                    vg mod -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
#                    vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
#                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                    # Additional inference - cleaning complex genes + allele inference
                    careful=false
                    if [ "${careful}" = true ]; then
                        echo "Careful mode - Additional round of flow-based inference"
                        if [ $(grep ${gene} ${bigfoot_dir}../custom_beds/complex_genes.txt | wc -l) -gt 0 ];then
                            echo "Potential overlap with complex gene - additional round of flow-based inference"
                            complex_gene=true
                        fi
                        if [ "${complex_gene}" = true ]; then
                            mkdir -p ${outdir}/${gene}_careful; cd ${outdir}/${gene}_careful
                            echo "Complex gene - additional round of flow-based inference"
                            gene_aug_orig_graph=${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                            vg convert -g ${gene_aug_orig_graph} -p > ${gene}.index.pg
                            vg index -t 16 -L -x ${gene}.index.xg ${gene}.index.pg;
                            vg gbwt -x ${gene}.index.xg -o ${gene}.index.gbwt -P --pass-paths
                            vg prune -u -g ${gene}.index.gbwt -k 31 -m ${gene}.index.node_mapping ${gene}.index.pg > ${gene}.index.pruned.vg
                            vg index -g ${gene}.index.gcsa -f ${gene}.index.node_mapping ${gene}.index.pruned.vg
                            vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam -d ${gene}.index -1 ${gene}.index.gbwt > ${gene}.aug.gam 
                            vg filter -r 0.95 -P -s 1 -x ${gene}.index.xg -D 0 -fu -t 4 ${gene}.aug.gam -v > ${gene}.aug.filt.gam
                            vg paths -Lx ${gene}.index.pg | grep -v "grch\|chm" > ${gene}.index.alleles
                            vg paths -r -p ${gene}.index.alleles -x ${gene}.index.pg > ${gene}.index.careful.vg
                            vg view -a ${gene}.aug.filt.gam > ${gene}.index.careful.aln.json
                            vg convert -fW ${gene}.index.careful.vg > ${gene}.index.careful.gfa
                            python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${gene}.index.careful -m 0
                            python3 ${bigfoot_dir}/vg-flow_immunovar.py --ilp --max_strains $(wc -l ${gene}.index.alleles | cut -f1 -d' ') --careful --optimization_approach ${opt} --min_depth 0 --trim 0 -m 0 -c 0.1 --remove_included_paths 0 ${gene}.index.careful.node_abundance.txt ${gene}.index.careful.final.gfa
                            mv trimmed_contigs.fasta ${gene}.index.careful.contigs.fasta; mv haps.final.fasta ${gene}.index.careful.haps.final.fasta ; mv genome_graph.gfa ${gene}.index.careful.genome_graph.gfa
                            Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${gene}_careful/${gene}.index.careful.contigs.fasta
                            echo "Careful allele-level abundance estimation completed for ${gene} ::"
                            echo "Saving previous sample-specific allele fasta file"
                            sed s/:path.*:path/:path/g ${gene}.index.careful.haps.final.annot.fasta > ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.careful.fasta
# add reference paths from original augmented graph - reconstruct augmented graph
                            echo "Reconstructing augmented graph using cleaned set of alleles"
                            vg paths -Fx ${gene_aug_orig_graph} | seqkit grep -r -p "chm|grch" > ${gene}.index.ref_paths.fasta
                            cat ${gene}.index.ref_paths.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.careful.fasta | sed s/' '/_/g > ${gene}.careful.w_ref_paths.fasta                            
                            minimap2 -x asm20 -t 16 -c -X ${gene}.careful.w_ref_paths.fasta ${gene}.careful.w_ref_paths.fasta > ${gene}.careful.w_ref_paths.paf
                            seqwish -s ${gene}.careful.w_ref_paths.fasta -p ${gene}.careful.w_ref_paths.paf -g ${gene}.careful.w_ref_paths.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${gene}.careful.w_ref_paths.gfa -o ${gene}.careful.w_ref_paths.tmp; mv ${gene}.careful.w_ref_paths.tmp ${gene}.careful.w_ref_paths.gfa
                            vg mod -n -U 10 -c ${gene}.careful.w_ref_paths.gfa -X 256 > ${gene}.careful.w_ref_paths.vg
                            vg convert -p ${gene}.careful.w_ref_paths.vg > ${gene}.careful.w_ref_paths.pg
                            vg convert -fW ${gene}.careful.w_ref_paths.vg > ${gene}.careful.w_ref_paths.gfa
                            gfaffix ${gene}.careful.w_ref_paths.gfa -o ${gene}.careful.w_ref_paths.gfa.tmp; mv ${gene}.careful.w_ref_paths.gfa.tmp ${gene}.careful.w_ref_paths.gfa
                            vg index -t 16 -L -x ${gene}.careful.w_ref_paths.xg ${gene}.careful.w_ref_paths.pg;
                            vg gbwt -x ${gene}.careful.w_ref_paths.xg -o ${gene}.careful.w_ref_paths.gbwt -P --pass-paths
                            vg gbwt -x ${gene}.careful.w_ref_paths.xg -g ${gene}.careful.w_ref_paths.gbz --gbz-format -P --pass-paths;
                            vg prune -u -g ${gene}.careful.w_ref_paths.gbwt -k 31 -m ${gene}.careful.w_ref_paths.node_mapping ${gene}.careful.w_ref_paths.pg > ${gene}.careful.w_ref_paths.pruned.vg
                            vg index -g ${gene}.careful.w_ref_paths.gcsa -f ${gene}.careful.w_ref_paths.node_mapping ${gene}.careful.w_ref_paths.pruned.vg
                        #   map locus-associated reads to augmented/annotated graph
                            vg map -N ${sample_id}.${graph}.${gene} -G ${gene}.aug.filt.gam -x ${gene}.careful.w_ref_paths.xg -g ${gene}.careful.w_ref_paths.gcsa -1 ${gene}.careful.w_ref_paths.gbwt -t 4 -M 1 > ${gene}.careful.w_ref_paths.gam
                            vg filter -r 0.975 -P -s 1 -q 60 -x ${gene}.careful.w_ref_paths.xg -D 0 -fu -t 4 ${gene}.careful.w_ref_paths.gam -v > ${gene}.careful.w_ref_paths.filt.gam
                            vg depth --gam ${gene}.careful.w_ref_paths.filt.gam ${gene}.careful.w_ref_paths.xg > ${gene}.careful.w_ref_paths.filt.depth;
                            depth_aug=$(awk -F ' ' '{print $1}' ${gene}.careful.w_ref_paths.filt.depth)
                            aug_depth=$(bc -l <<< "scale=2;${depth_aug}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
                            if [ "${aug_depth}" -gt 3 ]; then
                                augment_cov=${aug_depth}
                            else
                                augment_cov=3
                            fi
                            echo "Minimum coverage to add breakpoint: ${augment_cov} (3 <--> 10% of strain depth)"
                            vg augment -m ${augment_cov} -q 5 -Q 60 ${gene}.careful.w_ref_paths.pg ${gene}.careful.w_ref_paths.filt.gam -A ${gene}.careful.w_ref_paths.augmented.gam > ${gene}.careful.w_ref_paths.augmented.vg;
                            vg convert -p ${gene}.careful.w_ref_paths.augmented.vg > ${gene}.careful.w_ref_paths.augmented.pg
                            vg mod -c ${gene}.careful.w_ref_paths.augmented.pg > ${gene}.careful.w_ref_paths.augmented.tmp && mv ${gene}.careful.w_ref_paths.augmented.tmp ${gene}.careful.w_ref_paths.augmented.pg
                            vg convert -fW ${gene}.careful.w_ref_paths.augmented.pg > ${gene}.careful.w_ref_paths.augmented.gfa
                            gfaffix ${gene}.careful.w_ref_paths.augmented.gfa -o ${gene}.careful.w_ref_paths.augmented.gfa.tmp; mv ${gene}.careful.w_ref_paths.augmented.gfa.tmp ${gene}.careful.w_ref_paths.augmented.gfa
                            cp ${gene}.careful.w_ref_paths.augmented.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.careful.gfa
                            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.careful.gfa ];then
                                echo "Using careful augmented graph for association testing"
                                mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.careful.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
                            fi
                            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.careful.fasta ];then
                                echo "Reporting careful fasta alleles"
                                mv ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.careful.fasta ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta
                            fi
                            if [ "${de_novo}" = true ]; then
                                echo "De-novo allele inference... WiP"
                            cd ${outdir}
                            else
                                echo "No novel allele inference requested - to perform novel allele inference set de_novo=true"
                                cd ${outdir}
                            fi
                            rm -rf ${outdir}/${gene}_careful
                        else
                            echo "No complex gene overlap detected - skipping additional round of fslow-based inference"
                        fi
                    else
                        echo "No additional round of flow-based inference requested - to perform additional round of flow-based inference set careful=true"
                        if [ "${de_novo}" = true ]; then
                            echo "De-novo allele inference... WiP"
                            # map - remove nodes with 0/low coverage - unitigs - stitch together using contiguous abundance...
                            # at most number of alleles = ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta + number of nodes without a path
                            vg view -X ${gene}.careful.w_ref_paths.augmented.gam | seqkit fq2fa - > ${gene}.careful.w_ref_paths.augmented.fasta
                            odgi unitig -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -t 151 > ${gene}.careful.unitigs.fa 
                            Bifrost build -r ${gene}.careful.unitigs.fa -k 31 -o ${gene}.careful.bifrost
                            gunzip ${gene}.careful.bifrost*gfa.gz
                            ~/tools/SPAdes-4.0.0-Linux/bin/pathracer ${gene}.careful.w_ref_paths.augmented.fasta ${gene}.careful.bifrost.gfa --nt --annotate-graph --output ${gene}.careful.pathracer

                        else
                            echo "No novel allele inference requested - to perform novel allele inference set de_novo=true"
                        fi
                    fi
                else
                    echo "Insuffient coverage to infer alleles for $gene - not trying to embed variation"
                fi
            else
                echo "No association testing performed for ${gene} - set assoc_testing=true to embed reads/prep locus graph for association testing"
            fi
            # remove files we dont need anymore
            ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v "${gene}.genome_graph_ref.augmented.*gfa\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|.vcf\|node_abundance\|depth" > ${outdir}/${sample_id}_${gene}_files.txt
            #ls ${outdir}/${gene}\.*  | grep "haps.fasta\|alleles" >> ${outdir}/${sample_id}_${gene}_files.txt
            ls ${outdir}/*${gene}* | grep "asc_" >> ${outdir}/${sample_id}_${gene}_files.txt
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

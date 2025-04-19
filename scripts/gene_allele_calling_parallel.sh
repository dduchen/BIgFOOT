#!/bin/bash

cd ${datadir}
if [ "${prep_locus_graphs}" = "parse_complex" ]; then
    echo "Parsing local graphs across complex genes - ensuring graphs contain unique/gene-specific local haplotypes"
    if [ -s ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.succinct_fix.txt ]; then
        echo "Complex local graphs already parsed - skipping";
    else
        ls ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.*.succinct_locus.xg | grep -v "custom" | grep -f ${bigfoot_dir}/../custom_beds/complex_genes.txt > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.txt;
        for compgene in $(cat ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.txt);do vg convert -fW ${compgene} > ${compgene%.xg}.gfa; 
        done
        sed -i s/".xg$"/".gfa"/g ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.txt;
        odgi squeeze -f ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.txt -O -o ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.og;
        odgi paths -fi ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.og | seqkit grep -v -n -rp "IMGT|OGRDB|HLA|KIR" | seqkit rmdup -n -D ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt > test_tmp.fasta;
        rm test_tmp.fasta;
        cut -f2 ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt | sed 's/ //' | tr , '\n' | sort | uniq > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt.tmp && mv ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt.tmp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt;
        sed s/".*wg_immunovar."//g ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.txt | sed s/".succinct.*"/""/g > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.succinct_fix.txt;
        for tmp_graph in $(ls ${genotyping_nodes_dir}/gene_graphs/*.xg | grep -f ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.succinct_fix.txt | grep -v "custom");do vg convert -fW ${tmp_graph} > ${tmp_graph%.xg}.tmp.gfa;
            vg paths -Lx ${tmp_graph%.xg}.tmp.gfa | grep -f ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.complex_genes.dup_haps.txt > ${tmp_graph%.xg}.dup_paths.txt;
            if [[ $(wc -l <${tmp_graph%.xg}.dup_paths.txt) -ge 1 ]]; then
                cp ${tmp_graph%.xg}.tmp.gfa ${tmp_graph%.xg}_legacy.gfa;
                gfaffix ${tmp_graph%.xg}.tmp.gfa -o ${tmp_graph%.xg}.gfa;
                echo "removing duplicate paths from $(echo ${tmp_graph%.xg} | sed s/".*wg_immunovar."//g)";
                vg paths -x ${tmp_graph%.xg}.gfa -d -p ${tmp_graph%.xg}.dup_paths.txt | vg mod -N -nU 10 -X 32 - > ${tmp_graph%.xg}.vg;
                vg convert -p ${tmp_graph%.xg}.vg > ${tmp_graph%.xg}.pg;
                vg convert -fW ${tmp_graph%.xg}.pg > ${tmp_graph%.xg}.gfa;
                vg index -t 16 -L -x ${tmp_graph} ${tmp_graph%.xg}.pg;
                vg gbwt -x ${tmp_graph%.xg}.xg -o ${tmp_graph%.xg}.gbwt -P --pass-paths;
                vg prune -u -g ${tmp_graph%.xg}.gbwt -k 31 -m ${tmp_graph%.xg}.node_mapping ${tmp_graph%.xg}.pg > ${tmp_graph%.xg}.pruned.vg;
                vg index -g ${tmp_graph%.xg}.gcsa -f ${tmp_graph%.xg}.node_mapping ${tmp_graph%.xg}.pruned.vg;
                rm ${tmp_graph%.xg}.node_mapping; rm ${tmp_graph%.xg}.pruned.vg;
                rm ${tmp_graph%.xg}.tmp.gfa;
            else
                rm ${tmp_graph%.xg}.tmp.gfa;
            fi
        done
    fi
    return
fi
#
gene=${each%.nodes.txt}
gene=${gene%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
#
loci=$(vg paths -Lv ${graph_base}.xg | grep "#1#${gene}\*" | cut -f1 -d"#" | sort | uniq | grep "IMGT\|OGRDB\|HLA\|KIR")
echo "$gene --> $(echo $loci)"
if [[ $(echo $loci | grep "OGRDB" | wc -l) -ge 1 ]]; then
    echo "OGRDB alleles"
    loci="IMGT"
fi
#
unset asc_cluster #ensure no carryover of variables from previous gene analysis
unset complex_gene #ensure no carryover of variables from previous gene analysis
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
        echo "IMGT: Immunoglobulin/T-cell receptor gene inference";
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
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gbwt ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.pg ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gcsa ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa
                cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa.lcp
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
            # could be many more alleles outside the genes we're most focused on - constrain to specific alleles spanning the relevant ASC clusters
            asc_cluster=$(grep "${gene}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f4 | sed s/'\*.*'//g | sort | uniq);
            if [[ -n ${asc_cluster} ]]; then
                echo "limiting candidate alleles to ASC-constrained set"
                for asc_specific in ${asc_cluster[@]}; do echo ${asc_specific};
                    grep "${asc_specific}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f1 | sed s/'*'/'\\*'/g >> ${outdir}/asc_relevant_allele_for_${gene}
                    grep "${asc_specific}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f1 | sed s/'*.*'/''/g | sed s/".*#1#"/""/g >> ${outdir}/asc_relevant_genes_for_${gene}
                done
                # retain alleles from specific genes
                seqkit grep -r -n -f <(cut -f1 ${outdir}/asc_relevant_genes_for_${gene} | sort | uniq) ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.fasta.tmp 
                # retain specific alleles
                #seqkit grep -r -n -f <(cut -f1 ${outdir}/asc_relevant_allele_for_${gene} | sed s/".*#1#"/""/g | sort | uniq) ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.fasta.tmp 
                if [ $(grep ">" ${outdir}/${gene}.alleles.fasta.tmp | wc -l) -gt 1 ]; then
                    mv ${outdir}/${gene}.alleles.fasta.tmp ${outdir}/${gene}.alleles.fasta
                else 
                    rm ${outdir}/${gene}.alleles.fasta.tmp
                fi
            else 
                echo "No ASC table entry for ${gene}"
            fi
            seqkit grep -r -p ${gene}"\*" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.exact.fasta
            seqkit grep -r -v -p ${gene}"\*" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.offtarget.fasta
            # set min length of haplotypes = 100, max length should scale with allele length --> then downsample haplotypes
            seqkit stats ${outdir}/${gene}.alleles.exact.fasta > ${outdir}/${gene}.alleles.stats
            gene_min_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f6 -d' ' | sed s/","//g)
            gene_max_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f8 -d' ' | sed s/","//g)
            ##############################
            # use local haplotypes #
            ##############################
            seqkit grep -v -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.haps.fasta
            # trim terminal repeats in the local haplotype sequences
            tr-trimmer -i ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.trimmed.fasta;
            if [ -s ${outdir}/${gene}.haps.trimmed.fasta ]; then
               mv ${outdir}/${gene}.haps.trimmed.fasta ${outdir}/${gene}.haps.fasta
            fi
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
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                            # Alleles might align best to a speciifc haplotype - but need ensure the haplotype belongs in the chromosome-specific subgraph
                            if [ $(grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep ${chr_loc} | wc -l) -ge 1 ]; then
                                seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                            fi
                            if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta | wc -l ) -gt 2 ]; then
                                vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.rough.gfa
                                odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa;
                                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.tmp.gfa;
                                mv ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.tmp.gfa ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa
                                rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.rough.gfa;
                            else
                                # all-to-all graph induction (can lead to multi-component subgraphs)
                                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf
                                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa
                            fi
                            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.paths
                            paths_count=$(wc -l < ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.paths | awk '{print $1}')
                            fasta_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta | wc -l | awk '{print $1}')
                            if [ "$paths_count" -lt "$fasta_count" ]; then
                                echo "paths missing - retrying minimap2 all-to-all construction"
                                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf
                                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.${chr_graph_tmp}.gfa
                            fi
                        done
                        ls ${outdir}/${sample_id}.${graph}.${gene}.${clust}.*.gfa | grep -v "haplotypes" - > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt
                        odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og 
                        odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa
                        #delete temp files
                        rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.graphs_to_squeeze.txt; rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.chr*
                    else
                        echo "Single locus for this cluster of ${gene}"
                        seqkit grep -r -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta
                        if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt ]; then
                            grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt | cut -f2 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotype_assignments.patterns.txt
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta
                        fi
                        if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta  | wc -l ) -gt 2 ]; then
                            vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.rough.gfa
                            odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa;
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa.tmp;
                            mv ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa;
                            rm ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.rough.gfa;
                        else 
                            minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf
                            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa
                        fi
                        vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.paths
                        paths_count=$(wc -l < ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.paths | awk '{print $1}')
                        fasta_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta | wc -l | awk '{print $1}')
                        if [ "$paths_count" -lt "$fasta_count" ]; then
                            echo "paths missing - retrying minimap2 all-to-all construction"
                            minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf
                            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${clust}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${clust}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${clust}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${clust}.haplotypes.gfa
                        fi
                    fi
                done
                echo "Constructing pan-cluster graph for ${gene}"
                ls ${outdir}/${sample_id}.${graph}.${gene}.*.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -s "#" -O -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og 
                odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
                # final inference step after filtering out spurious reads - consolidated gene-specific graph when alleles spread across clusters:
                # --> add something during inference: if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
                cat ${outdir}/${gene}.alleles.exact.fasta ${outdir}/${gene}.haps.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta
                cat ${outdir}/${gene}.alleles.exact.fasta ${outdir}/${gene}.alleles.offtarget.fasta ${outdir}/${gene}.haps.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta
                # progressive alignment approach:
                if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta | wc -l ) -gt 2 ]; then
                    vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.rough.gfa;
                    odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa;
                    rm ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.rough.gfa;
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa
                    # confirmation / rescue if needed
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paths
                    paths_count=$(wc -l < ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paths | awk '{print $1}')
                    fasta_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta | wc -l | awk '{print $1}')
                    if [ "$paths_count" -lt "$fasta_count" ]; then
                        echo "paths missing - retrying minimap2 all-to-all construction"
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp;
                        mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa
                    fi                    
                    # for read filtering - rough version including off-target alleles
                    vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.rough.gfa;
                    odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa;
                    rm ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.rough.gfa;
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa;
                    # confirmation / rescue if needed
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paths
                    paths_count=$(wc -l < ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paths | awk '{print $1}')
                    fasta_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta | wc -l | awk '{print $1}')
                    if [ "$paths_count" -lt "$fasta_count" ]; then
                        echo "paths missing - retrying minimap2 all-to-all construction"
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paf -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp;
                        mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa
                    fi                    
                    #
                else 
                    # all-to-all graph induction (can lead to multi-component subgraphs)
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.paf -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gfa
                    # for read filtering - rough version including off-target alleles
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.fasta > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.paf -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa
                fi
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
                #
                vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg
                vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.vg
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.vg > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa
                vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg
                vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg;
                vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg -o ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt -P --pass-paths
                vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg > ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pruned.vg
                vg index -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pruned.vg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gfa ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gfa
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gbwt
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.xg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.pg
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gcsa
                cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa.lcp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}.succinct_pancluster.gcsa.lcp
                #delete temp files
                xargs rm < ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
            else
                echo "Single locus/cluster for ${gene}"
                if [ -s ${outdir}/${gene}_haps/haps.matching.txt ]; then
                    if [ -s ${outdir}/${gene}_haps/haps.matching.exact.txt ]; then
                        echo "Haplotypes exist with exact matches to alleles"
                        seqkit grep -n -f ${outdir}/${gene}_haps/haps.matching.exact.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    else
                        echo "Haplotypes exist with exact matches to ASC-related genes only"
                        seqkit grep -n -f ${outdir}/${gene}_haps/haps.matching.txt ${outdir}/${gene}.haps.fasta > $outdir/${gene}_haps/${gene}.haps.fasta
                    fi
                    seqkit split --quiet -i ${outdir}/${gene}_haps/${gene}.haps.fasta -f
                    mkdir -p ${outdir}/${gene}_haps_final
                    ${tools_dir}/Assembly-dereplicator/dereplicator.py --distance 0.000001 --count $(grep ">" ${outdir}/${gene}.alleles.exact.fasta | wc -l) $outdir/${gene}_haps/${gene}.haps.fasta.split/ $outdir/${gene}_haps_final
                    cat ${outdir}/${gene}_haps_final/*fasta > ${outdir}/${gene}.haps.fasta
                    seqkit stats ${outdir}/${gene}.haps.fasta
                    rm -rf $outdir/${gene}_haps_final ; rm -rf $outdir/${gene}_haps
                    cat ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                    seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
                else
                    echo "No haplotypes containing the alleles"
                    rm -rf $outdir/${gene}_haps
                fi
               # this still checks for multiple chromosomes - if single chromosome, then this shouldnt change things
                for chr_loc in $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                    chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                    grep ${chr_loc%:*} ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                    seqkit grep -n -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    if [ $(cut -f2 ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | sort | uniq | wc -l) -le 1 ]; then
                        echo "adding all alleles to subgraph"
                        seqkit grep -r -p "${gene}|IMGT|OGRDB|IGv2" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    fi
                    seqkit rmdup -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    #grep -r ${chr_graph_tmp} ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep -v "${gene}\|IMGT\|OGRDB\|IGv2" - | cut -f1 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.haplotype_assignments.patterns.txt
                    #seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta | wc -l ) -gt 2 ]; then
                        vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -w 5000 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.rough.gfa;
                        odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                        rm ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.rough.gfa;
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                    else
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf -g ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                    fi
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paths
                    vg paths -Fv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa | seqkit grep -s -p 'N' - > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paths.fasta
                    paths_count=$(wc -l < ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paths | awk '{print $1}')
                    fasta_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta | wc -l | awk '{print $1}')
                    fasta_N_count=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paths.fasta | wc -l | awk '{print $1}')
                    if [ "$paths_count" -lt "$fasta_count" ] || [ "$fasta_N_count" -gt 0 ]; then
                        echo "paths missing or contain ambiguous bases - retrying minimap2 all-to-all construction"
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.paf -g ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.tmp ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.gfa
                    fi
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
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fix.gfa && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fix.gfa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg;
            vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg;
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -P --pass-paths;
            vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz --gbz-format -P --pass-paths;
            if [ "${complex_gene}" = true ]; then
                vg prune ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg;
                vg index -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg;
            else
                vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg;
                vg index -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg;
            fi
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
      # Combined succinct locus graph for most complex genes - remove haplotypes containing exact matches to more than a single gene - could odgi squeeze all the complex gene succinct files to detect multiply occurring local haplotype sequences (seqkit rmdup)
      # If this improves things incorporate for other multi-asc genes (IGHV4-34/4-30-4, IGKV1-17, IGLV2-14/2-23)
        #if [[ $(echo $gene | grep "IGHV4-4$\|IGHV4-61$\|IGHV4-59$" | wc -l) -ge 1 ]]; then
        #    ls ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGH*.succinct_locus.gfa | grep -v "custom" | grep "V4-4.s\|V4-59.s\|V4-61.s" > ${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt
        #    if [[ $(wc -l <${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt) -ge 3 ]]; then
#       #         if [ -s ${genotyping_nodes_dir}/gene_graphs/*${gene}_custom_pancluster.vg ]; then
        #            echo "Custom pancluster graph already exists for ${gene}"
#                else
#                    odgi squeeze -f ${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt -O -o ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.og;  
#                    odgi sort -i ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.og --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.gfa
#                    grep "^P" ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.gfa | cut -f2 | sort | grep -v "IGHV" | uniq -c | grep " 2 " | sed s/".* 2 "//g > ${genotyping_nodes_dir}/gene_graphs/dup_paths.txt
#                    grep "^P"  ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-4.succinct_locus.gfa | cut -f2 | sort | grep -v "IGHV" > ${genotyping_nodes_dir}/gene_graphs/ighv44_hap_paths.txt
#                    grep "^P"  ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-59.succinct_locus.gfa | cut -f2 | sort | grep -v "IGHV" > ${genotyping_nodes_dir}/gene_graphs/ighv459_hap_paths.txt
#                    grep "^P"  ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-61.succinct_locus.gfa | cut -f2 | sort | grep -v "IGHV" > ${genotyping_nodes_dir}/gene_graphs/ighv461_hap_paths.txt
#                    vg paths -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.gfa -d -p ${genotyping_nodes_dir}/gene_graphs/dup_paths.txt | vg mod -N -nU 10 -X 32 - > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.vg
#                    vg paths -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.vg -d -p <(cat ${genotyping_nodes_dir}/gene_graphs/ighv459_hap_paths.txt ${genotyping_nodes_dir}/gene_graphs/ighv44_hap_paths.txt) | vg mod -N -nU 10 -X 32 - > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-61_custom_pancluster.vg
#                    vg paths -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.vg -d -p <(cat ${genotyping_nodes_dir}/gene_graphs/ighv459_hap_paths.txt ${genotyping_nodes_dir}/gene_graphs/ighv461_hap_paths.txt) | vg mod -N -nU 10 -X 32 - > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-4_custom_pancluster.vg
#                    vg paths -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.vg -d -p <(cat ${genotyping_nodes_dir}/gene_graphs/ighv44_hap_paths.txt ${genotyping_nodes_dir}/gene_graphs/ighv461_hap_paths.txt) | vg mod -N -nU 10 -X 32 - > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV4-59_custom_pancluster.vg
#                    for tmp_graph in $(ls ${genotyping_nodes_dir}/gene_graphs/*custom_pancluster.vg);do echo ${tmp_graph}; 
#                        vg convert -p ${tmp_graph} > ${tmp_graph%_pancluster*}.succinct_locus.pg
#                        vg convert -fW ${tmp_graph%_pancluster*}.succinct_locus.pg > ${tmp_graph%_pancluster*}.succinct_locus.gfa
#                        vg index -t 16 -L -x ${tmp_graph%_pancluster*}.succinct_locus.xg ${tmp_graph%_pancluster*}.succinct_locus.pg;
#                        vg gbwt -x ${tmp_graph%_pancluster*}.succinct_locus.xg -o ${tmp_graph%_pancluster*}.succinct_locus.gbwt -P --pass-paths
#                        vg prune -u -g ${tmp_graph%_pancluster*}.succinct_locus.gbwt -k 31 -m ${tmp_graph%_pancluster*}.succinct_locus.node_mapping ${tmp_graph%_pancluster*}.succinct_locus.pg > ${tmp_graph%_pancluster*}.succinct_locus.pruned.vg
#                        vg index -g ${tmp_graph%_pancluster*}.succinct_locus.gcsa -f ${tmp_graph%_pancluster*}.succinct_locus.node_mapping ${tmp_graph%_pancluster*}.succinct_locus.pruned.vg
#                    done
#                fi
#                unset tmp_graph
#                for tmp_graph in $(ls ${genotyping_nodes_dir}/gene_graphs/*${gene}*.xg | grep -v "custom");do vg convert -fW ${tmp_graph} > ${tmp_graph%.xg}.tmp.gfa;
#                    vg paths -Lx ${tmp_graph%.xg}.tmp.gfa | grep -f ${genotyping_nodes_dir}/gene_graphs/dup_paths.txt > ${genotyping_nodes_dir}/gene_graphs/dup_paths_${gene}.txt;
#                    if [[ $(wc -l <${genotyping_nodes_dir}/gene_graphs/dup_paths_${gene}.txt) -ge 1 ]]; then
#                        cp ${tmp_graph%.xg}.gfa ${tmp_graph%.xg}_legacy.gfa;
#                        gfaffix ${tmp_graph%.xg}.tmp.gfa -o ${tmp_graph%.xg}.gfa;
#                        echo "removing duplicate paths from ${tmp_graph%.xg}"
#                        vg paths -x ${tmp_graph%.xg}.gfa -d -p ${genotyping_nodes_dir}/gene_graphs/dup_paths.txt | vg mod -N -nU 10 -X 32 - > ${tmp_graph%.xg}.vg;
#                        vg convert -p ${tmp_graph%.xg}.vg > ${tmp_graph%.xg}.pg;
#                        vg convert -fW ${tmp_graph%.xg}.pg > ${tmp_graph%.xg}.gfa;
#                        vg index -t 16 -L -x ${tmp_graph} ${tmp_graph%.xg}.pg;
#                        vg gbwt -x ${tmp_graph%.xg}.xg -o ${tmp_graph%.xg}.gbwt -P --pass-paths;
#                        vg prune -u -g ${tmp_graph%.xg}.gbwt -k 31 -m ${tmp_graph%.xg}.node_mapping ${tmp_graph%.xg}.pg > ${tmp_graph%.xg}.pruned.vg;
#                        vg index -g ${tmp_graph%.xg}.gcsa -f ${tmp_graph%.xg}.node_mapping ${tmp_graph%.xg}.pruned.vg;
#                    else
#                        rm ${tmp_graph%.xg}.tmp.gfa;
#                    fi
#                done
#            fi
#            custom_succinct=false
#            if [ "${custom_succinct}" = true ]; then
#                if [[ $(echo $gene | grep "IGHV4-4$\|IGHV4-61$\|IGHV4-59$" | wc -l) -ge 1 ]]; then
#                    ls ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGH*.succinct_locus.gfa | grep -v "custom" | grep "V4-4.s\|V4-59.s\|V4-61.s" > ${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt
#                    vg paths -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.gfa -d -p ${genotyping_nodes_dir}/gene_graphs/dup_paths.txt | vg mod -N -nU 10 -X 32 - > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44_459_461_graphs.vg
#                    if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.gcsa ];then
#                        echo "Trying cross-complex gene succinct graph for ${gene}"
#                        cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt
#                        # cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt
#                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt
#                        cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg
#                        # cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
#                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg
#                        cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg
#                        # cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
#                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.pg
#                        cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa
#                        # cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa
#                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa
#                        cp ${genotyping_nodes_dir}/gene_graphs/${graph}.${gene}_custom.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp
#                        # cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp
#                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa.lcp
#                    else
#                        echo "Need existing graphs for all members of ${gene} co-cluster - rerun inference on all genes"
#                    fi    
#                fi
 #           fi
        #    fi
        #fi
        if [ "${prep_locus_graphs}" = true ]; then
            echo "$gene graphs already prepped"
            return
        fi
        if [ "${prep_locus_graphs}" = false ]; then
            # custom graphs
            if [ -s ${genotyping_nodes_dir}/gene_graphs/${graph}.IGHV44-459-461_custom.succinct_locus.gcsa ];then
            else 
                echo "Need to parse IGHV4-4, IGHV4-59, and IGHV4-61 graphs for custom locus graph"
                ls ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGH*.succinct_locus.gfa | grep -v "custom" | grep "V4-4.s\|V4-59.s\|V4-61.s" > ${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt;
                odgi squeeze -f ${genotyping_nodes_dir}/gene_graphs/ighv44_459_461_graphs.txt -O -o ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.og;
                odgi sort -i ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.og --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gfa;
                vg convert -p ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gfa > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pg;
                vg convert -fW ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pg > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gfa;
                vg index -t 16 -L -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.xg ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pg;
                vg gbwt -x ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.xg -o ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gbwt -P --pass-paths;
                vg prune -u -g ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gbwt -k 31 -m ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.node_mapping ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pg > ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pruned.vg;
                vg index -g ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gcsa -f ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.node_mapping ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pruned.vg;
                rm ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.og; rm ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.node_mapping; rm ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pruned.vg;
            fi
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
                vg filter -r 0.0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam;
                cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam;
                complex_gene=true;
                #if [ $(vg paths -Lv  ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep -v "${gene}\*" | grep "IMGT\|IGv2\|OGRDB" | wc -l) -gt 0 ]; then
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ]; then
                    echo "Complex locus detected for ${gene} - with multiple other genes in overlapping ASC clusters - filtering out reads aligning to non-gene nodes";
                    sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes;
                    Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes;
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg;
                    vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam;
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt;
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq;
                    # obvious off-target reads removed - pairwise similarity filtering to succinct locus graph now
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.prefilt.gam;
                    vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gam;
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gam > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq;
                    # redo read-based filtering using succinct version of the graph
                    secondary_filtering=false;
                    if [ "${secondary_filtering}" = true ]; then
                        vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg -g ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                        vg filter -r 0.0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                        vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
                        rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes;
                        rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt;
                        vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa;
                        sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                        Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                        vg find -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg
                        vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                        vg find -x ${outdir}/${sample_id}.${graph}.${gene}.succinct_pancluster.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt
                        vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
                    fi
                    # realign to locus-specific + succinct version of the graph
                    if [[ $(echo $gene | grep "IGHV4-4$\|IGHV4-61$\|IGHV4-59$" | wc -l) -ge 1 ]]; then
                        cp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
                        vg convert -fW ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa;
                        cp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
                        cp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt
                        cp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa
                        cp ${genotyping_nodes_dir}/gene_graphs/wg_immunovar.IGHV44-459-461_custom.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp
                    else
                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
                        vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa;
                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt
                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa
                        cp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp
                    fi
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                    vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                else
                    echo "Complex locus detected for ${gene} - multiple ASC clusters of target gene - retaining all reads"
                    vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                fi
            elif [ $(vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep "IMGT\|IGv2\|OGRDB" | grep -v "${gene}\*\|${gene}_\|${gene}" | wc -l) -gt 0 ]; then
                echo "Complex locus detected for ${gene} - single ASC cluster but additional genes present - filtering out reads aligning to non-gene nodes"
                sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                # add all orphon nodes to the list - if gene isn't an orphon gene
                if [ $(echo ${gene} | grep -v "/" | wc -l) -gt 0 ] && [ $(grep "/" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes | wc -l) -gt 0 ]; then
                    echo "Gene is not an orphon gene - and orphon alleles present - adding orphon nodes to filtering list";
                    complex_gene=true
                    grep "/" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes | cut -f3 | sed s/"\\+\|-"//g | tr ',' '\n' | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.orphon.filteringnodes;
                    cat ${outdir}/${sample_id}.${graph}.${gene}.orphon.filteringnodes ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes;
                    # ensure relevant gene-specific nodes are not used to filter out reads (if orphon alleles in same component)
                    grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes | grep ${gene} - | cut -f3 | sed s/"\\+\|-"//g | tr ',' '\n' | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.nonorphon.filteringnodes;
                    grep -Fvx -f ${outdir}/${sample_id}.${graph}.${gene}.nonorphon.filteringnodes ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes;
                fi
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes ]; then
                    echo "Filtering out reads aligning to non-gene nodes"
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg
                    vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                    if [ $(vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep "IMGT\|IGv2\|OGRDB" | wc -l) -gt 20 ]; then
                        vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam;
                    else
                        vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam;
                    fi
                else
                    vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                fi
            else
                echo "Non-ASC-based analysis for ${gene}"
                vg filter -r 0.0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            fi
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            echo "1) Performing inference on haplotype graph labled only with alleles of interest";
            # if variable gene -- limit to OGRDB validated set of alleles?
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "IMGT\|IGv2\|OGRDB" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
            valid_alleles=true
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
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg;
            # Avoid parse_graph_vgflow.py  script completely - bug compressing duplicate nodes of the same sequence content, and not needed - or update script to just use tmp1.gfa
            ############################
            if [ "${complex_gene}" = true ]; then
                echo "Complex gene: performing sequence to valid allele-specific graph-based filtering"
                ############################
                allele_graph_filt=false
                if [ "${allele_graph_filt}" = true ]; then
                    vg paths -Fx ${outdir}/${sample_id}.${graph}.${gene}.vg | seqkit seq -g '-' | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta;
                    if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta  | wc -l ) -gt 2 ]; then
                        vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta -N -a | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa
                        odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.alleles.gfa;
                        rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa;
                    else
                        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta > ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.paf
                        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.alleles.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.paf -g ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.alleles.gfa;
                        rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.rough.gfa;
                    fi
                    vg convert -g -p ${outdir}/${sample_id}.${graph}.${gene}.alleles.gfa > ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg;
                    # augment ids before mapping - otherwise gaf is off by 1!
        #            vg ids -i -1 ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg > ${outdir}/${sample_id}.${graph}.${gene}.alleles.base0.pg;
        #            mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.base0.pg ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg;
                    vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg;
                    vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg -o ${outdir}/${sample_id}.${graph}.${gene}.alleles.gbwt -P --pass-paths
                    vg index -g ${outdir}/${sample_id}.${graph}.${gene}.alleles.gcsa ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg
                    vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg -g ${outdir}/${sample_id}.${graph}.${gene}.alleles.gcsa -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.alleles.gam
                    vg filter -r 0 -P -s 1 -q 0 -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.alleles.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam;
                    vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam;
                    vg filter -r 0 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam
                else
                    echo "Complex gene - sequence to graph alignment filtering (90% pairwise identity) - for allele graph-based filtering set allele_graph_filt=true"
                    vg filter -r 0.9 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam
                fi
                #vg ids -i -1 ${outdir}/${sample_id}.${graph}.${gene}.alleles.pg | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa # ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa #
            else
                echo "sequence to graph alignment-based filtering (95% pairwise identity to local haplotype graph)"
                vg filter -r 0.95 -P -q 0 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam
            fi
            ############################################################################################################################
            #vg map --gaf -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gaf
            #vg convert -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.strict.gam ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gaf
            # or 
#            vg filter -r 0.95 -P -s 1 -q 5 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam -v | vg convert -G - ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gaf
            vg filter -r 0.9 -P -s 1 -q 0 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            ############################################################################################################################
            # gene-specific read depth for minimum strain-level coverage
            depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
            # min_strain_depth=$(bc -l <<< "scale=3;${depth_locus}/$gene_min_len")
            min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.05"| awk '{printf("%d\n",$1 + 0.5)}')
#            min_strain_depth=$(bc -l <<< "scale=2;(${depth_locus}*0.05)*$(echo "${asc_cluster[@]}" | wc -l)"| awk '{printf("%d\n",$1 + 0.5)}')
            if (( $(awk 'BEGIN {print ("'"$min_strain_depth"'" < 0.1) ? "1" : "0"}') )); then
                min_strain_depth=0.1
            fi
            echo "Minimum strain depth required: $min_strain_depth";
            cd ${outdir};
            # allele inference:
            mkdir -p ${outdir}/${gene}_allele_inference;
            #vg ids -i -1 ${outdir}/${sample_id}.${graph}.${gene}.vg | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa #need to increase by -1, cant decrease by 1 lol
            #vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json;
            #vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa;
            #python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0;
            #vg paths -Ev ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa > vgflow_paths.txt;
            #orig_len=$(cut -f2 vgflow_paths.txt | sort | uniq);
            #vg paths -Ev ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa > vgflow_paths.final.txt;
            #parsed_len=$(cut -f2 vgflow_paths.final.txt | sort | uniq);
            #if [ $(printf "%d\n" ${orig_len} | sort -n | head -1) -gt $(printf "%d\n" ${parsed_len} | sort -n | head -1) ]; then
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
                echo "Compress graph bug - allele truncated - using alternative graph prep"
                vg filter -r 0.9 -P -s 1 -q 0 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam -v | vg convert -G - ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gaf
                vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
                gafpack --gfa ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa --gaf ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gaf -lc | grep -v "#" | awk '{print NR-1 ":" $0}' > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt
            #fi
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
            if [ -s ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta ];then
                echo "parsing for association testing";
            else
                echo "No alleles > depth threshold"
                assoc_testing=false
            fi
            if [ "${assoc_testing}" = true ]; then
                echo "2) Augmenting annotated post-flow inference graph with reads for association testing";
                sed s/' '/_/g ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta > ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta
                # should ensure the full-length alleles are used moving forward - not extended/truncated - re-extract from haplotype graph
                #...
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
                            grep ${chr_loc%:*} ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                            seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.adding.fasta > ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta
                            seqkit grep -p $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta | cat - ${outdir}/${sample_id}.${graph}.${gene}.$chr_graph_tmp.fasta | seqkit sort --quiet - > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta
                            if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta  | wc -l ) -gt 2 ]; then
                                echo "Progressive graph"
                                vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.rough.gfa;
                                odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa
                                rm ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.rough.gfa;
                                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa
                            else
                                echo "All-to-all graph"
                                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf
                                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp; 
                                mv ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.genome_graph_ref.gfa
                            fi
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
                        cat ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta <(seqkit grep -r -p $ref_backbone ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta) > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                        echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local CHM13 reference sequence"
                        sed s/' '/_/g -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                        seqkit sort --quiet ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                        if [ $(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta  | wc -l ) -gt 2 ]; then
                            vg msga -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -w 500 | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.rough.gfa;
                            odgi sort -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.rough.gfa --threads 16 -p bgs -O -o - -P | odgi chop -i - -c 32 -o - | odgi view -i - -g > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa;
                            rm ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.rough.gfa;
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.rough.gfa;
                        else
                            minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
                            seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                            gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
                        fi
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
#                    vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg -g ${outdir}/${sample_id}.${graph}.${gene}.alleles.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.alleles.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
#                    vg filter -r 0.5 -P -s 1 -q 60 -x ${outdir}/${sample_id}.${graph}.${gene}.alleles.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
#                   re-align reads for variant calling
#                    vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
                    vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
                    vg filter -r 0.95 -P -s 1 -q 0 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
                    vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth;
                    depth_aug=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth)
                    aug_depth=$(bc -l <<< "scale=2;${depth_aug}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
                    if [ "${aug_depth}" -gt 3 ]; then
                        augment_cov=${aug_depth}
                    else
                        augment_cov=3
                    fi
                    echo "Minimum coverage to add breakpoint: ${augment_cov} (3 <--> 10% of strain depth)"
                    vg augment -m ${augment_cov} -q 5 -Q 20 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg;
                    vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
                    #
                    vg index -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.xg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg;
                    vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gbwt -P --pass-paths
                    #vg gbwt -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gbz --gbz-format -P --pass-paths;
                    vg prune -u -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pruned.vg
                    vg index -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pruned.vg
                #   retain only 100% identity aligned reads
                    vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented2.gam
                    vg filter -r 1 -P -s 1 -q 0 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented2.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.exact.gam
                    # gaf of reads aligned to augmented graph - gafpack to get coverage/depth
                    #vg convert -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf
                    vg convert -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.exact.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf
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
                    gafpack --gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa --gaf ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf -l > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage
                    #Rscript ${bigfoot_dir}/augment_graph_wdepth.R ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage >> ${outdir}/${sample_id}.putative_variants.csv
                    Rscript ${bigfoot_dir}/augment_graph_wdepth_asc.R ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage >> ${outdir}/${sample_id}.putative_variants.csv;
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa;
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa;
                    vg mod -N -X 256 -n -U10 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.tmp.gfa;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.tmp.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa;
                    vg mod -N -X 256 -n -U10 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa | vg convert -fW - > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.tmp.gfa;
                    mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.tmp.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa;
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths;
                    grep "chm\|grch\|${gene}" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths;
                    vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.pg;
                    vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa;
                    grep "^P" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.plines;
                    grep "^P" -v ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.noplines;
                    cat ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.noplines ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.plines > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa;
                    # potentially re-align reads and then pack/call variants so we can incorporate depth info
                    vg deconstruct -a -n -L 1 -P "chm" -P "grch" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.coding.vcf;
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
                    if [ "${de_novo}" = true ]; then
                        echo "De-novo allele inference...";
                        depth_graph=${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa;
                        head -1 ${depth_graph} > ${depth_graph%.gfa}.header;
                        grep "^P" ${depth_graph} > ${depth_graph%.gfa}.plines;
                        grep "^S" ${depth_graph} > ${depth_graph%.gfa}.slines;
                        grep "^L" ${depth_graph} > ${depth_graph%.gfa}.llines;
#                        Rscript ${bigfoot_dir}/phase_locus_depth_graph.R ${depth_graph%.gfa}.plines;
                        Rscript ${bigfoot_dir}/phase_locus_depth_graph_asc.R ${depth_graph%.gfa}.plines;
                        grep -P "\tIG|\tOG|\tTR|\tgrch|\tchm" ${depth_graph%.gfa}.cleaned_paths > ${depth_graph%.gfa}.noreads.cleaned_paths;
                        cat ${depth_graph%.gfa}.header ${depth_graph%.gfa}.slines ${depth_graph%.gfa}.llines ${depth_graph%.gfa}.cleaned_paths > ${depth_graph%.gfa}_cleaned.gfa;
                        cat ${depth_graph%.gfa}.header ${depth_graph%.gfa}.slines ${depth_graph%.gfa}.llines ${depth_graph%.gfa}.noreads.cleaned_paths > ${depth_graph%.gfa}_cleaned.filt.gfa;
                        rm ${depth_graph%.gfa}.*cleaned_paths;
                        rm ${depth_graph%.gfa}.*lines;
                        rm ${depth_graph%.gfa}.*header;
                        # map - remove nodes with 0/low coverage - unitigs - stitch together using contiguous abundance...
                        # at most number of alleles = ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta + number of nodes without a path
#                            vg view -X ${gene}.careful.w_ref_paths.augmented.gam | seqkit fq2fa - > ${gene}.careful.w_ref_paths.augmented.fasta
#                            odgi unitig -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -t 151 > ${gene}.careful.unitigs.fa 
#                            Bifrost build -r ${gene}.careful.unitigs.fa -k 31 -o ${gene}.careful.bifrost
#                            gunzip ${gene}.careful.bifrost*gfa.gz
#                           ~/tools/SPAdes-4.0.0-Linux/bin/pathracer ${gene}.careful.w_ref_paths.augmented.fasta ${gene}.careful.bifrost.gfa --nt --annotate-graph --output ${gene}.careful.pathracer
                        #
                        # re-align to *cleaned.filt.gfa --> use the minimum length-adjusted node abundance value for each inferred allele to obtain a ~confidence~ score for the allele
                        variant_graph=${outdir}/${sample_id}.${graph}.${gene}.variant;
                        vg mod -X 32 ${depth_graph%.gfa}_cleaned.filt.gfa | vg convert -fW - > ${variant_graph}.gfa;
                        # creating a new graph from scratch - limited to inferred alleles
                        #vg paths -Fx ${depth_graph%.gfa}_cleaned.filt.gfa | seqkit seq -g '-' | seqkit grep -nrv -p "grch|chm" > ${variant_graph}.fasta;
                        #if [ $(grep ">" ${variant_graph}.fasta | wc -l ) -gt 2 ]; then
                        #    vg msga -f ${variant_graph}.fasta -N -a | vg convert -fW - > ${variant_graph}.rough.gfa;
                        #    odgi sort -i ${variant_graph}.rough.gfa --threads 16 -p Ygs -O -o - -P | odgi chop -i - -c 3 -o - | odgi view -i - -g > ${variant_graph}.rough2.gfa;
                        #    gfaffix ${variant_graph}.rough2.gfa -o ${variant_graph}.gfa;
                        #    rm ${variant_graph}.rough.gfa; rm ${variant_graph}.rough2.gfa;
                        #    vg mod -X 3 ${variant_graph}.gfa | vg convert -fW - > ${variant_graph}2.gfa;
                        #    mv ${variant_graph}2.gfa ${variant_graph}.gfa;
                        #else
                        #    minimap2 -x asm20 -t 16 -c -X ${variant_graph}.fasta ${variant_graph}.fasta > ${variant_graph}.paf;
                        #    seqwish -s ${variant_graph}.fasta -p ${variant_graph}.paf -g ${variant_graph}.rough.gfa -b ${outdir}/seqwish_${sample_id}.${graph};
                        #    odgi sort -i ${variant_graph}.rough.gfa --threads 16 -p Ygs -O -o - -P | odgi chop -i - -c 3 -o - | odgi view -i - -g > ${variant_graph}.rough2.gfa;
                        #    gfaffix ${variant_graph}.rough2.gfa -o ${variant_graph}.gfa;
                        #    vg mod -X 3 ${variant_graph}.gfa | vg convert -fW - > ${variant_graph}2.gfa;
                        #    mv ${variant_graph}2.gfa ${variant_graph}.gfa;
                        #fi
                        vg convert -g -p ${variant_graph}.gfa > ${variant_graph}.pg;
                        vg index -t 16 -L -x ${variant_graph}.xg ${variant_graph}.pg;
                        vg gbwt -x ${variant_graph}.xg -o ${variant_graph}.gbwt -P --pass-paths;
                        #vg prune -u -g ${variant_graph}.gbwt -k 31 -m ${variant_graph}.node_mapping ${variant_graph}.pg > ${variant_graph}.pruned.vg
                        #vg index -g ${variant_graph}.gcsa -X 2 -f ${variant_graph}.node_mapping ${variant_graph}.pruned.vg
                        vg prune ${variant_graph}.pg > ${variant_graph}.pruned.vg
                        vg index -g ${variant_graph}.gcsa -X 4 ${variant_graph}.pruned.vg
                        rm ${variant_graph}.pruned.vg;
                        rm ${variant_graph}.node_mapping;
                    #   #retain only 100% identity aligned reads
                        # this one likely the best - low segment abundance for decent proportion of allele 
                        #vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.exact.gam -x ${variant_graph}.xg -g ${variant_graph}.gcsa -1 ${variant_graph}.gbwt -t 4 -M 1 > ${variant_graph}.gam;
                        # see if this blows up the coverage estimates - should still be low for 3-30, but maybe not 100% of the gene?
                        vg map -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -x ${variant_graph}.xg -g ${variant_graph}.gcsa -1 ${variant_graph}.gbwt -t 4 -M 1 > ${variant_graph}.gam
                        vg filter -r 1 -P -s 1 -x ${variant_graph}.xg -D 0 -fu -t 4 ${variant_graph}.gam -v > ${variant_graph}.filt.gam;
                        #if [ -s ${variant_graph}.filt.gam ];then
                        #    vg convert -G ${variant_graph}.filt.gam ${variant_graph}.xg > ${variant_graph}.gaf;
                        #    vg ids -i -1 ${variant_graph}.pg | vg convert -fW - > ${variant_graph}.gafpack.gfa;
                        #fi
                        #gafpack -l -g ${variant_graph}.gafpack.gfa -a ${variant_graph}.gaf
                        vg pack -x ${variant_graph}.pg -g  ${variant_graph}.filt.gam -e -o ${variant_graph}.pack;
                        vg pack -x ${variant_graph}.pg -g  ${variant_graph}.filt.gam -e -d > ${variant_graph}.pack.tsv;
                        # use coverage.tsv to infer minimum depth at each allelic position
                        vg depth -d -c -k ${variant_graph}.pack ${variant_graph}.pg -m 0 -b 1 | grep -v "grch\|chm" > ${variant_graph}.coverage.tsv;
                        echo "Coverage file: $(wc -l ${variant_graph}.coverage.tsv)"
                    else
                        echo "No novel allele inference requested - to perform novel allele inference set de_novo=true"
                    fi
                else
                    echo "Insuffient coverage to infer alleles for ${gene} - not trying to embed variation"
                fi
            else
                echo "No association testing performed for ${gene} - set assoc_testing=true to embed reads/prep locus graph for association testing"
            fi
            # remove files we dont need anymore
            ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v "${gene}.genome_graph_ref.augmented.*gfa\|coverage.tsv\|pack.tsv\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|.vcf\|node_abundance\|depth" > ${outdir}/${sample_id}_${gene}_files.txt
            #ls ${outdir}/${gene}\.*  | grep "haps.fasta\|alleles" >> ${outdir}/${sample_id}_${gene}_files.txt
            ls ${outdir}/*${gene}* | grep "asc_" >> ${outdir}/${sample_id}_${gene}_files.txt;
            #ls ${outdir}/*${gene}.alleles*fasta >> ${outdir}/${sample_id}_${gene}_files.txt;
            xargs rm < ${outdir}/${sample_id}_${gene}_files.txt;
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

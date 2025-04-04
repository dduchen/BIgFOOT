#!/bin/bash

cd ${datadir}
gene=${each%.nodes.txt}
gene=${gene%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
# ASC-based clustering only for IG/TR genes - loci=IMGT
loci=IMGT
echo "$gene --> $loci"
#
unset asc_cluster #ensure no carryover of variables from previous gene analysis
each=ig_asc/${each}
if [[ "$loci" =~ ^(IMGT)$ ]]; then
    echo "IG/TR inference"
    outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference
    mkdir -p ${outdir}/seqwish_${sample_id}.${graph}
    # ASC cluster - if more than 1 gene included in ASC cluster - use clustering-based gene ID for graph construction:
    grep "${gene}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv > ${outdir}/potential_asc_for_${gene}
    asc_cluster=$(cut -f1 ${outdir}/potential_asc_for_${gene} | sed s/'\*.*'//g | sed s/.*#1#//g | sort | uniq)
    if [ $(echo "${asc_cluster[@]}" | wc -l) -gt 1 ]; then
        echo "Complex locus - multiple V genes in this cluster";
        complex_gene=true
    fi
else
    echo "ERROR: perform gene-based analysis for ${gene}"
fi
# check if graph exists in subdir - if it does use it rather than reconstruct things
# check if we've completed analysis for this gene already
if [ -s ${outdir}/${sample_id}_${gene}_files.txt ]; then
    echo "Analysis completed for ${gene} - did you restart or fail to cleanup? (File=${outdir}/${sample_id}_${gene}_files.txt)"
else
    if [[ "$loci" =~ ^(IMGT)$ ]]; then
        echo "IMGT: Immunoglobulin/T-cell receptor gene inference"
        if [ -s ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gcsa ]; then
            echo "Using stored locus specific graph/indexes";
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.genotyping.immune_subset.vg ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gcsa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gbwt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.pg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
            cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
            vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            if [ -s ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.gbwt ];then
                cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.gbwt ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gbwt
                cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.xg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.xg
                cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.pg ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.pg
                cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.gcsa ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa
                cp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.succinct_locus.gcsa.lcp ${outdir}/${sample_id}.${graph}.${gene}.succinct_locus.gcsa.lcp
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
            seqkit grep -r -p "\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.match.fasta
            seqkit grep -r -v -p "\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.unmatch.fasta
            cat ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.match.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.unmatch.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            seqkit grep -r -p "\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.alleles.fasta
            # remove multi-mapping alleles since we're grabbing from a range of ASC-informed loci
            seqkit grep -r -p "#0$" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.fasta.tmp && mv ${outdir}/${gene}.alleles.fasta.tmp ${outdir}/${gene}.alleles.fasta
            seqkit grep -r -p "/" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.offtarget.fasta
            # remove co-localized alleles from genes outside the ASC cluster of interest:
            cut -f1 ${outdir}/potential_asc_for_${gene} | sed s/'\*.*'/"\\\\*"/g | sort | uniq > ${outdir}/potential_asc_for_${gene}_ascs.refined.txt
            seqkit grep -r -f ${outdir}/potential_asc_for_${gene}_ascs.refined.txt ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.fasta.tmp && mv ${outdir}/${gene}.alleles.fasta.tmp ${outdir}/${gene}.alleles.fasta
            if [ -s ${outdir}/${gene}.alleles.offtarget.fasta ]; then
                echo "Orphon alleles in this ${gene} - removing any allele with dist=0 to orphon alleles"
                cp ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.alleles.order.fasta
                #Rscript ${bigfoot_dir}/Remove_dup_allele.R ${outdir}/${gene}.alleles.order.fasta
                Rscript ${bigfoot_dir}/Remove_dup_allele.R ${outdir}/${gene}.alleles.order.fasta
            fi
            seqkit grep -r -v -p "/" ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.exact.fasta
            seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
            # set min length of haplotypes = 100, max length should scale with allele length --> then downsample haplotypes
            seqkit stats ${outdir}/${gene}.alleles.fasta > ${outdir}/${gene}.alleles.stats
            gene_min_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f6 -d' ' | sed s/","//g)
            gene_max_len=$(sed -n 2p ${outdir}/${gene}.alleles.stats | tr -s ' ' | cut -f8 -d' ')
            ##############################
            # use local haplotypes #
            ##############################
            seqkit grep -v -r -p "\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.haps.fasta
            if [ -s ${outdir}/${gene}.haps.fasta ]; then
                echo "Local haplotypes found for: ${gene_actual}"
#                seqkit seq --min-len $(bc -l <<< "scale=2;${gene_min_len}*0.9"| awk '{printf("%d\n",$1 + 0.5)}') --max-len 15000  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                seqkit seq --min-len $(bc -l <<< "scale=2;${gene_min_len}*0.9"| awk '{printf("%d\n",$1 + 0.5)}') --max-len $(bc -l <<< "scale=2;${gene_max_len}*3"| awk '{printf("%d\n",$1 + 0.5)}')  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
                # retain haplotypes with perfect match to one of our alleles
                mkdir -p $outdir/${gene}_haps
                echo "Exact allele:haplotype matching"
                # remove duplicate alleles - arrange orphons first to remove ~real~ alleles with an exact match
                minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
                minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.exact.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.exact.txt
                if [ -s ${outdir}/${gene}.alleles.offtarget.fasta ]; then
                    minimap2 -x sr --secondary=no -c ${outdir}/${gene}.alleles.offtarget.fasta ${outdir}/${gene}.haps.fasta | grep "NM:i:0" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.offtarget.txt
                    cp ${outdir}/${gene}_haps/haps.matching.offtarget.txt ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget.txt
                fi
                cp ${outdir}/${gene}_haps/haps.matching.exact.txt ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt
                echo "$(cut -f1 ${outdir}/${gene}_haps/haps.matching.txt | wc -l) out of $(grep ">" ${outdir}/${gene}.haps.fasta | cut -f1 | wc -l) haplotypes with a $gene (or ASC cluster member) allele mapping"
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
                seqkit grep -r -p "\*" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | minimap2 -x asm20 --secondary=no -c ${outdir}/${gene}.haps.fasta - | cut -f1,6 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt
            fi
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
            grep "${gene}\*" ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv | cut -f1 | sed s/'*'/'\\*'/g > ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt
            # check if >1 chromsome -- (cluster specific) is so, split again
            if [ $(grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f2 | sort | uniq | wc -l) -gt 1 ]; then
                echo "More than 1 chromosome represented in ${gene} -- additionally splitting by chr"
                for chr_loc in $(grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | cut -f2 | sort | uniq); do echo "subgraph for alleles specific to: $chr_loc";
                    chr_graph_tmp=$(echo $chr_loc | sed s/":.*"//g | sed s/".*#"//g)
                    grep $chr_loc ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt - | cut -f1 > ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt
                    grep -f <(sed s/'\*'/'\\*'/g ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt) ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt | cut -f2 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt
                    seqkit grep -r -f <(sed s/'\*'/'\\*'/g ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.patterns.txt) ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    # Alleles might align best to a speciifc haplotype - but need ensure the haplotype belongs in the chromosome-specific subgraph
                    if [ $(grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.chr_assignments.txt | grep ${chr_loc} | wc -l) -ge 1 ]; then
                        seqkit grep -r -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta
                    fi
                    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf
                    seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                    gfaffix ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.${chr_graph_tmp}.gfa
                done
                ls ${outdir}/${sample_id}.${graph}.${gene}.*.gfa | grep -v "vgflow.final" - | grep -v "haplotypes" - > ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt
                odgi squeeze -f ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt -O -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og 
                odgi view -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og -g > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
                #delete temp files
                rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.og;rm ${outdir}/${sample_id}.${graph}.${gene}.graphs_to_squeeze.txt; rm ${outdir}/${sample_id}.${graph}.${gene}.chr*
            else
                echo "Single locus for this cluster of ${gene}"
                seqkit grep -r -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.fasta
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt ]; then
                    grep -f ${outdir}/potential_asc_for_${gene}_ascs.iterative_assignment.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.txt | cut -f2 | sort | uniq > ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt
                    seqkit grep -f ${outdir}/${sample_id}.${graph}.${gene}.haplotype_assignments.patterns.txt ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta >> ${outdir}/${sample_id}.${graph}.${gene}.fasta
                fi
                minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.fasta ${outdir}/${sample_id}.${graph}.${gene}.fasta > ${outdir}/${sample_id}.${graph}.${gene}.paf
                seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.paf -g ${outdir}/${sample_id}.${graph}.${gene}.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
                gfaffix ${outdir}/${sample_id}.${graph}.${gene}.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
            fi
#
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
            cp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.genotyping.immune_subset.vg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.xg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gcsa
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa.lcp ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gcsa.lcp
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.gbwt
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes.pg
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haplotypes_ref.fasta
            cp ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.exact.txt ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haps.matching.exact.txt
            cp ${outdir}/${sample_id}.${graph}.${gene}_haps.matching.offtarget.txt ${genotyping_nodes_dir}/ig_asc/gene_graphs/${graph}.${gene}.haps.matching.offtarget.txt
            rm ${outdir}/${gene}.alleles.stats
            rm ${outdir}/${gene}.haps.fasta
        fi
        # parse alignments - probably too permissive currently with ..genotyping.immune_subset.vg for read extraction for complex genes e.g. 3-30
        vg find -x ${graph_base}.xg -l ${workdir}/${gam_file%.gam}.sorted.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam
        vg view -a ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -X | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt
        if [[ $(wc -l <${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam.txt) -ge 1 ]]; then
            echo "Re-aligning reads to locus-specific haplotype graph"
            vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            # copy mapped reads -- use this for augmentation?
            cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
            # extra read filtering for complex genes
#                echo "Complex locus detected for ${gene} - however the multiple ASC clusters contain only taget gene alleles - retaining all reads"
#                vg filter -r 0 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#                fi
# use valid_alleles=true here as well
# Check for orphons / multi-chromosomes, reads mapping to other chromosome we remove
            if [ $(vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa | grep "IMG\|IGv2\|OGR" | grep "/\|__" | wc -l) -gt 0 ]; then
                echo "Complex locus detected for ${gene} - single ASC cluster but additional genes present - filtering out reads aligning to non-gene nodes"
                sed -n '/^#1:/p;/^P/p' ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                Rscript ${bigfoot_dir}/identify_non_gene_nodes_complex_locus.R ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pathnodes
                if [ -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes ]; then
                    echo "Filtering out reads aligning to non-gene nodes"
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -c 0 -N ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filteringnodes > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg
                    vg gamsort ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.gai -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                    vg find -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -l ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -A ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.vg | vg view -X - | seqkit seq -n - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt
                    vg view -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam | seqkit grep -v -n -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.txt - > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq
                    vg map -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.filter.fastq -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
                    vg filter -r 0.0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                else
                    vg filter -r 0.0 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
                fi
            else
                echo "Minimal pairwise sequence-to-graph alignment-based filtering for ${gene}"
                vg filter -r 0 -P -q 5 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
            fi
            vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
            echo "1) Performing inference on haplotype graph labled only with alleles of interest";
            # if variable gene -- limit to OGRDB validated set of alleles?
            vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "IMGT\|IGv2\|OGRDB" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
            if [ "${valid_alleles}" = true ]; then
                echo "Where possible - retaining path labels only for OGRDB validated set of ${gene_actual} alleles";
                if [[ ${gene} == *"IGHV"* ]]; then
                    echo "IGHV";
                    ogrdb_refs="${bigfoot_dir}/../custom_beds/ogrdb_IGH_*fasta"
                    seqkit grep -r -f <(cut -f1 ${outdir}/potential_asc_for_${gene} | sed s/".*#1#"//g | sed s/'\*'/'\\*'/g ) ${ogrdb_refs} | grep ">" | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        if [ -s ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ];then 
                            mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                        else 
                            rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        fi
                    fi
                elif [[ ${gene} == *"IGKV"* ]]; then
                    echo "IGKV";
                    ogrdb_refs=${bigfoot_dir}/../custom_beds/ogrdb_IGK_*fasta
                    seqkit grep -r -f <(cut -f1 ${outdir}/potential_asc_for_${gene} | sed s/".*#1#"//g | sed s/'\*'/'\\*'/g ) ${ogrdb_refs} | grep ">" | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        if [ -s ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ];then 
                            mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                        else 
                            rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        fi
                    fi
                elif [[ ${gene} == *"IGLV"* ]]; then
                    echo "IGLV";
                    ogrdb_refs=${bigfoot_dir}/../custom_beds/ogrdb_IGL_*fasta
                    seqkit grep -r -f <(cut -f1 ${outdir}/potential_asc_for_${gene} | sed s/".*#1#"//g | sed s/'\*'/'\\*'/g ) ${ogrdb_refs} | grep ">" | sed s/'>'//g | sed s/'\*'/'\\*'/g > ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles
                    if [ -s  ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ]; then
                        grep -f ${outdir}/${sample_id}.${graph}.${gene}.ogrdb_alleles ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else
                        grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        if [ -s ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ];then 
                            mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                        else 
                            rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                        fi
                    fi
                else
                    echo "Retaining path labels only for non-orphon genes for: ${asc_cluster[@]}";
                    grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                    if [ -s ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ];then 
                        mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
                    else 
                        rm ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp
                    fi
                fi
            else
                echo "Retaining path labels only for non-orphon ${asc_cluster[@]} If you want to limit variable gene inference to OGRDB set - set 'valid_alleles=true'"
                grep -v "/" ${outdir}/${sample_id}.${graph}.${gene}.alleles > ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.alleles.tmp ${outdir}/${sample_id}.${graph}.${gene}.alleles
            fi
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
            echo "ASC-based allele-specific inference - removing duplicated alleles (e.g., IGHV1-68*02 ~ IGHV1-68D*02)"
            vg paths -F -x ${outdir}/${sample_id}.${graph}.${gene}.vg | seqkit rmdup -s - | grep ">" | sed s/'>'//g > ${outdir}/${sample_id}.${graph}.${gene}.alleles.rmdup
            vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles.rmdup -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
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
                        cat ${outdir}/${sample_id}.${graph}.${gene}.rel.haps.final.annot.fasta <(seqkit grep -r -p $ref_backbone ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta) > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
                        echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local ${ref_backbone} reference sequence"
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
                    gafpack --gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa --gaf ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gaf -l > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage
                    Rscript ${bigfoot_dir}/augment_graph_wdepth_asc.R ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.coverage >> ${outdir}/${sample_id}.putative_variants.csv
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa
                    sed -i 's/^[ \t]*//;s/[ \t]*$//' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa
                    vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.gfa > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths
                    grep "chm\|grch\|*" ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.filt.paths
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
                    if [ "${de_novo}" = true ]; then
                        echo "De-novo allele inference...";
                        depth_graph=${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.depth.gfa;
                        head -1 ${depth_graph} > ${depth_graph%.gfa}.header;
                        grep "^P" ${depth_graph} > ${depth_graph%.gfa}.plines;
                        grep "^S" ${depth_graph} > ${depth_graph%.gfa}.slines;
                        grep "^L" ${depth_graph} > ${depth_graph%.gfa}.llines;
                        Rscript ${bigfoot_dir}/phase_locus_depth_graph_asc.R ${depth_graph%.gfa}.plines;
                        grep -P "\tIG|\tOG|\tTR|\tgrch|\tchm" ${depth_graph%.gfa}.cleaned_paths > ${depth_graph%.gfa}.noreads.cleaned_paths;
                        cat ${depth_graph%.gfa}.header ${depth_graph%.gfa}.slines ${depth_graph%.gfa}.llines ${depth_graph%.gfa}.cleaned_paths > ${depth_graph%.gfa}_cleaned.gfa;
                        cat ${depth_graph%.gfa}.header ${depth_graph%.gfa}.slines ${depth_graph%.gfa}.llines ${depth_graph%.gfa}.noreads.cleaned_paths > ${depth_graph%.gfa}_cleaned.filt.gfa;
                        rm ${depth_graph%.gfa}.*cleaned_paths;
                        rm ${depth_graph%.gfa}.*lines;
                        rm ${depth_graph%.gfa}.*header;
                    else
                        echo "No novel allele inference requested - to perform novel allele inference set de_novo=true"
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
    else
        echo "Unknown locus or incorrect script for ${gene}"
    fi
fi

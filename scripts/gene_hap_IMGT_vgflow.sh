#!/bin/bash

# for slurm-based batch job
#SBATCH --ntasks=1
#SBATCH --job-name=run_bigfoot
#SBATCH --time=8:00:00
#SBATCH --mem=36GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out
#
#
export DATAPATH=${workdir}
export OUTPATH=${workdir}
export datadir=${workdir}
cd ${datadir}
# make $bigfoot_dir python scripts findable
export PYTHONPATH=$PYTHONPATH:$bigfoot_dir
export gam_file=${i}
echo "Performing VG-Flow at gene level for the ${gam_file} alignment"
#
if [[ $graph == "minimal" ]]; then
    echo "Using minimal SV-graph - Deprecated"
#    graphdir=~/project/sv_graphs/sv_graph_minimal
#    graph_base=${graphdir}/wg_sv_immunovar_simple.wDups
    genotyping_nodes_dir=${graphdir}/genotyping_wDups/
elif [[ $graph == "hirh" ]]; then
    echo "Using minimal SV-graph - Deprecated"
#    graphdir=~/project/sv_graphs/sv_graph_minimal
#    graph_base=${graphdir}/wg_sv_immunovar_simple.wDups
    genotyping_nodes_dir=${graphdir}/genotyping_wDups/
elif [[ $graph == "franken" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph - Deprecated"
#    graphdir=~/project/grch38_chm13_immunovar
#    graph_base=${graphdir}/grch38_chm13_immunovar.sort
    genotyping_nodes_dir=${graphdir}/genotyping_nodes/
elif [[ $graph == "wg_immunovar" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG/MHC Haplotypes + IMGT/OGRDB/IPD alleles"
    export graphdir=${bigfoot_source}
    export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar
    export genotyping_nodes_dir=${bigfoot_source}/wg_ig_hla_kir_immunovar_genotyping_nodes/
elif [[ $graph == "ig_hla_kir" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IMGT-IPD alleles - Deprecated"
#    graphdir=~/project/grch38_chm13_immunovar
#    graph_base=${graphdir}/grch38_chm13_immunovar.hla_kir
    genotyping_nodes_dir=${graphdir}/ig_hla_kir_genotyping_nodes/
else
    echo "Define graph"
fi
#
#
# suffix - based on ngmerged output or qc-passed PE reads
if [[ $gam_file == *"ngmerge"* ]]; then
    echo "Alignment of NGMerged fastq data"
    export aln_type="ngmerged";
else echo "Alignment of non-merged reads";
    export aln_type="pe"
fi
# 
# sort and index only once...
if [ -s ${gam_file%.gam}.sorted.gam ]; then
    echo "Sorted gam file for ${sample_id} already exists - using it (${gam_file%.gam}.sorted.gam)"
    if [ -s "${gam_file%.gam}.sorted.gam.gai" ]; then
        echo "index exists"
    else
        vg index -l ${gam_file%.gam}.sorted.gam
    fi
else
    vg gamsort ${gam_file} -i ${gam_file%.gam}.sorted.gam.gai -p > ${gam_file%.gam}.sorted.gam
fi
# run HLA / KIR inference -- only uncomment this line in IMGT script
#sbatch --export=i=$i,sample_id=$sample_id,graph=$graph,aln_type='pe',workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} ${bigfoot_dir}/gene_hap_HLA_vgflow.sh
#sbatch --export=i=$i,sample_id=$sample_id,graph=$graph,aln_type='pe',workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} ${bigfoot_dir}/gene_hap_KIR_vgflow.sh
#
#-- global initial pairwise percent identity filtering: 0%
export perc_id_filter="0.00"; echo "Using a final pairwise filter of $perc_id_filter for reads aligned to locus";
#
export outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference
mkdir -p ${outdir}
if [ "${reprocess}" = true ]; then
    echo "Forcing reprocessing - deleting gene txt files"
    if [ -s $(ls ${outdir}/*files.txt | head -1) ]; then
        rm $outdir/*files.txt
    fi
fi
echo "Storing output here: ${outdir}"
mkdir -p ${outdir}/HLA
#mkdir -p ${outdir}/KIR
# note: coeliac succesptibility genes - DQA1/DQB1/HLA-C/IGHV gene
# prefetch locus-specific graphs from following directory if possible:
mkdir -p ${genotyping_nodes_dir}/gene_graphs/
mkdir -p ${genotyping_nodes_dir}/ig_asc/gene_graphs/
# IG/TR inference
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^IGH\|^IGLV\|^DQA1\|^DQB1\|^C\\." | grep -v "__\|IGHD\|IGHJ");do echo $each;
# potentially save this as script - run in parallel #https://stackoverflow.com/questions/25158583/exporting-the-full-environment-to-gnu-parallel
# parallel!
export assoc_testing=true
# should tie number of parallel jobs to the number of compute nodes + memory
ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^IGH\|^IGLV\|^IGKV\|TR" | grep -v "__\|IGHD\|IGHJ\|TR.J\|TR.D" > ${outdir}/gene_list.txt;
# asc cluster-based approach vs. gene-based approach
ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^IGH\|^IGLV\|^IGKV\|TR" | grep -v "__\|IGHD\|IGHJ\|TR.J\|TR.D" | grep -v "^IGHV\|^IGLV\|^IGKV" > ${outdir}/asc_gene_list.txt;
ls ${genotyping_nodes_dir}/ig_asc/ | grep "nodes.txt" | grep "^IGHV_\|^IGLV_\|^IGKV_" >> ${outdir}/asc_gene_list.txt;
#
export asc_inference=false
export de_novo=true
mkdir -p ${outdir}/tmp
export TMPDIR=${outdir}/tmp
export TEMPDIR=${outdir}/tmp
#
if [ "${asc_inference}" = true ]; then
    parallel -j 6 'export each={}; \
        if [[ ${each} =~ IG[H|L|K]V_ ]]; then
            echo "ASC-allele-focused inference"
            . ${bigfoot_dir}/gene_asc_allele_calling_parallel.sh
        else
            echo "Gene-focused inference"
            . ${bigfoot_dir}/gene_allele_calling_parallel.sh
        fi' :::: <(cat ${outdir}/asc_gene_list.txt );
else
    parallel -j 6 'export each={}; \
        . ${bigfoot_dir}/gene_allele_calling_parallel.sh' :::: <(cat ${outdir}/gene_list.txt ); #IGHE / IGHV2-70D / IGHV3-30-22 / 3-30-5 / IGHV4-30-2
fi
#    . ${bigfoot_dir}/gene_asc_allele_calling_parallel.sh' :::: <(cat ${outdir}/asc_gene_list.txt );
#parallel -j 6 'export each={}; \
#    . ${bigfoot_dir}/gene_allele_calling_parallel.sh' :::: <(cat ${outdir}/gene_list.txt );
#
variant_file=$(ls ${outdir}/*putative_variants.csv)
grep -v ",1_reads\|,2_reads" ${variant_file} > ${variant_file%.csv}_strict.csv
#
rm -rf ${outdir}/seqwish_${sample_id}.${graph}/
#ls ${outdir}/*_files.txt > ${outdir}/${sample_id}_files_rm.txt # keep in dir to avoid re-analyzing samples
#ls -d ${outdir}/seqwish* >> ${outdir}/${sample_id}_files_rm.txt
#xargs rm -rf < ${outdir}/${sample_id}_files_rm.txt
#
# parse results
if [ "${asc_inference}" = true ]; then
    echo -e ${sample_id} "mean" "sd" | sed s/" "/'\t'/g > ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.asc.txt;
    grep ">" ${outdir%haplotype_inference*}"haplotype_inference"/*annot.fasta | sed s/"^.*>"/''/ | sed s/' or '/'_or_'/g > "${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.results_raw.txt";
    if [ "${de_novo}" = true ]; then
        for novel_alleles in $(ls ${outdir}/*_cleaned.filt.gfa); do echo ${novel_alleles};
            vg paths -F -x ${novel_alleles} | seqkit grep -r -p "\*" | grep ">" | sed s/"^.*>"/''/ | sed s/' or '/'_or_'/g >> "${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.results_raw.denovo.txt";
        done
    fi
    for each in $(ls ${outdir%haplotype_inference*}"haplotype_inference"/*.filtered.depth); do echo -e $(echo -ne ${each%.filtered.depth} ' ' | 
        sed s/"__"/"\/"/; echo $(cut -f1,2 $each)) | sed s/" "/'\t'/g | sed s/".*wg_immunovar."//g >> ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.asc.txt;
    done
    # Revert ASC genes to average gene-based depths
    Rscript ${bigfoot_dir}/parse_ASC_clusters_to_genes.R ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.asc.txt ${bigfoot_dir}/../custom_beds/ASC_metadata.matching.tsv
    mv ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.asc.recode.txt ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.txt
else
    echo -e ${sample_id} "mean" "sd" | sed s/" "/'\t'/g > ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.txt;
    grep ">" ${outdir%haplotype_inference*}"haplotype_inference"/*annot.fasta | sed s/"^.*>"/''/ | sed s/' or '/'_or_'/g > "${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.results_raw.txt";
    for each in $(ls ${outdir%haplotype_inference*}"haplotype_inference"/*.filtered.depth); do echo -e $(echo -ne ${each%.filtered.depth} ' ' |
        sed s/"__"/"\/"/; echo $(cut -f1,2 $each)) | sed s/" "/'\t'/g | sed s/".*wg_immunovar."//g >> ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.txt;
    done
    if [ "${de_novo}" = true ]; then
        for novel_alleles in $(ls ${outdir}/*_cleaned.filt.gfa); do echo ${novel_alleles};
            vg paths -F -x ${novel_alleles} | seqkit grep -r -p "\*" | grep ">" | sed s/"^.*>"/''/ | sed s/' or '/'_or_'/g >> "${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.results_raw.denovo.txt";
        done
    fi
fi
# append HLA results if they exist
if [ -s ${outdir%haplotype_inference*}"haplotype_inference"/HLA/*annot.fasta ]; then
    grep ">" ${outdir%haplotype_inference*}"haplotype_inference"/HLA/*annot.fasta | sed s/"^.*>"/''/ >> ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.results_raw.txt;
    for each in $(ls ${outdir%haplotype_inference*}"haplotype_inference/HLA"/*.filtered.depth); do echo -e $(echo -ne ${each%.filtered.depth} ' ' |
        sed s/"__"/"\/"/; echo $(cut -f1,2 $each)) | sed s/" "/'\t'/g | sed s/".*\\/"//g | sed s/".*wg_immunovar."//g >> ${outdir%haplotype_inference*}"haplotype_inference"/../${sample_id}.depth_raw.txt;
    done
else
    echo "No HLA allelic inference for ${sample_id}";
fi
# clean up:
echo "All cleaned up!"
cd ${workdir}
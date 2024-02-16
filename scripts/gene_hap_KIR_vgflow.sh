#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=hap_KIR_vgflow
#SBATCH --time=8:00:00
#SBATCH --mem=36GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out

DATAPATH=$workdir
OUTPATH=$workdir
datadir=$workdir
cd ${datadir}
#
export gam_file="${i}.gam"
echo "Performing VG-Flow at gene level for the ${gam_file} alignment"
#
if [[ $graph == "minimal" ]]; then
    echo "Using minimal SV-graph"
    graph_base="/home/dd392/project/sv_graphs/sv_graph_minimal/wg_sv_immunovar_simple.wDups"
    graphdir="/home/dd392/project/sv_graphs/sv_graph_minimal/"
    genotyping_nodes_dir="/gpfs/ycga/project/kleinstein/dd392/sv_graphs/sv_graph_minimal/genotyping_wDups/"
elif [[ $graph == "hirh" ]]; then
    echo "Using HIRH-embedded SV-graph -- need to update this, resorting to minimal graph"
    graph_base="/home/dd392/project/sv_graphs/sv_graph_minimal/wg_sv_immunovar_simple.wDups"
    graphdir="/home/dd392/project/sv_graphs/sv_graph_minimal/"
    genotyping_nodes_dir="/gpfs/ycga/project/kleinstein/dd392/sv_graphs/sv_graph_minimal/genotyping_wDups/"
elif [[ $graph == "franken" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/grch38_chm13_immunovar.sort"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
    genotyping_nodes_dir="/home/dd392/project/grch38_chm13_immunovar/genotyping_nodes/"
elif [[ $graph == "wg_immunovar" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IPD HLA and KIR alleles"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/whole_genome_ig_hla_kir_immunovar"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
    genotyping_nodes_dir="/home/dd392/project/grch38_chm13_immunovar/wg_ig_hla_kir_immunovar_genotyping_nodes/"
elif [[ $graph == "ig_hla_kir" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph + IPD HLA and KIR alleles"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/grch38_chm13_immunovar.hla_kir"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
    genotyping_nodes_dir="/home/dd392/project/grch38_chm13_immunovar/ig_hla_kir_genotyping_nodes/"
else
    echo "Define graph"
fi
#
curenv=$(declare -p -x)
source /home/${USER}/.bashrc
eval "$curenv"
conda activate /gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env
export 
PYTHONPATH=/gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env/bin/python3.11:/gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env/lib/python3.11/site-packages:$PYTHONPATH
#
# suffix - based on ngmerged output or qc-passed PE reads
if [[ $gam_file == *"ngmerge"* ]]; then
    echo "Alignment of NGMerged fastq data"
    aln_type="ngmerged";
else echo "Alignment of non-merged reads";
    aln_type="pe"
fi
# 
# sort and index only once...
FILE=${gam_file%.gam}.sorted.gam
if test -f "$FILE"; then
    echo "Sorted gam file for ${sample} already exists - using it (${gam_file%.gam}.sorted.gam)"
else
    /home/dd392/tools/vg gamsort ${gam_file} > ${gam_file%.gam}.sorted.gam
    /home/dd392/tools/vg index -l ${gam_file%.gam}.sorted.gam
fi
# run HLA / KIR inference -- only uncomment this line in IMGT script
#sbatch --export=i=$i,sample_id=$sample_id,graph=$graph,aln_type='pe',workdir=${PWD} /home/dd392/project/gene_hap_HLA_vgflow.sh
#sbatch --export=i=$i,sample_id=$sample_id,graph=$graph,aln_type='pe',workdir=${PWD} /home/dd392/project/gene_hap_KIR_vgflow.sh
#
#-- global initial pairwise percent identity filtering: 0%
perc_id_filter="0.00"; echo "Using a final pairwise filter of $perc_id_filter for reads aligned to locus";
#
outdir=${workdir}/${sample_id}_${graph}_genotyping/familywise_${aln_type}_haplotype_inference
mkdir -p ${outdir}
echo "Storing output here: ${outdir}"
mkdir -p ${outdir}/HLA
mkdir -p ${outdir}/KIR
# note: coeliac succesptibility genes - DQA1/DQB1/HLA-C/IGHV gene
# IG/TR inference
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^TR\|^IGH\|^IGLV\|^IGKV" | grep -v "__\|IGHD");do echo $each;
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^IGH\|^IGLV\|^IGKV" | grep -v "__\|IGHD");do echo $each;
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "IGHV1-69" | grep -v "__\|IGHD");do echo $each;
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "IGHV4-4" | grep -v "__\|IGHD");do echo $each;
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^DQB\|^DRB" | grep -v "__\|IGHD" | grep "DRB1");do echo $each;
# KIR inference only:
for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "KIR");do echo $each;
# HLA inference only:
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep -v "TR\|IGH\|IGK\|IGL\|KIR");do echo $each;
# --
cd ${datadir}
sample=${each%.nodes.txt}
gene=${sample%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
#
loci=$(/home/dd392/tools/vg paths -Lv ${graph_base}.xg | grep "#1#${gene}\*" | cut -f1 -d"#" | sort | uniq | grep "IMGT\|HLA\|KIR")
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
else 
    echo "Unknown locus for ${gene}"
fi
#
echo "subsetting graph and alignment to region surrounding locus of interest: $gene"
/home/dd392/tools/vg find -x ${graph_base}.xg -c 150 -L -N ${genotyping_nodes_dir}/${each} > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
echo "Extracting local haplotypes reads --> ${sample_id}.${graph}.${gene}.haplotypes.gfa"
/home/dd392/tools/vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "grch\|chm" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt
/home/dd392/tools/vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt | 
/home/dd392/tools/vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg
/home/dd392/tools/vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg -F > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
/gpfs/gibbs/project/kleinstein/dd392/conda_envs/seqkit/bin/seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
/home/dd392/tools/vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "IGL\|hirh_H\|hirh_N\|MHC\|${gene_actual}" > 
${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt
/home/dd392/tools/vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt | 
/home/dd392/tools/vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp && mv 
${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
/home/dd392/tools/vg find -x ${graph_base}.xg -l ${gam_file%.gam}.sorted.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg > 
${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam
#
if [[ "$loci" =~ ^(KIR)$ ]]; then
    echo "KIR gene-based inference"
    time /home/dd392/tools/vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -g 
${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gcsa -1 ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#    time /home/dd392/tools/vg giraffe -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -Z ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbz -H 
${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
    cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
    /home/dd392/tools/vg filter -r 0 -P -s 1 -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam 
-v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
    /home/dd392/tools/vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
    /home/dd392/tools/vg convert -fW ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
    /home/dd392/tools/vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > 
${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
    depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
    min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
    /gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env/bin/python3 /gpfs/ycga/project/kleinstein/dd392/sv_graphs/parse_graph_vgflow.py --sample 
${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
#   /gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env/bin/python3 /home/dd392/tools/vg-flow/scripts/vg-flow_immunovar.py --min_depth 0 --trim 0 --greedy_mode all -m 0 
-c ${min_strain_depth} --remove_included_paths 0 --threads 16 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt 
${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
#   skip some overlapping path-related QC + wfmash instead of minimap
    gene_min_len=$(/home/dd392/tools/vg paths -Ex ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg | grep "KIR" | cut -f2 | sort | uniq | head -1)
    wfmash_param=$(bc -l <<< "scale=2;${gene_min_len}/10" | awk '{printf("%d\n",$1 + 0.5)}')
    if [ "${wfmash_param}" -lt 100 ]; then
        wfmash_param=100
    fi
    time /gpfs/gibbs/project/kleinstein/dd392/conda_envs/vg-flow-env/bin/python3 /home/dd392/tools/vg-flow/scripts/vg-flow_immunovar_long_contigs.py --min_depth 0 --trim 0 
--greedy_mode all --ilp --max_strains 2 -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 --wfmash_param ${wfmash_param} 
${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
    mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
    mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.haps.final.fasta
    mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph.gfa
    rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
    /vast/palmer/apps/avx2/software/R/4.2.0-foss-2020b/bin/Rscript ~/project/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
    echo "Allele-level abundance estimation completed for ${gene} ::"
    grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta
# 2) augment annotated post-flow inference graph with reads for association testing
# limit to chm13 reference...
    /gpfs/gibbs/project/kleinstein/dd392/conda_envs/seqkit/bin/seqkit grep -r -p "chm" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > 
${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta.tmp 
${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
    cat ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
    echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local CHM13 reference sequence"
    sed s/' '/_/g -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
    minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
    /gpfs/gibbs/project/kleinstein/dd392/conda_envs/pggb/bin/seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -p 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
    /home/dd392/tools/gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp; mv 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
    /home/dd392/tools/vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg
    /home/dd392/tools/vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg
    /home/dd392/tools/vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
    /home/dd392/tools/vg index -p -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg;
    /home/dd392/tools/vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -P --pass-paths
    /home/dd392/tools/vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbz --gbz-format -P 
--pass-paths;
    /home/dd392/tools/vg prune -p -u -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -k 31 -m 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
    /home/dd392/tools/vg index -p -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
#   map locus-associated reads to augmented/annotated graph
    /home/dd392/tools/vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -x 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
    /home/dd392/tools/vg filter -r 0 -P -s 1 -q 60 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -D 0 -fu -t 4 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
    /home/dd392/tools/vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth;
    depth_aug=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.depth)
    aug_depth=$(bc -l <<< "scale=2;${depth_aug}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
    if [ "${aug_depth}" -gt 3 ]; then
        augment_cov=${aug_depth}
    else 
        augment_cov=3
    fi
    echo "Minimum coverage to add breakpoint: ${augment_cov} (3 <--> strain depth @ gene)"
    /home/dd392/tools/vg augment -m ${augment_cov} -q 5 -Q 60 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gam > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg;
    /home/dd392/tools/vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
    /home/dd392/tools/vg mod -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp && 
mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg
    /home/dd392/tools/vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.pg > 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
    /home/dd392/tools/gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa -o 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa.tmp 
${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.augmented.gfa
# remove files we dont need anymore
    ls ${outdir}/${sample_id}.${graph}.${gene}\.* | grep -v 
"${gene}.genome_graph_ref.augmented.gfa\|${gene}.genome_graph_ref.gfa\|${gene}.haplotypes.xg\|${gene}.haplotypes.gam\|annot.fasta\|annot.gfa\|final.gfa\|node_abundance\|depth" > 
${outdir}/${sample_id}_${gene}_files.txt
    xargs rm < ${outdir}/${sample_id}_${gene}_files.txt
    rm -rf ${outdir}/seqwish_${sample_id}.${graph}/
else 
    echo "Unknown locus or incorrect script for ${gene}"
fi
#
done
#
ls ${outdir}/*_files.txt > ${outdir}/${sample_id}_files_rm.txt
xargs rm < ${outdir}/${sample_id}_files_rm.txt
# clean up:
echo "All cleaned up!"

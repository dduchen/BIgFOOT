#!/bin/bash

#
# for slurm-based batch job
#SBATCH --ntasks=1
#SBATCH --job-name=prep_bigfoot
#SBATCH --time=8:00:00
#SBATCH --mem=36GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out
#

echo "gam is ${i}"
echo "graph is ${graph}"

DATAPATH=$outdir
OUTPATH=$outdir
cd $OUTPATH
workdir=${PWD}
#
if [[ $graph == "minimal" ]]; then
    echo "Using minimal SV-graph - Deprecated"
#    graphdir=~/project/sv_graphs/sv_graph_minimal
#    graph_base=${graphdir}/wg_sv_immunovar_simple.wDups
elif [[ $graph == "hirh" ]]; then
    echo "Using HIRH-embedded SV-graph - Deprecated"
#    graphdir=~/project/sv_graphs/sv_graph_minimal
#    graph_base=${graphdir}/wg_sv_immunovar_simple.wDups
elif [[ $graph == "franken" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph - Deprecated"
#    graphdir=~/project/grch38_chm13_immunovar
#    graph_base=${graphdir}/grch38_chm13_immunovar.sort
elif [[ $graph == "wg_immunovar" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG/MHC Haplotypes + IMGT/OGRDB/IPD alleles"
    graphdir=${bigfoot_source}
    graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar
elif [[ $graph == "ig_hla_kir" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IMGT-IPD alleles - Deprecated"
#    graphdir=~/project/grch38_chm13_immunovar
#    graph_base=${graphdir}/grch38_chm13_immunovar.hla_kir
else
    echo "Define graph"
fi
#
immune_graph=${graph_base}.subgraph
#
gam_file=${i}
sample=${i%%.*}; sample=${sample##*/}
sample_id=${sample}
if [[ $gam_file == *"sorted.gam" ]]; then
    echo "gam input is already sorted" 
    sorted_gam=${gam_file}
else
    echo "Sorting ${gam_file}" 
    sorted_gam=${gam_file%.gam}.sorted.gam
fi

# suffix - based on ngmerged output or qc-passed PE reads
if [[ $gam_file == *"ngmerge"* ]]; then
    echo "Alignment of NGMerged fastq data"
    aln_type="ngmerged";
else echo "Alignment of non-merged reads";
    aln_type="pe"
fi

sample_immune_graph=${sample}.immune_subset.${graph}.${aln_type}

# sort and index only once...
FILE=${sorted_gam}
if test -f "$FILE"; then
    echo "Sorted gam file for ${sample} already exists - using it (${sorted_gam})"
    if test -f "${sorted_gam}.gai"; then
        echo "index exists"
    else
        vg index -l ${sorted_gam}
    fi
else
    echo "sorting $sample gam"
    vg gamsort ${gam_file} -i ${sorted_gam}.gai -p > ${sorted_gam}
fi

echo "finding ${sample}'s reads which align to nodes of interest - then embedding variation implied by those reads in immunovariation-focused subgraph"
vg find -c 0 -l ${sorted_gam} -x ${graph_base}.xg -A ${immune_graph}.vg > ${sample_immune_graph}.gam
# don't need prep_graph_for_vgflow - performing analysis at gene-level
# uncomment line below to construct whole locus flow graph
#sbatch --export=i=${sample_immune_graph}.gam,outdir=${PWD},graph=$graph /home/dd392/project/prep_graph_for_vgflow.sh;

# depth estimates for whole locus + start gene-level inference
sample_immune_aln=${sample}.${graph}.${aln_type}
immune_graph=${graph_base}.subgraph
cd ${workdir}
#
echo "Embedding sample-specific variation into the whole immune subgraph"
output_graph=${sample_id}_graph.${graph}.${aln_type}
# More conservative filtering at the gene level based on vg map - v liberal thresh here using pairwise identity to locus + graph augmentation using reads with base quality > 5 and breakpoints with coverage > 2
vg filter -r 0.9 -s 1 -fu -t 4 ${sample_immune_graph}.gam -D 0 -v > ${sample_immune_graph}.filtered.gam
vg augment -s -m 2 -q 5 ${immune_graph}.vg ${sample_immune_graph}.filtered.gam > ${output_graph}.vg
echo "normalizing, condensing, sorting augmented graph"
vg mod ${output_graph}.vg -n -U 100 -c -X 256 - | vg ids -c - > ${sample_id}.tmp && mv ${sample_id}.tmp ${output_graph}.vg
vg convert -fW ${output_graph}.vg > ${output_graph}.gfa
odgi build -g ${output_graph}.gfa -OP -o ${output_graph}.og
odgi normalize -i ${output_graph}.og -o ${output_graph}.norm.og -P
#odgi sort -i ${output_graph}.norm.og -p Ygs -P -O -o ${output_graph}.sort.og
odgi view -i ${output_graph}.norm.og -g -a -P >${output_graph}.norm.gfa
gfaffix ${output_graph}.norm.gfa -o ${output_graph}.fix.gfa
#
# cleaning assumes ngmerged route takes sufficiently longer - pe workflow has completed
if [[ $aln_type == *"pe"* ]]; then
# perform some depth estimates:
    vg depth --gam ${sample_immune_graph}.gam ${immune_graph}.xg > ${sample}.${graph}_immune.pe.depth;
    depth_immune=$(awk -F ' ' '{print $1}' ${sample}.${graph}_immune.pe.depth)
    echo "${sample} depth across immunovariation loci: ${depth_immune}"
    # TTN exon 327 depth - immunotyper-SR reference
    #chr2:179424038-179441143
    grch_path_chr2=$(vg paths -Lv ${graph_base}.xg | grep "grch38#chr2:")
    vg depth --gam ${sample_immune_graph}.gam ${graph_base}.xg -p "${grch_path_chr2}:179424038[-179441143]" > ${sample}.${graph}_TTN.pe.depth;
    depth_ttn=$(awk -F ' ' '{print $1}' ${sample}.${graph}_TTN.pe.depth)
    echo "${sample} depth across Titin (TTN) exon 327 (Immunotyper-sr ref): ${depth_ttn}"
# finished!
    echo "fin";
    #sbatch --export=i=${sample_immune_graph}.gam,sample_id=$sample,graph=$graph,aln_type='pe',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/gene_hap_IMGT_vgflow.sh
    export i=${sample_immune_graph}.gam workdir=${PWD} sample_id=${sample} graph=${graph} aln_type='pe' use_augmented=FALSE bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} merged=FALSE graph=wg_immunovar
    . ${bigfoot_dir}/gene_hap_IMGT_vgflow.sh

else echo "fin - performing depth estimates and cleaning up after myself";
# perform some depth estimates:
    vg depth --gam ${sample_immune_graph}.gam ${immune_graph}.xg > ${sample}.${graph}_immune.ngmerged.depth;
    depth_immune=$(awk -F ' ' '{print $1}' ${sample}.${graph}_immune.ngmerged.depth)
    echo "${sample} depth (NGmerged) across immunovariation loci: ${depth_immune}"
    # TTN exon 327 depth - immunotyper-SR reference
    #chr2:179424038-179441143
    grch_path_chr2=$(vg paths -Lv ${graph_base}.xg | grep "grch38#chr2:")
    vg depth --gam ${sample_immune_graph}.gam ${graph_base}.xg -p "${grch_path_chr2}:179424038[-179441143]" > ${sample}.${graph}_TTN.ngmerged.depth;
    depth_ttn=$(awk -F ' ' '{print $1}' ${sample}.${graph}_TTN.ngmerged.depth)
    echo "${sample} depth (NGmerged) across Titin (TTN) exon 327 (Immunotyper-sr ref): ${depth_ttn}"
# finished!
    ls ${sample}.${graph}* | grep -v "_1.fastq\|_2.fastq" | grep -v "filtered\|final\|immune_subset\|sorted.gam\|fix.gfa\|depth" > ${sample}.${graph}_files.txt
    xargs rm < ${sample}.${graph}_files.txt
# run final gene-wise vg-flow inference
#    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='ngmerged',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/familywise_vgflow.sh;
# need to update merged script - but deprecated for now
#    sbatch --export=i=${sample_immune_graph}.gam,sample_id=$sample,graph=$graph,aln_type='ngmerged',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/gene_hap_simple_vgflow.sh;
    export i=${sample_immune_graph}.gam workdir=${PWD} sample_id=$sample graph=$graph aln_type='ngmerged' use_augmented=FALSE bigfoot_dir=${bigfoot_dir} merged=FALSE graph=wg_immunovar
    . ${bigfoot_dir}/gene_hap_IMGT_vgflow.sh
fi
#
rm ${output_graph}.gfa
rm ${output_graph}.vg
rm ${output_graph}.og
rm ${output_graph}.norm*
echo "fin!"

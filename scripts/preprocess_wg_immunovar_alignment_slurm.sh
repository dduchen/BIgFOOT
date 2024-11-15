#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=FqBIgFOOT
#SBATCH --time=24:00:00
#SBATCH --mem=70GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out


# sample=$each # prefix of fastq files  ${each}_1.fastq.gz
# workdir=${PWD} # working directory
# merged=FALSE # pe alignment only
# graph=wg_immunovar # correct graph
#

DATAPATH=${workdir}
OUTPATH=${workdir}
cd ${OUTPATH}

bf_env_load=$(conda info --envs | grep "bigfoot" | awk '{print $NF}')
curenv=$(declare -p -x)
source ${HOME}/.bashrc;
eval "$curenv"
conda activate ${bf_env_load};

PATH=${tools_dir}:$PATH

# ffq approach fastq are bgzip'd and appended correctly
#reformat.sh in1=${sample}_1.fastq in2=${sample}_2.fastq out=${sample}.fastq addslash int=f underscore overwrite=t
#reformat.sh in=${sample}.fastq out=${sample}_#.fastq overwrite=t

if [ -s ${sample}_qc_1.fastq.gz ]; then
    echo "QC'd ${sample} fastq files exist, skipping fastp QC"
else
    echo "Performing QC via fastp - no base correction but minimal read length=36bp"
    fastp -i ${sample}*_*1.f*q.gz -I ${sample}*_*2.f*q.gz -o ${sample}_qc_1.fastq.gz -O ${sample}_qc_2.fastq.gz --detect_adapter_for_pe -w 12 --length_required 36 -h ${sample}_report_fastp.html -R '${sample} fastp Report';
fi
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
echo "Alignment to immunovariation graph-accounting for duplicate IMGT alleles"
#write protected distance file
#cp ${graph_base}.dist ${sample}_${graph}_pe_prep_aln.dist

echo "Aligning QC-passed PE reads"
if [ -s ${sample}.${graph_base##*/}.gam ]; then
    echo "Alignment completed - using: ${sample}.${graph_base##*/}.gam";
else
    time vg giraffe -Z ${graph_base}.gbz -H ${graph_base}.gbwt -m ${graph_base}.min -d ${graph_base}.dist -f ${sample}_qc_1.fastq.gz -f ${sample}_qc_2.fastq.gz -p --watchdog-timeout 60 > ${sample}.${graph_base##*/}.gam
fi
#rm ${sample}_${graph}_pe_prep_aln.dist
# feed completed gam to next stage of workflow
#sbatch --export=i=${sample}.${graph_base##*/}.gam,graph=$graph,outdir=${PWD} /home/dd392/project/filter_immune_subgraph.sh;
export i=${sample}.${graph_base##*/}.gam workdir=${PWD} graph=${graph} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} valid_alleles=${valid_alleles}
. ${bigfoot_dir}/filter_immune_subgraph.sh
#
if [[ $merged == "TRUE" ]]; then
    if [-s ${sample}.NGmerge.fastq.gz ]; then
        echo "${sample} fastq files have been merged, skipping NGmerge + processing steps prior to alignment"
    else
        echo "Merge paired-end reads using overlapping region via NGmerge"
        /home/dd392/tools/NGmerge/NGmerge -1 ${sample}_qc_1.fastq.gz -2 ${sample}_qc_2.fastq.gz -o ${sample}.merged.fastq -f ${sample}.unmerged -p 0.2 -z -n 4 -v -t'/' # 0.2 approximates PEAR parameter;
        echo "Concatenate merged + unmerged reads"
        seqkit seq -t DNA -r -p ${sample}.unmerged_2.fastq.gz -o ${sample}.unmerged_rc_2.fastq.gz;
        cat ${sample}.unmerged_1.fastq.gz ${sample}.unmerged_rc_2.fastq.gz ${sample}.merged.fastq.gz > ${sample}.NGmerge.fastq.gz; # unmerged first to aid with fragment size estimation?
        cat ${sample}.unmerged_1.fastq.gz ${sample}.unmerged_rc_2.fastq.gz > ${sample}.NGmerge_unpaired.fastq.gz;
        echo "final read-length filter (>100bp) and optional fastq header QC (default=no)"
# unmerged reads only
        seqkit seq -m 100 --remove-gaps ${sample}.NGmerge_unpaired.fastq.gz | gzip -c > ${sample}.tmp && mv ${sample}.tmp ${sample}.NGmerge_unpaired.fastq.gz;
# merged reads only
        seqkit seq -m 100 --remove-gaps ${sample}.merged.fastq.gz | gzip -c > ${sample}.tmp && mv ${sample}.tmp ${sample}.merged.edit.fastq.gz;
# combined fastq - merged+unmerged reads
        seqkit seq -m 100 --remove-gaps ${sample}.NGmerge.fastq.gz | gzip -c > ${sample}.tmp && mv ${sample}.tmp ${sample}.NGmerge.fastq.gz;
# uncomment line below for fastq header editing
#        reformat.sh in=${sample}.NGmerge.fastq.gz out=${sample}.tmp.fastq.gz int=f underscore fixheaders=t uniquenames=t overwrite=t; mv ${sample}.tmp.fastq.gz ${sample}.NGmerge.fastq.gz
        echo "Aligning ${sample}'s QC-passed merged + unmerged reads"
        echo "First aligning unmerged reads"
        cp ${graph_base}.dist ${sample}_prep_aln.dist
        time vg giraffe -Z ${graph_base}.gbz -H ${graph_base}.gbwt -m ${graph_base}.min -d ${sample}_prep_aln.dist -f ${sample}.NGmerge_unpaired.fastq.gz -p --watchdog-timeout 60 > ${sample}.unmerged.gam
# merged reads
        echo "Aligning merged reads"
        cp ${graph_base}.dist ${sample}_prep_aln.dist
        time vg giraffe -Z ${graph_base}.gbz -H ${graph_base}.gbwt -m ${graph_base}.min -d ${sample}_prep_aln.dist -f ${sample}.merged.edit.fastq.gz -p --watchdog-timeout 60 > ${sample}.merged.gam
# combine into single gam
        echo "Concatenating NGmerged + unmerged reads"
        cat ${sample}.unmerged.gam ${sample}.merged.gam > ${sample}.${graph_base##*/}.ngmerge.gam
        rm ${sample}_prep_aln.dist;
        rm ${sample}.merged*;
        rm ${sample}.NGmerge*;
        rm ${sample}.unmerged*;
        echo "feeding completed ngmerge-based gam to next stage of workflow - immunovariation filtering"
#        sbatch --export=i=${sample}.${graph_base##*/}.ngmerge.gam,graph=$graph,outdir=${PWD} /home/dd392/project/filter_immune_subgraph.sh;
        export i=${sample}.${graph_base##*/}.ngmerge.gam workdir=${PWD} graph=$graph vg_flow_dir=${vg_flow_dir} tools_dir=${tools_dir}
        . ${vg_flow_dir}/filter_immune_subgraph.sh
    fi
# carry on with rest of pipeline
# Additional optional editing of fastq headers
#    reformat.sh in1=${sample}.NGmerge.fastq.gz out=${sample}.NGmerge.fastq2.gz uniquenames=t overwrite=t
vg validate ${graph}.xg -a ${sample}.${graph_base##*/}.gam
#
else
    echo "If interested in merging PE reads prior to inference, set 'merged=TRUE'"
fi
echo "fin!"

#!/bin/bash

# for slurm-based batch job
#SBATCH --ntasks=1
#SBATCH --job-name=process_bam
#SBATCH --time=8:00:00
#SBATCH --mem=36GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out
#
#
DATAPATH=$workdir
OUTPATH=$workdir
datadir=$workdir
cd ${datadir}
workdir=${PWD}
# make ${bigfoot_dir} python scripts findable
export PYTHONPATH=$PYTHONPATH:$bigfoot_dir
echo "Processing BAM: ${bam_file}"
#
if [[ $graph == "minimal" ]]; then
    echo "Using minimal SV-graph - Deprecated"
elif [[ $graph == "hirh" ]]; then
    echo "Using minimal SV-graph - Deprecated"
elif [[ $graph == "franken" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph - Deprecated"
elif [[ $graph == "wg_immunovar" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG/MHC Haplotypes + IMGT/OGRDB/IPD alleles"
    graphdir=${bigfoot_source}
    graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar
    graph=${graph_base}".subgraph"
elif [[ $graph == "ig_hla_kir" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IMGT-IPD alleles - Deprecated"
else
    echo "Define graph"
fi
#
if [[ $ref_build == "grch38" ]]; then
    ref=${ref};
    immunovar_bed="${bigfoot_dir}/../custom_beds/grch38_custom_immunovar_coords.bed";
elif [[ $ref_build == "chm13" ]]; then
    ref="${ref}";
    immunovar_bed="${bigfoot_dir}/../custom_beds/chm13_custom_immunovar_coords.bed";
else
    echo "Currently BIgFOOT only supports GRCh38-based BAM input. Future release will have CHM13 as an option"
fi
echo "Extracting candidate reads from BAM"
input_aln=${bam_file##*/};
sample=${input_aln%.bam};
# check if BAM or CRAM
aln_linear=$(echo ${input_aln} | sed s/.*\\.//g)
FILE=${input_aln%.${aln_linear}}.unmapped.fastq.gz
if test -f "$FILE"; then
    echo "${i##*/} has been processed already - using it"
else
    echo "Processing: ${input_aln##*/}"
    if [[ ${bam_file} == *"ftp"* ]]; then
        wget -c ${bam_file};
        samtools index ${input_aln}
    else
        echo "working with a BAM/CRAM file";
        if test -f "${input_aln%.bam}*i"; then
            echo "${input_aln} already indexed";
        else
            echo "indexing ${input_aln}";
            samtools index ${input_aln};
        fi
    fi
    echo "Extracting reads via BAZAM";
    if [[ ${aln_linear} == *"bam"* ]]; then
        time bazam -bam ${input_aln} -L ${immunovar_bed} | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz
        time bazam -bam ${input_aln} -L chr2:179424038-179441143 | gzip > ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz
    else
        echo "Assuming we're working with a cram file";
        bazam_path=$(which bazam);bazam_path=${bazam_path%bin/*}share/bazam*;bazam_path=${bazam_path}"/bazam.jar"
        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${immunovar_bed} | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz
        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L chr2:179424038-179441143 | gzip > ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz
    fi
    samtools view -@8 -C ${input_aln} -T ${ref} -f 4 | samtools fastq | gzip > ${input_aln%.${aln_linear}}.unmapped.fastq.gz
fi
echo "Alignment to whole genome immunovariation graph";
cat $(ls ${input_aln%.${aln_linear}}.*.fastq.gz | grep -v "mapped") > ${input_aln%.${aln_linear}}.mapped.fastq.gz
time vg giraffe -i -f ${input_aln%.${aln_linear}}.mapped.fastq.gz -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam
time vg giraffe -f ${input_aln%.${aln_linear}}.unmapped.fastq.gz -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam
cat ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam > ${input_aln%.${aln_linear}}.bazam.grch38.combined.gam
echo "${input_aln%.${aln_linear}} ready for VG Flow filtering-->inference"
#######################
#
export i=${input_aln%.${aln_linear}}.bazam.grch38.combined.gam outdir=${PWD} graph=${graph} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir}
. ${bigfoot_dir}/filter_immune_subgraph.sh


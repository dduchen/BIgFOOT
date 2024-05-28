#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=BamBIgFOOT
#SBATCH --time=24:00:00
#SBATCH --mem=70GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out
#

export graphdir=${bigfoot_source}
export graph_base=whole_genome_ig_hla_kir_immunovar

bf_env_load=$(conda info --envs | grep "bigfoot" | awk '{print $NF}')
curenv=$(declare -p -x)
source ${HOME}/.bashrc;
eval "$curenv"
conda activate ${bf_env_load};

PATH=${tools_dir}:$PATH

input_aln=${i##*/};echo $input_aln
aln_linear=$(echo ${input_aln} | sed s/.*\\.//g)

if [ -s ${bigfoot_dir}/../custom_beds/custom_bed.bed ]; then
    echo "Custom bed exists";
else
    echo -e "chr2\t179424038\t179441143" | cat ${immunovar_bed} ${bigfoot_dir}/../custom_beds/grch38_FCGR_loci.bed - | bedtools sort -i - | bedtools merge -i - -d 100 > ${bigfoot_dir}/../custom_beds/custom_bed.bed
fi

if [ -s ${input_aln%.${aln_linear}}*combined.gam ]; then
    echo "Pre-processing completed for ${input_aln}"
else
    if [ -s ${input_aln%.${aln_linear}}.unmapped.fastq.gz ]; then
        echo "${i##*/} has been processed already - using it"
    else
        echo "Processing: ${i##*/}"
        if test -f "${input_aln}"; then
            echo "downloaded cram already - using it"
        else
            echo "Downloading ${i}"
            wget -c -nv ${i};
        fi
        input_aln=${i##*/}
        aln_linear=$(echo ${input_aln} | sed s/.*\\.//g)
        samtools index ${input_aln}
        bazam_path=$(which bazam);bazam_path=${bazam_path%bin/*}share/bazam*;bazam_path=${bazam_path}"/bazam.jar"
        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${bigfoot_dir}/../custom_beds/custom_bed.bed | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz
#        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${immunovar_bed} | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz
#        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L chr2:179424038-179441143 | gzip > ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz
#        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${bigfoot_dir}/../custom_beds/grch38_FCGR_loci.bed | gzip > ${input_aln%.${aln_linear}}.bazam.FCGR.fastq.gz
        samtools view -@8 -C ${input_aln} -T ${ref} -f 4 | samtools fastq | gzip > ${input_aln%.${aln_linear}}.unmapped.fastq.gz
        rm ${input_aln}
    fi
#    cat ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz ${input_aln%.${aln_linear}}.bazam.FCGR.fastq.gz ${input_aln%.${aln_linear}}.bazam.fastq.gz > ${input_aln%.${aln_linear}}.mapped.fastq.gz
    cat ${input_aln%.${aln_linear}}.bazam.fastq.gz > ${input_aln%.${aln_linear}}.mapped.fastq.gz
    time vg giraffe -i -f ${input_aln%.${aln_linear}}.mapped.fastq.gz -x ${graphdir}/${graph_base}.xg -H ${graphdir}/${graph_base}.gbwt -d ${graphdir}/${graph_base}.dist -m ${graphdir}/${graph_base}.min -p > ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam
    time vg giraffe -f ${input_aln%.${aln_linear}}.unmapped.fastq.gz -x ${graphdir}/${graph_base}.xg -H ${graphdir}/${graph_base}.gbwt -d ${graphdir}/${graph_base}.dist -m ${graphdir}/${graph_base}.min -p > ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam
    cat ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam > ${input_aln%.${aln_linear}}.bazam.grch38.combined.gam
    echo "${input_aln%.${aln_linear}} ready for bigfoot inference"
    if [ -s ${input_aln%.${aln_linear}}.bazam.grch38.combined.gam ]; then
#        rm ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz 
#        rm ${input_aln%.${aln_linear}}.bazam.FCGR.fastq.gz 
        rm ${input_aln%.${aln_linear}}.bazam.fastq.gz
        rm ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam
        rm ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam
    fi
fi
if [ -s ${input_aln%%\.*}*genotyping/*results_cleaned.txt ]; then
    echo "BIgFOOT Pipeline completed for ${input_aln}"
else
    echo "Progressing into BIgFOOT inference pipeline"
    export workdir=${PWD}; 
    export PATH=${tools_dir}:$PATH ;
    export bigfoot_dir=${bigfoot_dir};
    export bigfoot_source=${bigfoot_source};
    export graphdir=${bigfoot_source}; 
    export graph="wg_immunovar";
    export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar;
    export immune_graph=${graph_base}".subgraph"; 
    export valid_alleles=true;
    #
    export i=${input_aln%.${aln_linear}}.bazam.grch38.combined.gam
    . ${bigfoot_dir}/filter_immune_subgraph.sh
    rm ${input_aln%.${aln_linear}}.bazam.grch38.combined.gam
fi
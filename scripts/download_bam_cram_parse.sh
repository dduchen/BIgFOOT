#!/bin/bash

#
# for slurm-based batch job
#SBATCH --ntasks=1
#SBATCH --job-name=prep_bigfoot
#SBATCH --time=8:00:00
#SBATCH --mem=65GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out
#

export graphdir=${bigfoot_source}
export graph_base=whole_genome_ig_hla_kir_immunovar

input_aln=${i##*/};echo $input_aln
aln_linear=$(echo ${input_aln} | sed s/.*\\.//g)

if [[ ${input_aln} == *"gam" ]]; then
    echo "input is GAM file - reprocessing ${input_aln}"
    export reprocess=true
    if [ -s ${input_aln%%\.*}*filtered.gam ]; then
        echo "Reperforming analysis for ${input_aln}"
        export workdir=${PWD}; 
        export PATH=${tools_dir}:$PATH;
        export bigfoot_dir=${bigfoot_dir};
        export bigfoot_source=${bigfoot_source};
        export graphdir=${bigfoot_source}; 
        export graph="wg_immunovar";
        export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar;
        export immune_graph=${graph_base}".subgraph"; 
        export valid_alleles=true;
        export i=${input_aln%.${aln_linear}}.filtered.gam;
        export aln_type='pe';
        sample=${i%%.*}; sample=${sample##*/};
        sample_immune_graph=${sample}.immune_subset.${graph}.${aln_type};
        output_graph=${sample}_graph.${graph}.${aln_type}
        export i=${sample_immune_graph}.gam workdir=${PWD} sample_id=${sample} aln_type='pe' tools_dir=${tools_dir} use_augmented=FALSE bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} merged=FALSE graph=${graph} valid_alleles=${valid_alleles};
        if [ -s ${sample_immune_graph}.gam ]; then
            . ${bigfoot_dir}/gene_hap_IMGT_vgflow.sh;
            cd ${workdir}/${sample_id}_${graph}_genotyping/
            mkdir -p ${sample_id}.${graph}_vcf
            ls ./familywise_pe_haplotype_inference/*depth_cleaned.gfa > ${sample_id}.${graph}_allelic_inference.txt
            rm ${sample_id}.${graph}.alleles.fasta;
            for allele_fasta in $(cat ${sample_id}.${graph}_allelic_inference.txt); do echo "Extracting alleles: ${allele_fasta}"
                vg paths -Fv ${allele_fasta} | seqkit grep -r -p "IG|TR" >> ${sample_id}.${graph}.alleles.fasta
            done
            mv ${outdir}/*vcf ${sample_id}.${graph}_vcf
            rm $outdir/*plines; rm $outdir/*paths;
            rm -rf ${outdir}/tmp;
            Rscript ${bigfoot_dir}/clean_genewise_results.R
            tar -czvf ${sample_id}.${graph}_allelic_inference.tar.gz -T ${sample_id}.${graph}_allelic_inference.txt --remove-files
            ls ${outdir}/* > ${sample_id}.${graph}_other_materials.txt
            tar -czvf ${sample_id}.${graph}_other_materials.tar.gz -T ${sample_id}.${graph}_other_materials.txt --remove-files
            cd ${workdir}
            echo "fin!"
            ###########################
            # -- association explo -- #
            odgi unitig -i ${output_graph}.sub.gfa -t 31 -t 4 -P > ${output_graph}.sub.unitigs.fa
        elif [[ ${input_aln} == *"bazam.grch38.combined.gam" ]]; then
            echo "${input_aln} as input - preprocessing + bigfoot analysis";
            export i=${input_aln} workdir=${PWD} reprocess=true graph=${graph} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} tools_dir=${tools_dir} valid_alleles=${valid_alleles}
            . ${bigfoot_dir}/filter_immune_subgraph.sh
        else
            echo "${input_aln} inappropriate input format - rerun";
        fi
    else
        if [ -s ${input_aln%%\.*}*.bazam.grch38.combined.sorted.gam ]; then
            echo "Looks like we're re-running BIgFOOT - using ${input_aln%.${aln_linear}}.bazam.grch38.combined.sorted.gam"
            export i=$(ls ${input_aln%%\.*}*.bazam.grch38.combined.sorted.gam);
            . ${bigfoot_dir}/filter_immune_subgraph.sh
        elif [ -s ${input_aln%%\.*}*.bazam.grch38.combined.gam ]; then
            export i=$(ls ${input_aln%%\.*}*.bazam.grch38.combined.gam)
            . ${bigfoot_dir}/filter_immune_subgraph.sh
        else 
            echo "Restart from scratch for this sample"
        fi
    fi
else
    echo "Processing (B/CR)AM: ${input_aln}"
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
        immune_graph=${graph_base}".subgraph"
    elif [[ $graph == "ig_hla_kir" ]]; then
        echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IMGT-IPD alleles - Deprecated"
    else
        echo "Forgot to define a graph"
        exit
    fi
    #
    if [ -s ${bigfoot_dir}/../custom_beds/custom_bed.bed ]; then
        echo "Custom bed exists";
    else
        echo -e "chr2\t179424038\t179441143" | cat ${immunovar_bed} ${bigfoot_dir}/../custom_beds/grch38_FCGR_loci.bed - | bedtools sort -i - | bedtools merge -i - -d 100 > ${bigfoot_dir}/../custom_beds/custom_bed.bed
        echo -e "chrM\t0\t16569" | cat ${bigfoot_dir}/../custom_beds/custom_bed.bed - | bedtools sort -i - | bedtools merge -i - -d 100 > ${bigfoot_dir}/../custom_beds/custom_bed.tmp.bed && mv ${bigfoot_dir}/../custom_beds/custom_bed.tmp.bed ${bigfoot_dir}/../custom_beds/custom_bed.bed
    fi
    ref_build=grch38
    if [[ ${input_aln} == *.cram ]];then
        if [[ $ref_build == "grch38" ]]; then
            ref=${ref};
            if [ -s $ref ]; then
                echo "using GRCh38 reference";
            else 
                echo "Parsing CRAM require an indexed reference"
                wget -P ${bigfoot_source} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa;
                ref=${bigfoot_source}/GRCh38_full_analysis_set_plus_decoy_hla.fa
                samtools faidx ${ref};
            fi
        elif [[ $ref_build == "chm13" ]]; then
            ref="${ref}";
            immunovar_bed="${bigfoot_dir}/../custom_beds/chm13_custom_immunovar_coords.bed";
        else
            echo "Currently BIgFOOT only supports GRCh38-based BAM input. Future release will have CHM13 as an option"
        fi
    fi
    if [ -s ${input_aln%%\.*}*combined.sorted.gam ]; then
        echo "Pre-processing + sorting completed for ${input_aln}"
    elif [ -s ${input_aln%%\.*}*combined.gam ]; then
        echo "Pre-processing completed for ${input_aln}"
    else
        if [ -s ${input_aln%.${aln_linear}}.unmapped.fastq.gz ]; then
            echo "${i##*/} has been processed already - using it"
        else
            echo "Processing: ${i##*/}"
            if test -f "${input_aln}"; then
                echo "downloaded cram already - using it"
            else
                wget -c -nv ${i};
            fi
            input_aln=${i##*/}
            aln_linear=$(echo ${input_aln} | sed s/.*\\.//g)
            if [ -s ${input_aln}.bai ]; then
                echo "${input_aln} already indexed"
            else
                echo "indexing ${input_aln}"
                samtools index ${input_aln}
            fi
            bazam_path=$(which bazam);bazam_path=${bazam_path%bin/*}share/bazam*;bazam_path=${bazam_path}"/bazam.jar";
            time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${bigfoot_dir}/../custom_beds/custom_bed.bed | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz;
    #        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${immunovar_bed} | gzip > ${input_aln%.${aln_linear}}.bazam.fastq.gz
    #        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L chr2:179424038-179441143 | gzip > ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz
    #        time java -Xmx36g -Dsamjdk.reference_fasta=${ref} -jar ${bazam_path} -bam ${input_aln} -L ${bigfoot_dir}/../custom_beds/grch38_FCGR_loci.bed | gzip > ${input_aln%.${aln_linear}}.bazam.FCGR.fastq.gz
            if [[ ${input_aln} == *.cram ]];then
                echo "Cram input - ref required"
                samtools view -@8 -b ${input_aln} -T ${ref} -f 4 | samtools fastq | gzip > ${input_aln%.${aln_linear}}.unmapped.fastq.gz;
            else 
                samtools view -@8 -b ${input_aln} -f 4 | samtools fastq | gzip > ${input_aln%.${aln_linear}}.unmapped.fastq.gz;
            fi
            rm ${input_aln};
        fi
    #    cat ${input_aln%.${aln_linear}}.bazam.TTN.fastq.gz ${input_aln%.${aln_linear}}.bazam.FCGR.fastq.gz ${input_aln%.${aln_linear}}.bazam.fastq.gz > ${input_aln%.${aln_linear}}.mapped.fastq.gz
        cat ${input_aln%.${aln_linear}}.bazam.fastq.gz > ${input_aln%.${aln_linear}}.mapped.fastq.gz;
        time vg giraffe -i -f ${input_aln%.${aln_linear}}.mapped.fastq.gz -x ${graphdir}/${graph_base}.xg -H ${graphdir}/${graph_base}.gbwt -d ${graphdir}/${graph_base}.dist -m ${graphdir}/${graph_base}.min -p > ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam;
        time vg giraffe -f ${input_aln%.${aln_linear}}.unmapped.fastq.gz -x ${graphdir}/${graph_base}.xg -H ${graphdir}/${graph_base}.gbwt -d ${graphdir}/${graph_base}.dist -m ${graphdir}/${graph_base}.min -p > ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam;
        cat ${input_aln%.${aln_linear}}.bazam.grch38.wg.gam ${input_aln%.${aln_linear}}.unmapped.grch38.wg.gam > ${input_aln%.${aln_linear}}.bazam.grch38.combined.gam;
        echo "${input_aln%.${aln_linear}} ready for bigfoot inference";
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
        if [ -s ${input_aln%%\.*}*.bazam.grch38.combined.sorted.gam ]; then
            echo "Looks like we're re-running BIgFOOT - using ${input_aln%.${aln_linear}}.bazam.grch38.combined.sorted.gam"
            export i=${input_aln%.${aln_linear}}.bazam.grch38.combined.sorted.gam
        else 
            export i=${input_aln%.${aln_linear}}.bazam.grch38.combined.gam
        fi
        . ${bigfoot_dir}/filter_immune_subgraph.sh
    fi
fi
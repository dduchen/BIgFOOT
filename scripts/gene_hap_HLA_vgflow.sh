#!/bin/bash


# for slurm-based batch job
#SBATCH --ntasks=1
#SBATCH --job-name=run_hla_flow
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
workdir=${PWD}
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
    graphdir=${bigfoot_source}
    graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar
    genotyping_nodes_dir=${bigfoot_source}/wg_ig_hla_kir_immunovar_genotyping_nodes/
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
    vg gamsort ${gam_file} -i ${gam_file%.gam}.sorted.gam.gai -p > ${gam_file%.gam}.sorted.gam
fi
#
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
for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "^TR\|^IGH\|^IGLV\|^IGKV" | grep -v "__\|IGHD");do echo $each;
# KIR inference only:
#for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep "KIR");do echo $each;
# HLA inference only:
for each in $(ls ${genotyping_nodes_dir} | grep "nodes.txt" | grep -v "TR\|IGH\|IGK\|IGL\|KIR");do echo $each;
# --
cd ${datadir}
sample=${each%.nodes.txt}
gene=${sample%.immune_subset}
gene_actual=$(echo $gene | sed 's!__!/!g')
#
loci=$(vg paths -Lv ${graph_base}.xg | grep "#1#${gene}\*" | cut -f1 -d"#" | sort | uniq | grep "IMGT\|HLA\|KIR")
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
vg find -x ${graph_base}.xg -c 150 -L -N ${genotyping_nodes_dir}/${each} > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
echo "Extracting local haplotypes reads --> ${sample_id}.${graph}.${gene}.haplotypes.gfa"
vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "grch\|chm" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt
vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.txt | 
vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg
vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.vg -F > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta
vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg | grep "IGL\|hirh_H\|hirh_N\|MHC\|${gene_actual}" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt
vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt | 
vg mod - -N > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg.tmp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg
vg find -x ${graph_base}.xg -l ${gam_file%.gam}.sorted.gam -A ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam
#
if [[ "$loci" =~ ^(HLA)$ ]]; then 
    downsampled=FALSE
    if [[ "$downsampled" =~ ^(TRUE)$ ]]; then
        echo "basing inference on haplotypes + alleles embedded in graph, for inference based on full set of alleles set downsampled=FALSE (this might be much slower for certain HLA genes)"
    # haplotype graph:
        vg paths -v ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -F -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.txt > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
        seqkit rmdup -s < ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta
#        seqkit grep -r -p ${gene} ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${gene}.alleles.fasta
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
        if [ -s  ${outdir}/${gene}.haps.fasta ]; then
            echo "Local haplotypes found for: ${gene_actual}"
            seqkit seq --min-len ${gene_min_len} --max-len 1000000  ${outdir}/${gene}.haps.fasta > ${outdir}/${gene}.haps.fasta.tmp && mv ${outdir}/${gene}.haps.fasta.tmp ${outdir}/${gene}.haps.fasta
            # retain haplotypes with perfect match to one of our alleles
            mkdir -p $outdir/${gene}_haps
            # allow for an edit distance less than 10bp
#            time minimap2 -x asm20 -c ${outdir}/${gene}.alleles.fasta ${outdir}/${gene}.haps.fasta | cut -f1,13 | egrep "NM:i:0|NM:i:[0-9]$" | cut -f1 | sort | uniq > ${outdir}/${gene}_haps/haps.matching.txt
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
#        time minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf
        blend -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf
        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.paf -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
        ${tools_dir}/gfaffix ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.tmp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
    #odgi explode -i ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -b 1 -s P -O -P -p "${sample_id}.${gene}.exp" -g
    #mv ${sample_id}.${gene}.exp*gfa ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
        vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
        vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg
        vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.vg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa
        echo "Augmenting haplotype graph with sample-specific variation for association testing"
        vg convert -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gfa -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg
        vg index -p -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg;
        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -P --pass-paths
        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz --gbz-format -P --pass-paths;
        vg index -p -j ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz
        vg minimizer -p -d ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist -o ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.min ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz
#        vg stats -b ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.dist | sort -nr -k 3 | head -n 5
        rm ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_paths.gbwt
        vg prune -p -u -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -k 45 -m ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
        vg index -p -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pruned.vg
# novel personalized graph-based mapping:
#        vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg | grep "HLA" > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_personalized.txt
#        vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_personalized.txt -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.personalized.xg
#        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.personalized.xg -o ${outdir}/${sample_id}.${graph}.${gene}.personalized_paths.gbwt -E;
#        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.personalized.xg -o ${outdir}/${sample_id}.${graph}.${gene}.personalized.gbwt -a ${outdir}/${sample_id}.${graph}.${gene}.personalized_paths.gbwt -n 64 --pass-paths;
#        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.personalized.xg -g ${outdir}/${sample_id}.${graph}.${gene}.personalized.gbz --gbz-format -a ${outdir}/${sample_id}.${graph}.${gene}.personalized_paths.gbwt -n 64 --pass-paths;
#        vg gbwt -p --num-threads 16 -r ${outdir}/${sample_id}.${graph}.${gene}.personalized.ri -Z ${outdir}/${sample_id}.${graph}.${gene}.personalized.gbz
#        vg index -j ${outdir}/${sample_id}.${graph}.${gene}.personalized.dist ${outdir}/${sample_id}.${graph}.${gene}.personalized.gbz
#        export TMPDIR2=${outdir}/tmp
#        mkdir -p ${TMPDIR2}
#        vg view -X ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq
#        kmc -k29 -m128 -okff -t16 -hp ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq ${TMPDIR2}/${gene}.personalized ${TMPDIR2}
#        vg haplotypes -v 2 -t 16 -H ${outdir}/${sample_id}.${graph}.${gene}.personalized.hapl -k ${TMPDIR2}/${gene}.personalized.kff -g ${TMPDIR2}/${gene}.personalized.gbz ${outdir}/${sample_id}.${graph}.${gene}.personalized.gbz --subchain-length 1000 --diploid-sampling
# subsampling approach
#        bbnorm.sh in=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq out=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz target=30 mindepth=1 threads=8
#        vg giraffe -p -t 16 -Z ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz --haplotype-name ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.hapl --kff-name ${TMPDIR2}/sample.kff -N sample -f ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#        time vg map -f ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
# use VG map if lots of reads...
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
#        vg view -X ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq
#        bbnorm.sh in=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq out=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz target=50 mindepth=2 threads=8
#        ~/tools/rasusa/rasusa -i ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq --coverage 100 -s 123 --genome-size ${gene_min_len} -o ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz
#        seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz
#        time vg map -f ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
# if bed input - # make ifelse statement using "bwa" match or something
#         # bbnorm rather than extreme filtering based on giraffe alignment
#        vg filter -r 0.66 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.vg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.filt.gam
# giraffe-mapping to downsampled haplotype+allele graph --> filtering again?
#        time vg giraffe -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.filt.gam -Z ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbz -H ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#       vg filter -r 0.9 -P -q 60 -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${gene}.tmp.gam && mv ${gene}.tmp.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
# k-mer-based coverage downsampling
#        vg view -X ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.filt.gam > ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq
#        bbnorm.sh in=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq out=${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz target=100 mindepth=2 threads=8
#        seqkit stats ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz
#        time vg map -f ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.fastq.gz -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#        time vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.filt.gam -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -g ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#        time vg giraffe -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -Z ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbz -H ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
    # copy mapped reads -- use this for augmentation?
        cp ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam
#        allele_adaptive_filt=$(grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.fasta | grep "HLA" | wc -l)
#        allele_adaptive_filt=$(bc -l <<< "scale=2;((100-${allele_adaptive_filt}/2)/100)")
        #echo "Filtering reads using pairwise similarity of: ${allele_adaptive_filt} to aligned region"
#       vg filter -r ${allele_adaptive_filt} -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#        vg filter -r 0.98 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
        vg filter -r 0.8 -P -s 1 -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.prefilt.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
#       filter here ^^^^^^ makes big difference
        vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
# 1) should perform inference on haplotype graph with non-allele paths removed?
        echo "Limiting inference to alleles of interest"
        vg paths -Lv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg | grep "HLA" > ${outdir}/${sample_id}.${graph}.${gene}.alleles
        vg paths -r -p ${outdir}/${sample_id}.${graph}.${gene}.alleles -x ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.pg > ${outdir}/${sample_id}.${graph}.${gene}.vg
        vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
        vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.vg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
        depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
        min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.05"| awk '{printf("%d\n",$1 + 0.5)}')
        cd ${outdir}
        python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
#        time python3 ${bigfoot_dir}/vg-flow_immunovar.py --careful --optimization_approach 2 --min_depth 0 --trim 0 -m 0 -c 2 --remove_included_paths 0 --ilp --max_strains 2 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
        time python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --min_depth 0 --trim 0 --greedy_mode all --ilp --max_strains 2 -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 --wfmash_param ${wfmash_param} ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
#       python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --careful --min_depth 0 --trim 0 -m 0 -c ${min_strain_depth} --remove_included_paths 0 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
        mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta ; mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.haps.final.fasta ; mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph.gfa
        rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
        Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
        echo "Allele-level abundance estimation completed for ${gene} ::"
        grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta
    else 
        echo "basing inference on full set of alleles for ${gene}, for inference based on embedded haplotypes+alleles set downsampled=TRUE"
# todo 1) - in cleanup - need to remove wfmash... directories like seqwish
#      2) - extract reference backbone to embed variation along with called alleles
        echo "MHC/HLA locus :: Using pre-constructed IPD-based HLA graphs for inference"
        # re-align to haplotype graph + embed variation with adequate support
        time vg map -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -g ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gcsa -1 ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
    #    time vg giraffe -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -Z ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbz -H ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.gbwt -p > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
        vg filter -r 0 -P -s 1 -x ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.tmp.gam && mv ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.tmp.gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam
        vg depth --gam ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth;
        depth_locus=$(awk -F ' ' '{print $1}' ${outdir}/${sample_id}.${graph}.${gene}.filtered.depth)
        min_strain_depth=$(bc -l <<< "scale=2;${depth_locus}*0.10"| awk '{printf("%d\n",$1 + 0.5)}')
        echo "Minimum strain depth required: ${min_strain_depth}"
        vg convert -fW ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.gfa
        vg view -a ${outdir}/${sample_id}.${graph}.${gene}.haplotypes.gam > ${outdir}/${sample_id}.${graph}.${gene}.vgflow.aln.json
        python3 ${bigfoot_dir}/parse_graph_vgflow.py --sample ${outdir}/${sample_id}.${graph}.${gene}.vgflow -m 0
#   skip some overlapping path-related QC + wfmash instead of minimap - segment length 33% shortest allele, % pairwise identity = 95%
# -- premake overlaps file etc?
#   python3 ${bigfoot_dir}/vg-flow_immunovar.py --min_depth 0 --trim 0 --greedy_mode all -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
#    time python3 ${bigfoot_dir}/vg-flow_immunovar.py --min_depth 0 --trim 0 --greedy_mode all -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
# VG-Flow - minimap2 
# real    13m33.262s
# user    13m27.818s
# sys     0m2.291s
        gene_min_len=$(vg paths -Ex ${genotyping_nodes_dir}${loci}_${gene}.haplotypes.xg | grep "HLA" | cut -f2 | sort | uniq | head -1)
        wfmash_param=$(bc -l <<< "scale=2;${gene_min_len}/10" | awk '{printf("%d\n",$1 + 0.5)}')
        if [ "${wfmash_param}" -lt 100 ]; then
            wfmash_param=100
        fi
# can solve using ILP knowing max 2 alleles ('--ilp --max_strains 2')
        time python3 ${bigfoot_dir}/vg-flow_immunovar_long_contigs.py --min_depth 0 --trim 0 --greedy_mode all --ilp --max_strains 2 -m 0 -c ${min_strain_depth} --remove_included_paths 0 --threads 16 --wfmash_param ${wfmash_param} ${outdir}/${sample_id}.${graph}.${gene}.vgflow.node_abundance.txt ${outdir}/${sample_id}.${graph}.${gene}.vgflow.final.gfa
# VG-Flow - wfmash (also correct)
# segment param = 33% of shortest allele:
# real    3m3.114s
# user    2m59.941s
# sys     0m1.591s
        mv trimmed_contigs.fasta ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
        mv haps.final.fasta ${outdir}/${sample_id}.${graph}.${gene}.haps.final.fasta
        mv genome_graph.gfa ${outdir}/${sample_id}.${graph}.${gene}.genome_graph.gfa
        rm genome_graph.gt; rm haps.fasta; rm overlaps.minimap2.paf; rm trimmed_contigs.paths ; rm trimmed_contigs.gfa
        Rscript ${bigfoot_dir}/parse_vgflow_output.R ${outdir}/${sample_id}.${graph}.${gene}.contigs.fasta
        echo "Allele-level abundance estimation completed for ${gene} ::"
        grep ">" ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta
# 2) augment annotated post-flow inference graph with reads for association testing
# embed along with reference for variant calling
        cat ${outdir}/${sample_id}.${graph}.${gene}.haps.final.annot.fasta ${outdir}/${sample_id}.${graph}.${gene}.haplotypes_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
        echo "Embedding novel variation with adequate support (~strain depth) to inferred flow graph + local CHM13 reference sequence"
        sed s/' '/_/g -i ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta
#        minimap2 -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
        blend -x asm20 -t 16 -c -X ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf
        seqwish -s ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.fasta -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.paf -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -b ${outdir}/seqwish_${sample_id}.${graph}
        gfaffix ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp; mv ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.tmp ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
        vg mod -n -U 10 -c ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa -X 256 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg
        vg convert -p ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg
        vg convert -fW ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.vg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gfa
        vg index -p -t 16 -L -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg;
#        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -P --pass-paths -n 64
        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -o ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -P --pass-paths
        vg gbwt -p -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbz --gbz-format -P --pass-paths;
        vg prune -p -u -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -k 31 -m ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pg > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
        vg index -p -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -f ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.node_mapping ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.pruned.vg
#   map locus-associated reads to augmented/annotated graph
        vg map -N ${sample_id}.${graph}.${gene} -G ${outdir}/${sample_id}.${graph}.${gene}.genotyping.immune_subset.gam -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -g ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gcsa -1 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gbwt -t 4 -M 1 > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam
        vg filter -r 0 -P -s 1 -q 60 -x ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.xg -D 0 -fu -t 4 ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.gam -v > ${outdir}/${sample_id}.${graph}.${gene}.genome_graph_ref.filt.gam
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
#
done
#
ls ${outdir}/*_files.txt > ${outdir}/${sample_id}_files_rm.txt
xargs rm < ${outdir}/${sample_id}_files_rm.txt
# clean up:
echo "All cleaned up!"
done
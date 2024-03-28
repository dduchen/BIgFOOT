#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=vgflow_prep
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=dylan.duchen@yale.edu
#SBATCH --mail-type=FAIL
#SBATCH --partition=pi_kleinstein
#SBATCH --output=slurm-%x.%j.out

DATAPATH=$workdir
OUTPATH=$workdir
cd $OUTPATH

sample=${i%%.*}

#
if [[ $graph == "minimal" ]]; then
    echo "Using minimal SV-graph"
    graph_base="/home/dd392/project/sv_graphs/sv_graph_minimal/wg_sv_immunovar_simple.wDups"
    graphdir="/home/dd392/project/sv_graphs/sv_graph_minimal/"
elif [[ $graph == "hirh" ]]; then
    echo "Using HIRH-embedded SV-graph -- need to update this, resorting to minimal graph"
    graph_base="/home/dd392/project/sv_graphs/sv_graph_minimal/wg_sv_immunovar_simple.wDups"
    graphdir="/home/dd392/project/sv_graphs/sv_graph_minimal/"
elif [[ $graph == "franken" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/grch38_chm13_immunovar.sort"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
elif [[ $graph == "wg_immunovar" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + IG Haplotypes + IPD HLA and KIR alleles"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/whole_genome_ig_hla_kir_immunovar"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
elif [[ $graph == "ig_hla_kir" ]]; then
    echo "Using IGenotyper GRCh38 + CHM13 + HIRH-embedded graph + IPD HLA and KIR alleles"
    graph_base="/home/dd392/project/grch38_chm13_immunovar/grch38_chm13_immunovar.hla_kir"
    graphdir="/home/dd392/project/grch38_chm13_immunovar/"
else
    echo "Define graph"
fi
#
# suffix - based on ngmerged output or qc-passed PE reads
if [[ $i == *"ngmerge"* ]]; then
    echo "Alignment of NGMerged fastq data"
    aln_type="ngmerged";
else echo "Alignment of non-merged reads";
    aln_type="pe"
fi
#
sample_immune_aln=${sample}.${graph}.${aln_type}
immune_graph=${graph_base}.subgraph
#
/home/dd392/tools/vg find -G ${i} -c 0 -x ${graph_base}.xg > ${sample_immune_aln}.vg
/home/dd392/tools/vg augment -iSs ${sample_immune_aln}.vg ${immune_graph}.paths.gam > ${sample_immune_aln}.tmp; mv ${sample_immune_aln}.tmp ${sample_immune_aln}.vg
echo "need to normalize/re-order nodes for processing - realign the GAM to filtered/normalized graph"
/home/dd392/tools/vg mod ${sample_immune_aln}.vg -n -c -X 32 - | /home/dd392/tools/vg ids -c - > ${sample}.${aln_type}.tmp; mv ${sample}.${aln_type}.tmp ${sample_immune_aln}.vg
/home/dd392/tools/vg index -L -x ${sample_immune_aln}.xg -t 4 ${sample_immune_aln}.vg -p
#vg gbwt -x ${sample}.${aln_type}.xg -o ${sample}.${aln_type}.gbwt -P --pass-paths
/home/dd392/tools/vg gbwt -p -x ${sample_immune_aln}.xg -o ${sample_immune_aln}.gbwt -P --pass-paths
/home/dd392/tools/vg gbwt -p -x ${sample_immune_aln}.xg -g ${sample_immune_aln}.gbz --gbz-format -P --pass-paths
/home/dd392/tools/vg index -p -j ${sample_immune_aln}.dist ${sample_immune_aln}.gbz
/home/dd392/tools/vg minimizer -p -d ${sample_immune_aln}.dist -o ${sample_immune_aln}.min ${sample_immune_aln}.gbz
# time consuming to construct gcsa- prune first? 
/home/dd392/tools/vg prune -p -u -g ${sample_immune_aln}.gbwt -k 45 -m ${sample_immune_aln}.node_mapping ${sample_immune_aln}.vg > ${sample_immune_aln}.pruned.vg
/home/dd392/tools/vg index -p -g ${sample_immune_aln}.gcsa -f ${sample_immune_aln}.node_mapping ${sample_immune_aln}.pruned.vg
/home/dd392/tools/vg snarls ${sample_immune_aln}.xg -T >${sample_immune_aln}.snarls;
/home/dd392/tools/vg map -G ${i} -x ${sample_immune_aln}.xg -g ${sample_immune_aln}.gcsa -1 ${sample_immune_aln}.gbwt > ${sample_immune_aln}.gam
/home/dd392/tools/vg filter -r 0.1 -P -s 1 -fu -t 4 ${sample_immune_aln}.gam -v > ${sample_immune_aln}.filtered.gam
#/home/dd392/tools/vg giraffe -Z ${sample_immune_aln}.gbz -m ${sample_immune_aln}.min -d ${sample_immune_aln}.dist -G $i -p -N '$sample' > ${sample_immune_aln}.gam
#/home/dd392/tools/vg filter -r 0.75 -s 2 -fu -t 4  ${sample_immune_aln}.gam -v > ${sample_immune_aln}.filtered.gam
/home/dd392/tools/vg convert -fW ${sample_immune_aln}.vg > ${sample}.gfa
/home/dd392/tools/vg view -a ${sample_immune_aln}.filtered.gam > ${sample}.aln.json
# -- parse using vg-flow python code
#parse_graph_vgflow.py
# /home/dd392/project/parse_graph_vgflow.py
# conda activate vg-flow-env
/gpfs/ycga/project/kleinstein/dd392/conda_envs/vg-flow-env/bin/python /home/dd392/project/sv_graphs/parse_graph_vgflow.py --sample ${sample}
#
echo "VG-Flow files prepped - ${sample}.${graph}.${aln_type}.final..."
mv ${sample}.final.gfa ${sample}.${graph}.${aln_type}.final.gfa
grep -v "_alt_*\|grch38*" ${sample}.${graph}.${aln_type}.final.gfa > ${sample}.${graph}.${aln_type}.tmp && mv ${sample}.${graph}.${aln_type}.tmp 
${sample}.${graph}.${aln_type}.final_no_sv.gfa
mv ${sample}.final.gt ${sample}.${graph}.${aln_type}.final.gt
mv ${sample}.edge_abundances.txt ${sample}.${graph}.${aln_type}.final.edge_abundances.txt
mv ${sample}.node_abundance.txt ${sample}.${graph}.${aln_type}.final.node_abundance.txt
#
echo "ready for VG-Flow - vg-flow_immunovar.py - w/ ${sample}.${graph}.${aln_type}.final.node_abundance.txt  ${sample}.${graph}.${aln_type}.final_no_sv.gfa"
#
# cleaning assumes ngmerged route takes sufficiently loner - pe workflow has completed
if [[ $aln_type == *"pe"* ]]; then
# perform some depth estimates:
    /home/dd392/tools/vg depth --gam ${i} ${immune_graph}.xg > ${sample}.${graph}_immune.pe.depth;
    depth_immune=$(awk -F ' ' '{print $1}' ${sample}.${graph}_immune.pe.depth)
    echo "${sample} depth across immunovariation loci: ${depth_immune}"
    # TTN exon 327 depth - immunotyper-SR reference
    #chr2:179424038-179441143
    grch_path_chr2=$(/home/dd392/tools/vg paths -Lv ${graph_base}.xg | grep "grch38#chr2:")
    /home/dd392/tools/vg depth --gam ${i} ${graph_base}.xg -p "${grch_path_chr2}:179424038[-179441143]" > ${sample}.${graph}_TTN.pe.depth;
    depth_ttn=$(awk -F ' ' '{print $1}' ${sample}.${graph}_TTN.pe.depth)
    echo "${sample} depth across Titin (TTN) exon 327 (Immunotyper-sr ref): ${depth_ttn}"
# finished!
    echo "fin";
#    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='pe',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/familywise_vgflow.sh;
#    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='pe',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/gene_hap_simple_vgflow.sh;
    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='pe',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/gene_hap_IMGT_vgflow.sh
else echo "fin - performing depth estimates and cleaning up after myself";
# perform some depth estimates:
    /home/dd392/tools/vg depth --gam ${i} ${immune_graph}.xg > ${sample}.${graph}_immune.ngmerged.depth;
    depth_immune=$(awk -F ' ' '{print $1}' ${sample}.${graph}_immune.ngmerged.depth)
    echo "${sample} depth (NGmerged) across immunovariation loci: ${depth_immune}"
    # TTN exon 327 depth - immunotyper-SR reference
    #chr2:179424038-179441143
    grch_path_chr2=$(/home/dd392/tools/vg paths -Lv ${graph_base}.xg | grep "grch38#chr2:")
    /home/dd392/tools/vg depth --gam ${i} ${graph_base}.xg -p "${grch_path_chr2}:179424038[-179441143]" > ${sample}.${graph}_TTN.ngmerged.depth;
    depth_ttn=$(awk -F ' ' '{print $1}' ${sample}.${graph}_TTN.ngmerged.depth)
    echo "${sample} depth (NGmerged) across Titin (TTN) exon 327 (Immunotyper-sr ref): ${depth_ttn}"
# finished!
    ls ${sample}.${graph}* | grep -v "_1.fastq\|_2.fastq" | grep -v "filtered\|final\|immune_subset\|sorted.gam\|fix.gfa\|depth" > ${sample}.${graph}_files.txt
    xargs rm < ${sample}.${graph}_files.txt
# run final gene-wise vg-flow inference
#    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='ngmerged',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/familywise_vgflow.sh;
# need to update merged script - but deprecated for now
    sbatch --export=i=${i%.gam},sample_id=$sample,graph=$graph,aln_type='ngmerged',use_augmented=FALSE,workdir=${PWD} /home/dd392/project/gene_hap_simple_vgflow.sh;
fi

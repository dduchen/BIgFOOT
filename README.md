# BIgFOOT: inference of Biological I(g)mmunovariation via Graph FOOTprinting
### Current version: 0.0.1

This workflow infers alleles, calls novel variation, and constructs sample-specific sequence variation graphs for immunoglobulin(Ig)/other immune-related loci which can be used to perform genetic association tests. The workflow inolves a combination of various graph-construction steps, sequence-to-graph alignment, flow graph decomposition, and unitig calling.<br>
I hope to expand this workflow to enable genome-to-genome analyses/assessing genetic associations between host germline immunovariation and pathogen/metagenomic genetic variation/diversity (i.e., searching for immunological <b>FOOTprints</b> using joint host-pathogen genomic data).

<i>Genetic loci where BIgFOOT performs accurate allele calling:</i>
- IGH
- IGL
- HLA (DQA1/DQB1/... more to come)<br>

<i>Infers alleles - but, like bigoot, I have no evidence they're real (WiP):</i><br>
- IGK
- TR
- KIR

### Input: 
- Raw fastq(.gz)
- BAM/CRAM alignment
<i>Note: you'll need ~65GB of RAM to sucessfully perform sequence-to-graph alignment against the full genome immunovariation graph</i>

## Set up conda environment
<b>BIgFOOT is heavily influenced/relies on methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow (v0.0.4)</a>. </b>

1) Clone me! <code>git clone https://github.com/dduchen/BIgFOOT.git</code>
2) set up conda/mamba environment we'll be needing -- can move some of these after the '#' if they're already in your path (e.g., samtools, we assume you have R)<br>
<code>mamba create --name bigfoot -c bioconda -c conda-forge -c gurobi python=3 fastp graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit bbmap minimap2 seqwish blend-bio wfmash samtools pyseer unitig-caller parallel #fastq-dl kmc r-base cd-hit
conda activate bigfoot </code><br>
Ensure you have an active gurobi licence:<br>
<code>gurobi_cl</code><br>
<i>We also use the following R/bioconductor packages: </i><br>
- data.table;
- dplyr;
- Biostrings/DECIPHER </code>

3) We also use some external tools which need to be accessible in your PATH<br>
<code>tools_dir=~/tools;</code> # (wherever you normally install+store software)<br>
<code>PATH=$PATH:${tools_dir};</code><br>
<code>cd ${tools_dir};</code><br>

### Download BIgFOOT graph materials from zenodo <a href="https://doi.org/10.5281/zenodo.10869771"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10869771.svg" alt="DOI"></a><br>

<code>bigfoot_source=~/pi_kleinstein/bigfoot</code> # where are we storing all of the reference graph files? <br>
<code>bigfoot_source=${tools_dir}/bigfoot</code><br>
<code>mkdir -p ${bigfoot_source} </code><br>
<code>wget -P ${bigfoot_source} "https://zenodo.org/records/10869771/files/immunovar_graph_materials.tar.gz?download=1" </code><br>
<code>cd ${bigfoot_source} ; tar -xvf ${bigfoot_source}/immunovar_graph_materials.tar.gz* --keep-newer-files</code><br>
Make distance indexes read only<br>
<code>chmod 0444 *.dist</code><br>
We also need the variation graph toolkit (VG) executable<br>
- <code>wget -P ${tools_dir}/ https://github.com/vgteam/vg/releases/download/v1.56.0/vg; chmod +x ${tools_dir}/vg <br>
PATH=${tools_dir}:$PATH</code><br>

We use Ryan Wick's Assembly-dereplicator package during haplotype selection <a href="https://github.com/rrwick/Assembly-Dereplicator">Assembly-dereplicator</a>.<br>
- <code>git clone https://github.com/dduchen/Assembly-Dereplicator.git ${tools_dir}/Assembly-dereplicator </code><br>
We provide the option of using merged paired-end reads from NGmerge for alignment/inference (optional, not always recommended) <a href="https://github.com/harvardinformatics/NGmerge">NGmerge</a>.<br></code>
- <code>git clone https://github.com/dduchen/NGmerge.git ${tools_dir}/NGmerge </code><br>

### Running bigfoot - Example using sequencing/alignment files from ISGR: <a href="https://www.internationalgenome.org/data-portal/sample/NA19240">NA19240</a><br>
#### Yoruba in Ibadan, Nigeria, African Ancestry<br>
Set up example directory, download relevant files, and then run BIgFOOT pipeline<br>
- <code>conda activate bigfoot<br>
bigfoot_dir=${bigfoot_source}/scripts</code> (Change this if you've downloaded the github repo somewhere else/have the bigfoot analysis scripts saved elsewere)<br>
<code>bigfoot_dir=${tools_dir}/BIgFOOT/scripts ; immunovar_bed=${bigfoot_source}/grch38_custom_immunovar_coords.bed<br>
test_dir=${bigfoot_source}/example/ ; mkdir -p ${test_dir}; cd ${test_dir}<br></code>

#### Starting from raw reads (WES)<br>
<i>Illumina chemistry: V2, Array: Agilent Sure Select Whole exome capture 50 Mb</i><br>
- <code>#fastq-dl -a SRR507323 -o ${test_dir}/<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507323/SRR507323_1.fastq.gz<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507323/SRR507323_2.fastq.gz
- export sample="SRR507323" workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} merged="FALSE" graph="wg_immunovar" valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/preprocess_wg_immunovar_alignment.sh<br>
################################################################</code>

##### Starting from BAM/CRAM (WGS)<br>
- <code>wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram
- export bam_file="NA19240.final.cram" workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} ref_build="grch38" ref="${bigfoot_source}/GRCh38_full_analysis_set_plus_decoy_hla.fa" merged="FALSE" graph="wg_immunovar" valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/process_from_bam_wg_immunovar_alignment.sh
################################################################</code><br>
<i>Support for CHM13-based BAM/CRAM is planned</i>

##### Starting from subset of reads<br>
- <code>graphdir=${bigfoot_source};graph="wg_immunovar";graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar;immune_graph=${graph_base}".subgraph";<br>
bazam_reads=${i};
sample_id=${bazam_reads%.bazam.fastq.gz};sample_id=${sample_id##*\/};
vg giraffe -i -f ${bazam_reads} -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.bazam.grch38.wg.gam
vg giraffe -f ${sample_id}.unmapped.fastq.gz -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.unmapped.grch38.wg.gam
cat ${sample_id}.bazam.grch38.wg.gam ${sample_id}.unmapped.grch38.wg.gam > ${sample_id}.bazam.grch38.combined.gam
echo "${sample_id} ready for VG Flow filtering-->inference"
- export i=${sample_id}.bazam.grch38.combined.gam workdir=${PWD} graph=${graph} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/filter_immune_subgraph.sh
################################################################</code><br>

##### Parallel processing using GNU parallel
<code>conda activate bigfoot <br>
ls *bazam.fastq.gz > process_sample_ids.txt <br>
export workdir=${PWD}; export bigfoot_dir=${bigfoot_dir}; \ <br>
export graphdir=${bigfoot_source}; export graph="wg_immunovar"; \ <br>
export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar; \ <br>
export immune_graph=${graph_base}".subgraph"; export valid_alleles=true; <br>
for i in $(cat process_sample_ids.txt);do echo ${i}; <br>
    cd ${workdir}; <br>
    sample_id=${i%.bazam.fastq.gz};sample_id=${sample_id##*\/}; <br>
    cat ${sample_id}.bazam*fastq.gz > ${sample_id}.mapped.fastq.gz; <br>
    bazam_reads=${sample_id}.mapped.fastq.gz; <br>
    sample_id=${bazam_reads%.mapped.fastq.gz};sample_id=${sample_id##*\/}; <br>
    if [ -s ${sample_id}.bazam.grch38.wg.gam ]; then <br>
        echo "Alignment of linearly mapped reads completed"; <br>
    else <br>
        vg giraffe -i -f ${bazam_reads} -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.bazam.grch38.wg.gam <br>
    fi <br>
    if [ -s ${sample_id}.unmapped.grch38.wg.gam ]; then <br>
        echo "Alignment of unmapped reads completed"; <br>
    else <br>
        vg giraffe -f ${sample_id}.unmapped.fastq.gz -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.unmapped.grch38.wg.gam <br>
    fi <br>
    if [ -s ${sample_id}.bazam.grch38.combined.gam ]; then <br>
        echo "Graph alignment of unmapped reads completed"; <br>
    else <br>
        cat ${sample_id}.bazam.grch38.wg.gam ${sample_id}.unmapped.grch38.wg.gam > ${sample_id}.bazam.grch38.combined.gam <br>
    fi <br>
    echo "${sample_id} ready for VG Flow filtering-->inference" <br>
    i=${sample_id}.bazam.grch38.combined.gam; <br>
    . ${bigfoot_dir}/filter_immune_subgraph.sh <br>
done <br></code>

##### Process 1kGenomes individuals<br>

<code>ls *bazam.grch38.combined.gam > run_pipeline_sample_ids.txt
grep -f igh_samples.txt run_pipeline_sample_ids.txt > run_pipeline_sample_ids_igh.txt
grep -f igh_samples.txt -v run_pipeline_sample_ids.txt > run_pipeline_sample_ids_other.txt

parallel -j 8 'export workdir=${PWD}; export tools_dir=~/tools;
export PATH=${tools_dir}:$PATH ; \
export bigfoot_dir=${bigfoot_dir}; \
export bigfoot_source=${bigfoot_source}; \
export graphdir=${bigfoot_source}; export graph="wg_immunovar"; \
export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar; \
export immune_graph=${graph_base}".subgraph"; export valid_alleles=true;
export i={}; \
time . ${bigfoot_dir}/filter_immune_subgraph.sh > ${i%.final.bazam.*}_bigfootprint.txt' :::: <(cat run_pipeline_sample_ids_other.txt);
time . ${bigfoot_dir}/filter_immune_subgraph.sh > ${i%.final.bazam.*}_bigfootprint.txt' :::: <(cat run_pipeline_sample_ids_igh.txt);
time . ${bigfoot_dir}/filter_immune_subgraph.sh > ${i%.final.bazam.*}_bigfootprint.txt' :::: <(cat run_pipeline_sample_ids_other.txt);

#bigfoot_dir=~/tools/BIgFOOT/scripts/
bigfoot_dir=~/Documents/github/BIgFOOT/scripts/
#bigfoot_source=~/pi_kleinstein/bigfoot/
bigfoot_source=/home/dduchen/Documents/bigfoot/

split -l 20 run_pipeline_sample_ids_other.txt process_sample_split_
for i in $(ls process_sample_split_* );do echo $i;
    time parallel -j 4 'export workdir=${PWD}; export tools_dir=~/tools;
    export PATH=${tools_dir}:$PATH ; \
    export bigfoot_dir=${bigfoot_dir}; \
    export bigfoot_source=${bigfoot_source}; \
    export graphdir=${bigfoot_source}; export graph="wg_immunovar"; \
    export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar; \
    export immune_graph=${graph_base}".subgraph"; export valid_alleles=true;
    export i={}; \
    . ${bigfoot_dir}/filter_immune_subgraph.sh > ${i%.final.bazam.*}_bigfootprint.txt' :::: <(cat ${i});
done
</code>

#### Process samples in chunked parallel threads<br>

<code>ls *.final.bazam.grch38.combined.gam | sort | uniq > process_sample_ids.txt
grep -f 1kgenomes_samples.txt process_sample_ids.txt 
grep -f process_sample_ids.txt 1kgenomes_samples.txt 
split -l 20 process_sample_ids.txt process_sample_split_
for i in $(ls process_sample_split_* );do echo $i;
    time parallel -j 4 'export workdir=${PWD}; export tools_dir=~/tools;
    export PATH=${tools_dir}:$PATH ; \
    export bigfoot_dir=~/tools/BIgFOOT/scripts/; \
    export bigfoot_source=~/pi_kleinstein/bigfoot/; \
    export graphdir=${bigfoot_source}; export graph="wg_immunovar"; \
    export graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar; \
    export immune_graph=${graph_base}".subgraph"; export valid_alleles=true;
    export i={}; \
    . ${bigfoot_dir}/filter_immune_subgraph.sh > ${i%.final.bazam*}_bigfootprint.txt' :::: <(cat ${i});
done
</code>

<i>Assess completed samples - make a list of remaining files to process <br>
<code>grep "All cleaned up!" *_bigfootprint.txt | sed s/":All cleaned up!"//g | sed s/"_bigfootprint.txt"//g > completed_runs.txt
grep -v -f completed_runs.txt run_pipeline_sample_ids.txt > run_pipeline_sample_ids_remaining.txt
</code>

#### Download bam/cram, extract reads, and align in parallel! <br>

<code>cd /media/dduchen/Data/1kgenomes
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
time parallel -j 2 'export workdir=${PWD}; export tools_dir=~/tools;
    export PATH=${tools_dir}:$PATH ; \
    export bigfoot_dir=/home/dduchen/Documents/github/BIgFOOT/scripts;
    export bigfoot_source=/home/dduchen/Documents/bigfoot/;
    export ref=~/tools/refs/grch38_full_wHLA/GRCh38_full_analysis_set_plus_decoy_hla.fa;
    export immunovar_bed=${bigfoot_dir}/../custom_beds/grch38_custom_immunovar_coords.bed;
    export i={}; \
    . ${bigfoot_dir}/download_bam_cram_parse.sh' :::: <(cut -f1 1000G_2504_high_coverage.sequence.index | grep -v "^#");

end with *.combined.gam file for analysis

</code>
<i>To do: <br>
1) Explain default parameters (graph/valid alleles/pe...) <br>
2) Global/local ancestry inference
-- whole genome gam > immunovariation ...sub.gfa > gene graphs<br>
3) Association testing
- Haplotype/allele level association testing script
-- regression analyses in R
- Unitig-caller vs. reference-backbone VCF for association testing - script
-- ls *.sub.gfa > 1kGenomes_immunovar_sub.txt
-- unitig-caller --call --kmer 15 --refs 1kGenomes_immunovar_sub.txt --out 1kGenomes_immunovar_subgraph --write-graph

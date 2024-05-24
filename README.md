# BIgFOOT: Biomarkers of Immunovariation via Graph FOOTprinting
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

<code>bigfoot_source=${tools_dir}/bigfoot # where are we storing all of the reference graph files?</code><br>
<code>mkdir -p ${bigfoot_source} </code><br>
<code>wget -P ${bigfoot_source} "https://zenodo.org/records/10869771/files/immunovar_graph_materials.tar.gz?download=1"</code><br>
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
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507323/SRR507323_2.fastq.gz<br>
- export sample="SRR507323" workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} merged="FALSE" graph="wg_immunovar" valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/preprocess_wg_immunovar_alignment.sh<br>
################################################################</code>

##### Starting from BAM/CRAM (WGS)<br>
- <code>wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989410/NA19240.final.cram
tools_dir=${tools_dir}
PATH=${tools_dir}:$PATH<br>
- export bam_file="NA19240.final.cram" workdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} ref_build="grch38" ref="${bigfoot_source}/GRCh38_full_analysis_set_plus_decoy_hla.fa" merged="FALSE" graph="wg_immunovar" valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/process_from_bam_wg_immunovar_alignment.sh
################################################################</code><br>
<i>Support for CHM13-based BAM/CRAM is planned</i>

##### Starting from subset of reads, some manual pre-processing<br>
- <code>graphdir=${bigfoot_source};graph="wg_immunovar";graph_base=${graphdir}/whole_genome_ig_hla_kir_immunovar;immune_graph=${graph_base}".subgraph";<br>
bazam_reads=${i};
sample_id=${bazam_reads%.bazam.fastq.gz};sample_id=${sample_id##*\/};</code><br>
Sequence-to-graph alignment using VG-giraffe<br>
- <code>vg giraffe -i -f ${bazam_reads} -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.bazam.grch38.wg.gam</code><br>
- <code>vg giraffe -f ${sample_id}.unmapped.fastq.gz -x ${graph_base}.xg -H ${graph_base}.gbwt -d ${graph_base}.dist -m ${graph_base}.min -p > ${sample_id}.unmapped.grch38.wg.gam<br>
cat ${sample_id}.bazam.grch38.wg.gam ${sample_id}.unmapped.grch38.wg.gam > ${sample_id}.bazam.grch38.combined.gam</code><br>
Ready for BIgFOOT
<code>- export i=${sample_id}.bazam.grch38.combined.gam workdir=${PWD} graph=${graph} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} valid_alleles=true<br>
################################################################
. ${bigfoot_dir}/filter_immune_subgraph.sh
################################################################</code><br>

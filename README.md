# BIgFOOT: inference of Biological I(g)mmunovariation via Graph FOOTprinting
### Current version: 0.0.1

This workflow infers alleles, calls novel variation, and constructs sample-specific sequence variation graphs for immunoglobulin(Ig)/other immune-related loci which can be used to perform genetic association tests. The workflow inolves a combination of various graph-construction steps, sequence-to-graph alignment, flow graph decomposition, and unitig calling.<br>
I hope to expand this workflow to enable genome-to-genome analyses/assessing genetic associations between host germline immunovariation and pathogen/metagenomic genetic variation/diversity (i.e., searching for immunological <b>FOOTprints</b> using joint host-pathogen genomic data).

<i>Genetic loci where BIgFOOT performs accurate allele calling:</i>
- IGH
- IGL
- HLA (DQA1/DQB1/... more to come)<br>
<i>Infers alleles - but, like bigoot, I have no proof if they're real (WiP):</i><br>
- IGK
- TR
- KIR

### Input: 
- Raw fastq(.gz)
- BAM/CRAM alignment

## set up conda environment
<b> BIgFOOT is heavily influenced/relies on methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow (v0.0.4)</a>. </b>

1) Clone me! <code> git clone https://github.com/dduchen/BIgFOOT.git </code>
2) set up conda/mamba environment we'll be needing -- can move some of these after the '#' if they're already in your path (e.g., samtools, we assume you have R)<br>
<code>mamba create --name bigfoot -c bioconda -c conda-forge -c gurobi python=3 fastp graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit bbmap minimap2 seqwish blend-bio wfmash samtools pyseer unitig-caller #fastq-dl kmc r-base
conda activate bigfoot </code><br>
Ensure yuo have an active gurobi licence:<br>
<code>gurobi_cl</code><br>
<i> We also use the following R/bioconductor packages: </i><br>
- data.table;
- dplyr;
- Biostrings/DECIPHER </code>

3) We also use some external tools which need to be accessible in your PATH<br>
<code>tools_dir=~/tools;</code> # [wherever you normally install+store software]<br>
<code>PATH=$PATH:${tools_dir};</code><br>
<code>cd ${tools_dir};</code><br>

### download BIgFOOT graph materials from zenodo <a href="https://doi.org/10.5281/zenodo.10674696"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10674696.svg" alt="DOI"></a><br>

<code>bigfoot_source=~/pi_kleinstein/bigfoot</code><br> # where are we planning on storing all of the reference graph files?
<code>bigfoot_source=${tools_dir}/bigfoot</code><br>
<code>mkdir -p ${bigfoot_source} </code><br>

<code>wget -P ${bigfoot_source} "https://zenodo.org/records/10674696/files/immunovar_graph_materials.tar.gz?download=1" </code><br>
<code>cd ${bigfoot_source} ; tar -xvf ${bigfoot_source}/immunovar_graph_materials.tar.gz*</code><br>
#### make distance indexes read only<br>
chmod 0444 *.dist<br>
We also need the VG executable<br>
<code>wget -P ${tools_dir}/ https://github.com/vgteam/vg/releases/download/v1.55.0/vg; chmod +x ${tools_dir}/vg </code><br>
We use the Ryan Wick's Assembly-dereplicator package during haplotype selection <a href="https://github.com/rrwick/Assembly-Dereplicator">Assembly-dereplicator</a>.<br>
- <code> git clone https://github.com/rrwick/Assembly-dereplicator.git ${tools_dir}/Assembly-dereplicator </code><br>
We provide the option of using merged paired-end reads from NGmerge for alignment/inference (optional, not always recommended) <a href="https://github.com/harvardinformatics/NGmerge">NGmerge</a>.<br></code>
- <code>git clone https://github.com/harvardinformatics/NGmerge.git ${tools_dir}/ </code><br>

### Running bigfoot - Example using WES data from ISGR: HG00138
#### British in England and Scotland, European Ancestry

<code>bigfoot_dir=${bigfoot_source}/scripts # Change this if you've downloaded the github repo somewhere else/have the bigfoot analysis scripts saved elsewere<br>
bigfoot_dir=${tools_dir}/BIgFOOT/scripts<br>
immunovar_bed=${bigfoot_source}/grch38_custom_immunovar_coords.bed<br>
conda activate bigfoot<br>
test_dir=${bigfoot_source}/example/<br>
mkdir -p ${test_dir}<br>
cd ${test_dir}<br>
#fastq-dl -a ERR031935 -o ${test_dir}/<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031935/ERR031935_1.fastq.gz<br>
wget -P ${test_dir}/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031935/ERR031935_2.fastq.gz<br> wget -P ${test_dir}/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00138/exome_alignment/HG00138.unmapped.ILLUMINA.bwa.GBR.exome.20120522.bam<br>
wget -P ${test_dir}/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00138/exome_alignment/HG00138.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam<br></code>
samtools merge -o HG00138.bam HG00138.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam HG00138.unmapped.ILLUMINA.bwa.GBR.exome.20120522.bam
rm HG00138.*ILLUMINA*bam

#### Process + Align raw reads<br>
<code>export sample="ERR031935" outdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} merged="FALSE" graph="wg_immunovar"</br>
. ${bigfoot_dir}/preprocess_wg_immunovar_alignment.sh<br></code>

##### from indexed BAM:
<code>export bam_file="HG00138.bam" outdir=${PWD} bigfoot_source=${bigfoot_source} bigfoot_dir=${bigfoot_dir} ref_build="grch38" ref="/home/dd392/tools/refs/annots/GRCh38_full_analysis_set_plus_decoy_hla.fa" merged="FALSE" graph="wg_immunovar"<br>
. ${bigfoot_dir}/process_from_bam_wg_immunovar_alignment.sh

<i>Support for CHM13 available in the next release</i>

####################################################


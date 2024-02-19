# BIgFOOT
## software to infer Biological Immunovariation via Graph FOOTprinting
### Current version: 0.0.1

This workflow infers alleles and calls novel variation via a combination of sequence-to-graph alignment and flow decomposition.

Genetic loci we recommend using BIgFOOT for:
- IGHV
- IGLV
- HLA (DQA1/DQB1/... more to come)
 <i> possible but unvalidated/untuned parameters: </i>
- IGKV
- TR
- KIR

### Input: 
- Raw fastq(.gz)
- BAM/CRAM alignment

## set up conda environment
<b> BIgFOOT is heavily influenced/relies on certain methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow</a>. </b> Development of BIgFOOT based on VG-Flow v0.0.4

1) Clone me! https://github.com/dduchen/BIgFOOT.git
2) set up conda/mamba environment we'll be needing -- can move some of these after the '#' if they're already in your path (e.g., samtools, we assume you have R)
<code> mamba create --name bigfoot -c bioconda -c conda-forge -c gurobi python=3 graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit bbmap minimap2 seqwish blend-bio wfmash samtools pyseer unitig-caller #fastq-dl kmc fastp r-base
conda activate bigfoot </code>
Ensure yuo have an active gurobi licence:
gurobi_cl

<i> We also use the following R/bioconductor packages:
<code> data.table;
dplyr;
Biostrings/DECIPHER </code>

3) We also use some external tools which need to be accessible in your PATH
<code> tools_dir=~/tools; # [wherever you normally install+store software]
PATH=$PATH:${tools_dir};
cd ${tools_dir};
mkdir -p ${tools_dir}/bigfoot </code>
### download igfoot graph materials from zenodo <a href="https://doi.org/10.5281/zenodo.10674696"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10674696.svg" alt="DOI"></a>
<code> wget -P ${tools_dir}/bigfoot/ https://zenodo.org/records/10674696/files/immunovar_graph_materials.tar.gz?download=1 </code>
We also need the VG executable
<code> wget -P ${tools_dir}/ https://github.com/vgteam/vg/releases/download/v1.54.0/vg; chmod +x ${tools_dir}/vg </code>
We use the Ryan Wick's Assembly-dereplicator package during haplotype selection <a href="https://github.com/rrwick/Assembly-Dereplicator">Assembly-dereplicator</a>.
- <code> git clone https://github.com/rrwick/Assembly-dereplicator.git ${tools_dir}/ </code>
We provide the option of using merged paired-end reads from NGmerge for alignment/inference (optional, not always recommended) <a href="https://github.com/harvardinformatics/NGmerge">NGmerge</a>.
- <code> git clone https://github.com/harvardinformatics/NGmerge.git ${tools_dir}/ </code>

### Running bigfoot
- Assuming BAM input
bigfoot_dir=${tools_dir}/bigfoot/scripts
immunovar_bed=${tools_dir}/bigfoot/grch38_custom_immunovar_coords.bed

...


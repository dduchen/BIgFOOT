# Igfoot
### Scalably infer immune-specific allelic variation for genetic association testing 

igfoot is heavily influenced/relies on certain methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow</a>.

A computational method to perform immunoglobulin allele-level genotyping and identification of sample specific variation via flow variation graphs.

Note: This graph-based workflow and scripts currently allow for inference of IG/TR/KIR/HLA alleles from raw fastq sequencing data, BAM/CRAM-based input is possible, the automation of this will be reflected in the next update to this code.

- Current version: 0.0.1

# set up conda environment
Igfoot relies on the excellent VG-Flow method developed to infer the relative abundance of individual viral strains. Development of Igfoot builds upon VG-Flow v0.0.4
1) Please clone the environment from: https://github.com/dduchen/igfoot.git
2) follow the install instructions - with some additional conda packages we'll be needing -- can move some of these after the '#' if they're already in your path (e.g., samtools, we assume you have R)
#mamba create --name igfoot -c bioconda -c conda-forge -c gurobi python=3 graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit bbmap minimap2 seqwish blend-bio wfmash samtools pyseer unitig-caller #fastq-dl kmc fastp r-base
conda activate igfoot
### ensure active gurobi licence:
gurobi_cl
### We also use the following R/bioconductor packages:
--> data.table
--> dplyr
--> Biostrings
--> DECIPHER

3) Download additional tools needed to be installed/accessible in your PATH
tools_dir=~/tools # [choose directory where you normally install/store tools]
PATH=$PATH:${tools_dir}
cd ${tools_dir}
mkdir -p ${tools_dir}/igfoot
### download igfoot graph materials from zenodo:
   [.....]
We need the VG executable
wget https://github.com/vgteam/vg/releases/download/v1.54.0/vg; chmod +x vg
We use the Ryan Wick's Assembly-dereplicator package during haplotype selection <a href="https://github.com/rrwick/Assembly-Dereplicator">Assembly-dereplicator</a>.
- git clone https://github.com/rrwick/Assembly-dereplicator.git
We can use merged paired-end reads from NGmerge for alignment/inference (optional, not necessarily recommended) <a href="https://github.com/harvardinformatics/NGmerge">NGmerge</a>.
- git clone https://github.com/harvardinformatics/NGmerge.git #need to cd NGmerge && make - optional

igfoot_dir=${tools_dir}/igfoot_dir/scripts
immunovar_bed=${tools_dir}/igfoot/grch38_custom_immunovar_coords.bed



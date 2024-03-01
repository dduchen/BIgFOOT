# BIgFOOT: software to infer Biological Immunoglobulins/Immunovariation Graph FOOTprints
### Current version: 0.0.1

This workflow infers alleles, calls novel variation, and constructs sample-specific sequence variation graphs for immunoglobulin(Ig)/other immune-related loci which can be used to perform genetic association tests. The workflow inolves a combination of various graph-construction steps, sequence-to-graph alignment, flow graph decomposition, and unitig calling.<br>
I hope to expand this workflow to enable genome-to-genome analyses/assessing genetic associations between host germline immunovariation and pathogen/metagenomic genetic variation/diversity (i.e., searching for immunological <b>FOOTprints</b> using joint host-pathogen genomic data).

Genetic loci where BIgFOOT performs accurate allele calling:
- IGH
- IGL
- HLA (DQA1/DQB1/... more to come)
 <i>Calls possible - but unvalidated/untuned parameters - WiP:</i>
- IGKV
- TR
- KIR

### Input: 
- Raw fastq(.gz)
- BAM/CRAM alignment

## set up conda environment
<b> BIgFOOT is heavily influenced/relies on certain methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow</a>. </b> Development of BIgFOOT based on VG-Flow v0.0.4

1) Clone me! https://github.com/dduchen/BIgFOOT.git
2) set up conda/mamba environment we'll be needing -- can move some of these after the '#' if they're already in your path (e.g., samtools, we assume you have R)<br>
<code>mamba create --name bigfoot -c bioconda -c conda-forge -c gurobi python=3 graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit bbmap minimap2 seqwish blend-bio wfmash samtools pyseer unitig-caller #fastq-dl kmc fastp r-base
conda activate bigfoot </code><br>
Ensure yuo have an active gurobi licence:<br>
<code>gurobi_cl</code><br>
<i> We also use the following R/bioconductor packages: </i><br>
- data.table;
- dplyr;
- Biostrings/DECIPHER </code>

3) We also use some external tools which need to be accessible in your PATH<br>
<code>tools_dir=~/Documents/tools; # [wherever you normally install+store software]
PATH=$PATH:${tools_dir};
cd ${tools_dir};
mkdir -p ${tools_dir}/BIgFOOT </code><br>
### download BIgFOOT graph materials from zenodo <a href="https://doi.org/10.5281/zenodo.10674696"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10674696.svg" alt="DOI"></a><br>
#### note: this is a large compressed file: ~23GB
<code>wget -P ${tools_dir}/BIgFOOT/ https://zenodo.org/records/10674696/files/immunovar_graph_materials.tar.gz?download=1 </code><br>
We also need the VG executable<br>
<code>wget -P ${tools_dir}/ https://github.com/vgteam/vg/releases/download/v1.55.0/vg; chmod +x ${tools_dir}/vg </code><br>
We use the Ryan Wick's Assembly-dereplicator package during haplotype selection <a href="https://github.com/rrwick/Assembly-Dereplicator">Assembly-dereplicator</a>.<br>
- <code> git clone https://github.com/rrwick/Assembly-dereplicator.git ${tools_dir}/Assembly-dereplicator </code><br>
We provide the option of using merged paired-end reads from NGmerge for alignment/inference (optional, not always recommended) <a href="https://github.com/harvardinformatics/NGmerge">NGmerge</a>.<br>
- <code> git clone https://github.com/harvardinformatics/NGmerge.git ${tools_dir}/ </code><br>

### Running BIgFOOT
- Assuming BAM input
bigfoot_dir=${tools_dir}/BIgFOOT/scripts
immunovar_bed=${tools_dir}/BIgFOOT/grch38_custom_immunovar_coords.bed

...


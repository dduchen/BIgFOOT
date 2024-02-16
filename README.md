# Bigfoot
### Scalably infer immune-specific allelic variation for genetic association testing 

bigfoot is heavily influenced/relies on certain methods developed for <a href="https://bitbucket.org/jbaaijens/vg-flow/src/master/">VG-Flow</a>.

A computational method to perform immunoglobulin allele-level genotyping and identification of sample specific variation via flow variation graphs.

Note: This graph-based workflow and scripts currently allow for inference of IG/TR/KIR/HLA alleles from raw fastq sequencing data, BAM/CRAM-based input is possible, the automation of this will be reflected in the next update to this code.

- Current version: 0.0.1

# set up conda environment
IG flow relies on the excellent VG-Flow method developed to infer the relative abundance of individual viral strains. Development of IG Flow builds upon version: 0.0.4
1) Please clone the environment from: https://bitbucket.org/jbaaijens/vg-flow/src/master/
2) follow the install instructions - with some additional conda packages we'll be needing
	mamba create --name vg-flow-env
	mamba activate vg-flow-env
	mamba install -c bioconda -c conda-forge -c gurobi python=3 graph-tool minimap2 gurobi biopython numpy seqkit odgi bbmap bazam fastp gfaffix fastq-dl kmc blend-bio seqwish wfmash

3) Download the Iw scripts from this repo






# conda environment with everything we need + R and some tools
tools_dir=~/tools # [choose directory where you normally install/store tools]
PATH=$PATH:${tools_dir}
#
cd ${tools_dir}
mkdir -p ${tools_dir}/bigfoot
# you need an executable VG - potentially worth using latest version
wget https://github.com/vgteam/vg/releases/download/v1.54.0/vg; chmod +x vg
# initial setup
#git clone https://bitbucket.org/jbaaijens/vg-flow.git
#git clone https://github.com/rrwick/Assembly-dereplicator
#git clone https://github.com/harvardinformatics/NGmerge.git #need to cd NGmerge && make - optional
# if on AWS / cloud - need to also install gcc
# moving graph files to EC2
#scp -i ~/.ssh/dylangraph.pem /Volumes/DylanDuo/comp_immuno/immunovar_graph_materials.tar.gz dylan@ec2-3-226-244-128.compute-1.amazonaws.com:/home/dylan/tools/ig_flow/
# sign in --> ~/tools/ig_flow
#tar -xvf immunovar_graph_materials.tar.gz
# ig_flow conda environment:
#mamba 
#mamba create --name ig_flow -c bioconda -c conda-forge -c gurobi python=3 graph-tool bazam minimap2 gurobi biopython numpy odgi gfaffix seqkit minimap2 seqwish blend-bio wfmash r-base samtools pyseer unitig-caller
conda activate ig_flow
# ensure active gurobi licence via 'gurobi_cl'
# also should have the following R packages installed:
--> data.table
--> dplyr
--> DECIPHER

vg_flow_dir=${tools_dir}/vg-flow/scripts
#
tools_dir=~/tools
immunovar_bed=${tools_dir}/ig_flow/grch38_custom_immunovar_coords.bed



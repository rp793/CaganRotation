#Scripts are designed for running on the Cambridge CSD3 HPC Cluster using Slurm and an SBATCH system.
#Scripts are optimised for input using a family and/or sample ID and file structure is organised accordingly but is specific to this rotation project's fileset and would need adaptation for future work. 
#Scripts are separated into 3 steps: 
#Step 1 covering processing, alignment to a reference genome and base recalibration with final metrics capture
#Step 2 covering the haplotypecaller step split by chromosome 
#Step 3 covering SNP and Indel calling and hard-filtering step before more metrics capture 
#Step 4 covering denovo annotation, annotation against known variants and additional filtering

#Steps 1,3 and 4 can be run usine the masterscript in each subfolder. Take care to check script version names in case these haven't been updated e.g. v2 vs v3 as they went through several iterations 
#Step 2 can be ran as an individual script for the Haplotypecaller feature

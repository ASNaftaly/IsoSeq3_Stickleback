#PBS -S /bin/bash
#PBS -q batch
#PBS -N randomize_shared_all_tissues
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00
#PBS -M alice.shanfelter@uga.edu
#PBS -m abe

cd /scratch/afs16076/Iso_seq_RNA/Stickleback_Tissues/Compare_Combined_vs_single_datasets/Filtered_Single_Tissue_Analysis/GO_Term_Analyses/

ml Python

for ((n=0;n<10000;n++))
do
  python3 Randomize_GO_terms_isoseq.py Ensembl_98_allproteincoding_genes.csv Genes.Shared.Between.All.Samples.txt 205 Random.GO.Terms.Shared.All.Samples.Distribution.txt ${n}
  echo ${n}
done

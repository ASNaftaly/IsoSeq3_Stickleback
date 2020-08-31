#PBS -S /bin/bash
#PBS -q batch
#PBS -N Pull_40bp_upstream_TSS
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=100:00:00
#PBS -M alice.shanfelter@uga.edu
#PBS -m abe

cd /scratch/afs16076/Iso_seq_RNA/Stickleback_Tissues/TATA_Box_Confirmation/

ml Python

while IFS= read -r line
do
  id=$(echo $line | awk '{print $1}')
  chr=$(echo $line | awk '{print $2}')
  strand=$(echo $line | awk '{print $3}')
  start_pos=$(echo $line | awk '{print $4}')
  echo $id
  python3 PullSeq_Check_TATAbox.py Ensembl_unmasked.fa $id $chr $strand $start_pos All.Isoforms.TSS.plus.40bp.upstream.fasta
done < All.Isoforms.TSS.positions.txt

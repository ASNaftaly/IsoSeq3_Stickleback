#program: SQANTI 
#I was able to load this package locally and just move the scripts and dependencies to a local directory
#other programs needed:
#R 3.4.4; Python 2.7; GMAP 2018
#will run this script before SQANTI QC

#command:
#will need to get gtf from gff from previous steps (collapse redundant isoforms step)
#I used cufflinks 2.2.1 gffread
#will need to remove <> prior to running script

python sqanti_qc.py -g -o <output file prefix> -d <directory for output files> --fl_count <abundance file from step 8> <gtf of isoforms> <known GTF file from reference assembly> <reference assembly in fasta form>

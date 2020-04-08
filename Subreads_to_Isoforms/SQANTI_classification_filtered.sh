#program: SQANTI
#I was able to load this package locally and just move the scripts and dependencies to a local directory
#other programs needed:
#R 3.4.4; Python 2.7; GMAP 2018
#I also had to manually load some python packages and change the export line for those packages
#run this script after FilterSQANTI_input.py

#command:
#will need to get gtf from gff from previous steps (collapse redundant isoforms step)
#will need to remove <> prior to running script

python sqanti_qc.py -g -o <output file prefix> -d <directory for output files> <gtf of filtered isoforms> <known GTF file from reference assembly> <reference assembly in fasta form>

#program: SQANTI
#I was able to load this package locally and just move the scripts and dependencies to a local directory
#other programs needed:
#R 3.4.4; Python 2.7; GMAP 2018
#had to download some R packages and python packages for this script and alter the locations in the python file
#run after SQANTI_classification_Unfiltered.sh

#command:
#will need to remove <> prior to running script

python sqanti_filter.py -i <fasta file from SQANTI QC> <classification file from SQANTI QC>

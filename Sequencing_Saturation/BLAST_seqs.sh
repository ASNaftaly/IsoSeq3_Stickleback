#program: BLAST

#command:
#need to remove <> prior to running script
#-m 9 = tab delimited output format

blastall -p blastn -d <database> -m 9 -a 16 -i <input fasta file> -o <output file>

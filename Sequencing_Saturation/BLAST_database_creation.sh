#program: BLAST

#command:
#remove <> prior to running script

makeblastdb -in <in fasta> -parse_seqids -title "Title" -dbtype nucl -out <output database>

#program: Isoseq3 through conda (v3.1)

#command
#need to remove <> prior to running script
#Primers.fasta = fasta file with primer sequences where headers should be primer_5p or primer_3p (or mutliples of these if you have more than 2 primers)

lima --isoseq --dump-clips --no-pbi -j <# of threads to use> <output.ccs.bam> <Primers.fasta> <output.lima.bam>


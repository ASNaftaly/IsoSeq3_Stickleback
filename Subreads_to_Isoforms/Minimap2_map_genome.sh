#program: minimap2/2.13

#command: 
#need to remove <> prior to running scripts

minimap2 -t <# of threads to use> -ax splice -uf --secondary=no -C5 <reference.fasta> <output.all.hq.fastq> > <output.hq.sam> 2> <output.all.hq.sam.log>

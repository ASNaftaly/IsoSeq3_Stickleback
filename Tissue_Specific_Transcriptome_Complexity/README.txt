#Scripts descriptions and functions

#to remove duplicates from Combined Sexes Analysis (both sexes, all tissues analyzed together)
Script order:  Remove_duplicates_BLAST.py -> Pull_seqs_exons_dups_removed.py

Remove_duplicates_BLAST.py:
	function: identifies isoforms that BLASTed to both itself and another isoform (where the alignment is the same for both and covers a majority of the sequence)
	input files = BLAST results in output format 9, isoseq classification file from SQANTI
	output file = new classification file without the duplicates

Pull_seqs_exons_dups_removed.py:
	function: removes duplicates from the fasta file, faa file, and the gtf file produced from SQANTI (isoseq)
	input files = output file from Remove_duplicates_BLAST.py, nucleotide sequences from fasta file, protein sequences from faa file, exon positions from gtf
	output files = duplicates removed fasta file, protein faa file, exon gtf file


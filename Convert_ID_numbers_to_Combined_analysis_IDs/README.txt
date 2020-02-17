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

#after BLASTing individual tissues to combined sexes database without duplicates
BLAST_form9_pull_fullmatches.py
	function: filters BLAST results based on full matches (isoform exactly matches full length of query (individual tissue) sequences 
	input files = individual tissue classification file from SQANTI, BLAST results from individual tissue to combined sexes database
	output file = blast query matches with full matches marked by "*" and partial matches marked by "-"

#after BLAST_form9_pull_fullmatches.py; needed to further sort BLAST matches so single tissue isoform would only match with one combined tissues isoform 
Sort_best_matches_BLAST.py
	function: excludes BLAST matches except the longest alignment of query to subject, sorts through single tissue isoforms 		that match multiple combined sexes isoforms and keeps the match with the longest alignment
	input files, etc. = <output file from BLAST_form9_pull_fullmatches.py> <individual tissue classification file from SQANTI> <combined tissues classification file from SQANTI> <output sample name, i.e. female.liver.isoform>
	output files = 2 column tab delimited file with column 1 being the isoform ID from the single tissue analysis and column 2 being the isoform ID from the combined tissues analysis
	note: also included the option to expand this script to examine the BLAST results that were excluded from the final file

#Scripts descriptions and functions

#####to remove duplicates from Combined Sexes Analysis (both sexes, all tissues analyzed together)
Script order:  Remove_duplicates_BLAST.py -> Pull_seqs_exons_dups_removed.py

Remove_duplicates_BLAST.py:
	function: identifies isoforms that BLASTed to both itself and another isoform (where the alignment is the same for both and covers a majority of the sequence)
	input files = BLAST results in output format 9, isoseq classification file from SQANTI
	output file = new classification file without the duplicates (combined.sexes.filtered_classification_dupsremoved.txt)

Remove_duplicates_protein_BLAST.py:
	function: runs exactly as Remove_duplicates_BLAST.py, but is for BLAST results run on blastp (protein search); Also collapses isoforms that have the same predicted 	amino acid sequence (these isoforms will be printed as one line)
	input files = BLASTp results in output format 9, combined sexes amino acid sequence file (.faa)
	output files = collapsed isoforms file where each line has each isoform (or set of isoforms) that have a unique predicted protein sequence; removed isoforms where 	each isoform is on a separate line (these were removed if different genes were predicted to have the same amino acid sequence)
	*this script was not used in the final analysis

Pull_seqs_exons_dups_filtered.py:
	function: removes duplicates from the fasta file, faa file, and the gtf file produced from SQANTI (isoseq)
	input files = output file from Remove_duplicates_BLAST.py, nucleotide sequences from fasta file, protein sequences from faa file, exon positions from gtf
	output files = duplicates removed fasta file, protein faa file, exon gtf file (*.filtered_corrected.nodups.fa,fasta, gtf)


#####after BLASTing individual tissues to combined sexes database without duplicates
BLAST_form9_pull_fullmatches.py
	function: filters BLAST results based on full matches (isoform exactly matches full length of query (individual tissue) sequences 
	input files = individual tissue classification file from SQANTI, BLAST results from individual tissue to combined sexes database
	output file = blast query matches with full matches marked by "*" and partial matches marked by "-" (tissue*.BLAST.query.matches.txt)



######after BLAST_form9_pull_fullmatches.py; needed to further sort BLAST matches so single tissue isoform would only match with one combined tissues isoform 
Sort_best_matches_BLAST.py
	function: excludes BLAST matches except the longest alignment of query to subject, sorts through single tissue isoforms 		that match multiple combined sexes isoforms and keeps the match with the longest alignment
	input files, etc. = <output file from BLAST_form9_pull_fullmatches.py> <individual tissue classification file from SQANTI> <combined tissues classification file from SQANTI> <output sample name, i.e. female.liver.isoform>
	output files = 2 column tab delimited file with column 1 being the isoform ID from the single tissue analysis and column 2 being the isoform ID from the combined tissues analysis (tissue*.BLAST.best.matches.combined.sexes.isoforms.txt
	note: also included the option to expand this script to examine the BLAST results that were excluded from the final file


#####adding combined sexes isoform id to classification file
Convert_Classification_File.py
	function: adds combined sexes isoform to column 1 for isoforms in single tissue that matched one combined sexes isoform; also creates a file with all of the individual tissue isoforms that were excluded using the previous scripts
	input files = output from Sort_best_matches_BLAST.py, individual tissue classification file from SQANTI
	output files = new classification file with combined sexes isoforms added, classification file with excluded single tissue isoforms (*filtered.converted.isoforms.classification.txt or *.filtered.excluded.isoforms.classification.txt)
	*can then use Pull_seqs_exons_filtered.py to sort through individual tissue fasta, faa, and gtf files



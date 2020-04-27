#Scripts descriptions and functions

####After running scripts in Convert_ID_numbers_to_Combined_Analysis_IDs, these scripts can be run on single tissue analyses or joined tissue analyses

#To get overall summary counts for each analysis
Gene_Isoform_Counts.py
	function: provides counts for total genes, genes excluded because of mismatches 	between the combined sexes analysis and single/joined tissue analyses, number of 	annotated genes vs movel genes, breakdown of number of isoforms per gene, 		breakdown of number of exons per isoform
	input files = single/joined tissue classification file with converted isoform ID, 	combined sexes classification file with no duplicates
	output files = isoform counts with the format = Gene.ID\tIsoform.Counts
		       exon counts with the format = Isoform.ID\tExon.Counts

#To get counts of isoforms and genes shared between various tissues
#isoforms = Compare_single_tissues.py & Compare_joined_tissues.py
#genes = Compare_single_tissues_genes.py & Compare_joined_tissues_genes.py
#Compare_Single_Tissue_Nts_Proteins.py will likely not be used

Compare_single_tissues.py
	function: comparing isoforms between single tissue analyses; returns isoforms 		shared between the analyses
	input files = single tissue exon counts (output from Gene_Isoform_Counts.py) for 	all 8 single tissues
	output = stdout with counts of total isoforms, counts of shared isoforms between 	all tissues, all somatic tissues, all female tissues, all male tissues, liver, 		brain, pronephros, gonads, liver and brain, liver and pronephros, liver and 		gonads, brain and pronephros, brain and gonads, pronephros and gonads
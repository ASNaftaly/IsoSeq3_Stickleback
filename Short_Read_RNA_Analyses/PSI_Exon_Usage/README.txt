#Script descriptions and functions

####To create exon trio database to calculate PSI
#used Ensembl 97 build and isoseq combined sexes to create trios for this case

PSI_ExonUsage_Database_Creator.py
	function = combines ensembl and isoseq isoforms, collapses duplicates or similar isoforms as needed, creates exon trios, filters these trios for duplicates and collapses differences in start position of first exon or end of last exon, writes final output
	input files = isoseq classification file, isoseq gtf file, ensembl gtf file
	output file = tab delimited file with format:
	#Gene	Strand	Exon.1.Start	Exon.1.End	Exon.2.Start	Exon.2.End	Exon.3.Start	Exon.3.End


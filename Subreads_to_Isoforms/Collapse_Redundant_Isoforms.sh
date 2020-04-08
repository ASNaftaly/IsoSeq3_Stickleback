#program: cDNA_cupcake

#command:
#need to remove <> prior to running script
#collapse_isoforms_by_sam.py = source: https://github.com/Magdoll/cDNA_Cupcake/wiki#refminimap
#will need to download cDNA cupcake scripts/package

python collapse_isoforms_by_sam.py --input <polished hq fastq> --fq -s <sorted sam file> -o <output prefix>

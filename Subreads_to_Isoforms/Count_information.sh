#program = cDNA Cupcake

#command:
#need to remove <> prior to running script
#get_abundance_post_collapse.py = need to download from https://github.com/Magdoll/cDNA_Cupcake/wiki#refminimap
#cluster report = created in polish step not cluster (as that would make sense)

python get_abundance_post_collapse.py <prefix for output collaspe_redundant_isoforms.sh> <cluster report.csv>

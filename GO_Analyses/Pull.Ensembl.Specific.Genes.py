#Pulling genes from ensembl that are not represented in the Isoseq transcriptome
#will read in isoseq full transcriptome from classification file and the gene ids from Ensembl.Genes.GO.terms.txt
#note: this will only examine genes that have GO terms in ensembl (which is true for all of the GO analyses)
#to run script: python3 Pull.Ensembl.Specific.Genes.py <classification file> <ensembl genes GO terms file> <output file with 1 gene per line)
#Author: Alice Naftaly, May 2020


import sys


#read in classification file
#returns list of gene ids
def read_class():
    class_file = sys.argv[1]
    isoseq_genes = []
    final_isoseq_genes = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                gene_id = new_line[6]
                isoseq_genes.append(gene_id)
    final_isoseq_genes = list(set(isoseq_genes))
    return final_isoseq_genes


#read in ensembl gene ids
#file format: Gene.ID \t GO.terms (in list separated by comma)
#returns list of gene ids
def read_ensembl_genes():
    gene_file = sys.argv[2]
    ensembl_genes = []
    with open(gene_file, 'r') as genes:
        for line in genes:
            new_line = line.split()
            gene_id = new_line[0]
            ensembl_genes.append(gene_id)
    return ensembl_genes


#filter for genes found in ensembl, but not isoseq
#rreturns list of ensembl specific genes
def filter_genes():
    isoseq_genes = read_class()
    ensembl_genes = read_ensembl_genes()
    ensembl_specific_genes = []
    for gene in ensembl_genes:
        if gene not in isoseq_genes:
            ensembl_specific_genes.append(gene)
    return ensembl_specific_genes


#write ensembl specific genes to a file
#output format: 1 gene id per line
def write():
    ensembl_specific_genes = filter_genes()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for gene in ensembl_specific_genes:
            final = "%s\n" % str(gene)
            out.write(final)

write()

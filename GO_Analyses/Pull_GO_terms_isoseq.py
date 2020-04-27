#Pull GO terms for ensembl annotated genes
#need input files that list gene IDs 1 per line, csv or xlsx file from ensembl with go terms, genes, etc.
#will pull all go terms from each gene in file, if no go terms are available, will skip that gene
#final file should have the final format:
#Gene \t list of GO terms separated by ,
#to run script: Pull_GO_terms_isoseq.py <Gene id file> <ensembl table (csv or xlsx) with go terms associated with genes> <output file>
#Author: Alice Naftaly, March 2020

import sys

#reading gene ids from file
#returns list of gene IDs
#skips all novel genes
def read_gene_ids():
    gene_id_file = sys.argv[1]
    all_gene_ids = []
    with open(gene_id_file, 'r') as genes:
        for line in genes:
            if line.startswith("ENSG"):
                all_gene_ids.append(line.strip())
    return all_gene_ids

#read in ensembl information
#only need gene id and go terms
#returns dictionary with key = gene id and value == list of go terms
def read_ensembl_table():
    ensembl_table_file = sys.argv[2]
    go_term_dict = {}
    final_go_terms = {}
    with open(ensembl_table_file, 'r') as ensembl_table:
        for line in ensembl_table:
            new_line = line.split(",")
            gene_id = new_line[0]
            go_term = new_line[len(new_line)-1].strip("\n")
            if gene_id in go_term_dict:
                go_term_dict[gene_id].append(go_term)
            elif gene_id not in go_term_dict:
                go_term_dict.update({gene_id:[go_term]})
    #remove genes with no go terms or empty go term values
    for key in go_term_dict:
        single_key = go_term_dict[key]
        list_of_go_terms = []
        for value in single_key:
            if value.startswith("GO:"):
                list_of_go_terms.append(value)
        if len(list_of_go_terms) > 0:
            final_go_terms.update({key:list_of_go_terms})
    return final_go_terms


#create dictionary with go terms for the genes present in shared genes file
#returns dictionary with key = gene id and value == list of go terms
def match_genes_go_terms():
    shared_genes = read_gene_ids()
    go_terms = read_ensembl_table()
    shared_go_terms = {}
    for gene in shared_genes:
        if gene in go_terms:
            single_go_terms = go_terms[gene]
            shared_go_terms.update({gene:single_go_terms})
    return shared_go_terms

#write go terms to file with gene ids
def write():
    output = sys.argv[3]
    genes_and_go_terms = match_genes_go_terms()
    with open(output, 'a') as out:
        for gene in genes_and_go_terms:
            single_gene = genes_and_go_terms[gene]
            final_value = ",".join(single_gene)
            to_write = "%s\t%s\n" % (str(gene), final_value)
            out.write(to_write)

write()

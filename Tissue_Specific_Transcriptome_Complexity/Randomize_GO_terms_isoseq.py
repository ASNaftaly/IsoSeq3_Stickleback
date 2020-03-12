#compare go terms from shared gene subsets to random gene selections
#this is the script to create randomly selected genes and go terms
#1. randomly pull specific number of genes from all genes list and pull those GO terms (only choose genes that have annotated GO terms)
#2. repeat 10,000 times and write summary file for each iteration (bash script)
#will compare this distribution to go terms from shared gene subsets
#to run script: python3 Randomize_GO_terms_isoseq.py <ensembl table as csv> <shared genes file> <number of genes to randomize; integer> <output file> <iteration number>

import sys
import random

#first pull genes and go terms from ensembl table
#only need gene id and go terms
#returns dictionary with key = gene id and value == list of go terms
def read_ensembl_table():
    ensembl_table_file = sys.argv[1]
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


#pull genes from shared tissue genes
#need to remove these from the randomization
#returns a list with the shared genes
def pull_shared_genes():
    shared_gene_file = sys.argv[2]
    gene_list = []
    with open(shared_gene_file, 'r') as genes:
        for line in genes:
            if line.startswith("ENSGACG"):
                gene_list.append(line.strip())
    return gene_list

#removing shared genes from all genes list
def remove_genes():
    all_genes = read_ensembl_table()
    shared_genes = pull_shared_genes()
    for gene in shared_genes:
        if gene in all_genes.keys():
            del all_genes[gene]
    return all_genes

#randomly select a subset of genes based on the number of genes in the shared gene set
def randomize_genes():
    all_genes = remove_genes()
    gene_names = list(all_genes.keys())
    randomize_total = int(sys.argv[3])
    go_dict = {}
    #random.sample is sampling without replacement
    random_genes = random.sample(gene_names, randomize_total)
    for gene in random_genes:
        random_go_terms = all_genes[gene]
        #now iterate through go terms
        #creates dictionary with key == go term, and value == number of instances that go term was seen
        for value in random_go_terms:
            dict_value = 1
            if value in go_dict:
                go_dict[value].append(dict_value)
            elif value not in go_dict:
                go_dict.update({value:[dict_value]})
    return go_dict

#count occurrences of GO terms from dictionary
#writes go term counts to file
def count_GO():
    go_dict = randomize_genes()
    output = sys.argv[4]
    iteration_number = sys.argv[5]
    with open(output,'a') as out:
        for key in go_dict:
            counts = len(go_dict[key])
            final = "%s\t%s\t%s\n" % (str(iteration_number), str(key), str(counts))
            out.write(final)

count_GO()

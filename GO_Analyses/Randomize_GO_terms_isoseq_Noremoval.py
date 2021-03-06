#compare go terms from shared gene subsets to random gene selections
#this is the script to create randomly selected genes and go terms
#this script will not remove gene set being examined
#1. randomly pull specific number of genes from all genes list and pull those GO terms (only choose genes that have annotated GO terms)
#2. repeat 10,000 times and write summary file for each iteration (bash script)
#will compare this distribution to go terms from shared gene subsets
#to run script: python3 Randomize_GO_terms_isoseq_Noremoval.py <Ensembl genes and GO terms file in format of from Combined_Ensembl_Novel_GO_terms.py> <number of genes to randomize; integer> <output file> <iteration number>
#Author: Alice Naftaly, April 2020

import sys
import random

#read GO terms from combined and novel GO terms file
#returns dictionary with key == gene and value == list of go terms
def read_GO_terms():
    all_go_terms = sys.argv[1]
    go_term_dict = {}
    with open(all_go_terms, 'r') as go_terms:
        for line in go_terms:
            new_line = line.split("\t")
            gene = new_line[0]
            stripped_terms = new_line[1].strip("\n")
            list_of_go_terms = stripped_terms.split(",")
            go_term_dict.update({gene:list_of_go_terms})
    return go_term_dict


#randomly select a subset of genes based on the number of genes in the shared gene set
def randomize_genes():
    all_genes = read_GO_terms()
    gene_names = list(all_genes.keys())
    randomize_total = int(sys.argv[2])
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
    output = sys.argv[3]
    iteration_number = sys.argv[4]
    with open(output,'a') as out:
        for key in go_dict:
            counts = len(go_dict[key])
            final = "%s\t%s\t%s\n" % (str(iteration_number), str(key), str(counts))
            out.write(final)

count_GO()

#converting GO term files into WEGO input format
#WEGO input format = gene id/transcript id \t GO terms separated by tab
#will need to read in permutations distribution (before and after adjusting p values), file with gene ids or transcript ids with all GO terms for each gene/transcript (in this case the output from Combined_Ensembl_Novel_GO_terms.py gene or transcript version)
#to run script: python3 Make_WEGO_input_format.py <permutations distribution> <output from Combined_Ensembl_Novel_GO_terms to use as input> <final output file in WEGO format>
#Author: Alice Naftaly, April 2020

import sys


#read permutations output
#returns list of go terms
def read_permutations():
    permutations_file = sys.argv[1]
    GO_terms = []
    with open(permutations_file, 'r') as permutations:
        for line in permutations:
            if line.startswith("GO:"):
                new_line = line.split("\t")
                GO_terms.append(new_line[0])
    return GO_terms


#read id files
#returns dictionary with key == GO term and value = list of isoforms with that GO term
#will have to resort this dictionary so key == isoform after filtering with permutations output
def read_ids():
    ids_file = sys.argv[2]
    ids_dict = {}
    with open(ids_file, 'r') as ids:
        for line in ids:
            new_line = line.split()
            isoform_id = new_line[0]
            split_go_terms = new_line[1].split(",")
            for term in split_go_terms:
                if term in ids_dict:
                    ids_dict[term].append(isoform_id)
                elif term not in ids_dict:
                    ids_dict.update({term:[isoform_id]})
    return ids_dict

#sort ids_dict by removing any GO terms not in permutations file
#this is more important for the enriched sets that the observed sets
#returns dictionary with key == GO term and value = all isoforms with that GO term
def filter_ids():
    ids_dict = read_ids()
    go_terms = read_permutations()
    filtered_ids_dict = {}
    for term in go_terms:
        if term in ids_dict:
            single_term = ids_dict[term]
            filtered_ids_dict.update({term:single_term})
    return filtered_ids_dict

#create WEGO format dictionary
#returns a dictionary with key == isoform id and value == GO terms for that isoform
def create_WEGO_format():
    filtered_ids = filter_ids()
    final_dict = {}
    for go_term in filtered_ids:
        single_go_term = filtered_ids[go_term]
        for isoform in single_go_term:
            if isoform in final_dict:
                final_dict[isoform].append(go_term)
            elif isoform not in final_dict:
                final_dict.update({isoform:[go_term]})
    return final_dict


#write WEGO output format
def write_WEGO():
    ids_dict = create_WEGO_format()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for isoform in ids_dict:
            go_terms = ids_dict[isoform]
            final_go_terms = "\t".join(go_terms)
            final = isoform + "\t" + final_go_terms + "\n"
            out.write(final)

write_WEGO()

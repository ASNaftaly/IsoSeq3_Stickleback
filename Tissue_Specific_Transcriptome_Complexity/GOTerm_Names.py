#need to pull GO term names from significant GO terms
#will use significant go terms file for input and csv with all ensembl go terms (If one of the novel GO terms isn't in this list, will pull manually)
#to run script: python3 GOTerm_Names.py <sig go terms file> <tab delimited file with ensembl go terms> <output file; format = GO.Term\tGO.name>
#author: Alice Naftaly, March 2020

import sys

#read in significant go terms file
#returns list of GO terms
def read_sig_go_terms():
    sig_file = sys.argv[1]
    go_terms = []
    with open(sig_file, 'r') as sig:
        for line in sig:
            if line.startswith("GO:"):
                new_line = line.split()
                go_terms.append(new_line[0])
    return go_terms


#read ensembl csv file
#Need GO terms and go term names
#returns dictionary with key == go term and value == go term name
def read_ensembl():
    ensembl_file = sys.argv[2]
    ensembl_dict = {}
    with open(ensembl_file, 'r') as ensembl:
        for line in ensembl:
            new_line = line.split("\t")
            go_term = new_line[6].strip("\n")
            go_term_name = new_line[4]
            if len(go_term) > 0:
                ensembl_dict.update({go_term:go_term_name})
    return ensembl_dict


#pull sig go term names
#if significant go term is in the ensembl file, it will pull the name of the GO term
#if not, just the GO term will be added and I will manually fill in the name later
#returns dictionary with key = go term and value == go term name
def pull_names():
    sig_go_terms = read_sig_go_terms()
    ensembl = read_ensembl()
    final_dict = {}
    for go_term in sig_go_terms:
        if go_term in ensembl:
            single_ensembl = ensembl[go_term]
            final_dict.update({go_term:single_ensembl})
        else:
            final_dict.update({go_term:""})
    return final_dict

#write go names to file
def write():
    go_names = pull_names()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "GO.Term\tGO.Name\n"
        out.write(header)
        for go_term in go_names:
            single_name = go_names[go_term]
            final = "%s\t%s\n" % (str(go_term), str(single_name))
            out.write(final)

write()

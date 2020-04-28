#pull GO terms from interproscan to analyze in revigo and other GO term analyses
#want to get number of isoforms that do not have GO terms, no protein match, and ones with GO terms
#reads in tsv output from interproscan
#to run script: python3 Pull_GO_terms_Interproscan_Output.py <tsv file from interproscan> <output no GO terms> <output GO terms>
#Author: Alice Naftaly, March 2020, edited April 2020

import sys

#read in tsv file from interproscan
#returns dictionary with key = isoform id and value = go terms
def read_interproscan():
    interproscan_file = sys.argv[1]
    go_term_dict = {}
    final_go_terms = {}
    with open(interproscan_file, 'r') as interproscan:
        for line in interproscan:
            new_line = line.split("\t")
            isoform_id = new_line[0]
            if len(new_line) == 11:
                go_terms = "n/a"
                if isoform_id in go_term_dict:
                    go_term_dict[isoform_id].append(go_terms)
                elif isoform_id not in go_term_dict:
                    go_term_dict.update({isoform_id:[go_terms]})
            elif len(new_line) == 15:
                go_values = new_line[13]
                if go_values.startswith("GO:"):
                    if len(go_values) == 10:
                        go_terms = go_values
                    elif len(go_values) > 10:
                        go_terms = go_values.split("|")
                    if isoform_id in go_term_dict:
                        go_term_dict[isoform_id].append(go_terms)
                    elif isoform_id not in go_term_dict:
                        go_term_dict.update({isoform_id:[go_terms]})
        for key in go_term_dict:
            single_key = go_term_dict[key]
            all_go_terms = []
            for term in single_key:
                if term == "n/a":
                    all_go_terms.append(term)
                elif len(term) == 10 and isinstance(term, list) == False:
                    all_go_terms.append(term)
                else:
                    for single in term:
                        all_go_terms.append(single)
            set_go_terms = list(set(all_go_terms))
            final_go_terms.update({key:set_go_terms})
    return final_go_terms


#sort out isoforms with no go terms
#will print summary to stdout
#also returns lsit of isoforms with no go terms
def sort_no_go_terms():
    go_terms = read_interproscan()
    no_go_terms = []
    for isoform in go_terms:
        single_isoform = go_terms[isoform]
        term_count = 0
        for term in single_isoform:
            if term == "n/a":
                term_count += 1
        if term_count == len(single_isoform):
            no_go_terms.append(isoform)
    print("Number of Isoforms with no GO terms")
    print(len(no_go_terms))
    return no_go_terms


#sorts through isoforms with GO terms:
#returns dictionary with key = isoform and value ==
#there can be multiple go terms per isoform
def sort_go_terms():
    go_terms = read_interproscan()
    isoforms_with_go_terms = {}
    go_terms_count = []
    for isoform in go_terms:
        single_isoform = go_terms[isoform]
        final_go_terms = []
        for term in single_isoform:
            if term.startswith("GO:"):
                final_go_terms.append(term)
                go_terms_count.append(term)
        if len(final_go_terms) > 0:
            isoforms_with_go_terms.update({isoform:final_go_terms})
    print("Number of isoforms with GO terms")
    print(len(isoforms_with_go_terms))
    print("Number of GO terms identified")
    print(len(set(go_terms_count)))
    return isoforms_with_go_terms


#no go terms output format = 1 column with header (Isoform.ID), one isoform per line
def write_no_GO_terms():
    no_go_terms = sort_no_go_terms()
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "Isoform.ID"
        out.write(header)
        for value in no_go_terms:
            final = "%s\n" % str(value)
            out.write(final)

#isoforms with go terms output format = 2 columns with header, 1 = Isoform.ID, 3 = GO term
def write_GO_terms():
    go_terms = sort_go_terms()
    output = sys.argv[3]
    with open(output, 'w') as out:
        header = "Isoform.ID\tGO.Term\n"
        out.write(header)
        for isoform in go_terms:
            single_isoform = go_terms[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                final = "%s\t%s\n" % (str(isoform),str(single))
                out.write(final)
            elif len(single_isoform) > 1:
                all_go_terms = ",".join(single_isoform)
                final = "%s\t%s\n" % (str(isoform),all_go_terms)
                out.write(final)


#call all functions
def call():
    no_go_terms = write_no_GO_terms()
    go_terms = write_GO_terms()

call()

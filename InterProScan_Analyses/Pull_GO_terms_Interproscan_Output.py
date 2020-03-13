#pull GO terms from interproscan to analyze in revigo and other GO term analyses
#want to get number of isoforms that do not have GO terms, no protein match, and ones with GO terms
#1. read BLAST2GO table (created from loading interproscan results into BLAST2GO and then exporting table)


#to run script: python3 Pull_GO_terms_Interproscan_Output.py

import sys

#read in table from BLAST2GO
#returns dictionary with key = isoform id and value = go terms
def read_table():
    input_file = sys.argv[1]
    go_term_dict = {}
    with open(input_file, 'r') as blast2go_table:
        for line in blast2go_table:
            new_line = line.split("\t")
            isoform_id = new_line[2]
            go_terms = new_line[14]
            if isoform_id in go_term_dict:
                go_term_dict[isoform_id].append(go_terms)
            elif isoform_id not in go_term_dict:
                go_term_dict.update({isoform_id:[go_terms]})
    return go_term_dict

#sort out isoforms with no protein matches through interproscan
#will print summary to stdout
#also returns list of isoforms with no matches
def sort_no_matches():
    go_terms = read_table()
    no_protein_matches = []
    for isoform in go_terms:
        single_isoform = go_terms[isoform]
        if single_isoform[0] == "no IPS match":
            no_protein_matches.append(isoform)
    print("Number of Isoforms with no Interproscan match")
    print(len(no_protein_matches))
    return no_protein_matches


#sort out isoforms with no go terms
#will print summary to stdout
#also returns lsit of isoforms with no go terms
def sort_no_go_terms():
    go_terms = read_table()
    no_go_terms = []
    for isoform in go_terms:
        single_isoform = go_terms[isoform]
        if single_isoform[0] == "no GO terms":
            no_go_terms.append(isoform)
    print("Number of Isoforms with no GO terms")
    print(len(no_go_terms))
    return no_go_terms


#sorts through isoforms with GO terms:
#returns dictionary with key = isoform and value == [ category (F, C, P)*I think these refer to molecular function (F), cellular component (C), biological process (P), GO term]
#there can be multiple go terms per isoform
def sort_go_terms():
    go_terms = read_table()
    isoforms_with_go_terms = {}
    go_terms_count = []
    for isoform in go_terms:
        single_isoform = go_terms[isoform]
        if single_isoform[0] != "no IPS match" and single_isoform[0] != "no GO terms" and single_isoform[0] != "InterPro GO IDs":
            split_isoform = single_isoform[0].split(";")
            for value in split_isoform:
                no_spaces_value = value.strip(" ")
                split_value = no_spaces_value.split(":")
                final_go_term = split_value[1] + ":" + split_value[2]
                go_terms_count.append(final_go_term)
                final = [split_value[0], final_go_term]
                if isoform in isoforms_with_go_terms:
                    isoforms_with_go_terms[isoform].append(final)
                elif isoform not in isoforms_with_go_terms:
                    isoforms_with_go_terms.update({isoform:[final]})
    print("Number of isoforms with GO terms")
    print(len(isoforms_with_go_terms))
    print("Number of GO terms identified")
    print(len(set(go_terms_count)))
    return isoforms_with_go_terms


#Writing output for no IPS match, no go terms, and isoforms with go terms
#no matches output format = 1 column with header (Isoform.ID), one isoform per line
def write_no_matches():
    no_matches = sort_no_matches()
    output = sys.argv[2]
    with open(output, 'w') as out:
        header = "Isoform.ID"
        out.write(header)
        for value in no_matches:
            final = "%s\n" % str(value)
            out.write(final)

#no go terms output format = 1 column with header (Isoform.ID), one isoform per line
def write_no_GO_terms():
    no_go_terms = sort_no_go_terms()
    output = sys.argv[3]
    with open(output, 'w') as out:
        header = "Isoform.ID"
        out.write(header)
        for value in no_go_terms:
            final = "%s\n" % str(value)
            out.write(final)

#isoforms with go terms output format = 3 columns with header, 1 = Isoform.ID, 2 = Catgeory, 3 = GO term
def write_GO_terms():
    go_terms = sort_go_terms()
    output = sys.argv[4]
    with open(output, 'w') as out:
        header = "Isoform.ID\tCategory\tGO.Term\n"
        out.write(header)
        for isoform in go_terms:
            single_isoform = go_terms[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                final = "%s\t%s\t%s\n" % (str(isoform),str(single[0]),str(single[1]))
                out.write(final)
            elif len(single_isoform) > 1:
                for value in single_isoform:
                    final = "%s\t%s\t%s\n" % (str(isoform),str(value[0]),str(value[1]))
                    out.write(final)


#call all functions
def call():
    no_matches = write_no_matches()
    no_go_terms = write_no_GO_terms()
    go_terms = write_GO_terms()

call()

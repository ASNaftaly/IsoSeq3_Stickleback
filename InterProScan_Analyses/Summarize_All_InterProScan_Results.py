#summarizing interproscan output for all novel isoforms
#will use tsv output from interproscan as input
#to run script: python Summarize_All_InterProScan_Results.py <tsv from interproscan>

import sys

#read in tsv file
#columns that are needed:
#isoform id (0), sequence length (2), analysis (3), signature assession (4), signature description (5), start location (6), end location (7), score (8), interproscan accession if available (11), interproscan description if available (12), GO term if available (13), Pathway info if available (14)
def read_tsv():
    tsv_file = sys.argv[1]
    domain_dict = {}
    with open(tsv_file, 'r') as tsv:
        for line in tsv:
            new_line = line.split("\t")
            #if line is only 11 values long, this means all of the optional columns (12-15) are not present)
            if len(new_line) == 11:
                isoform = new_line[0]
                seq_length = new_line[2]
                program = new_line[3]
                signature_identifier = new_line[4]
                signature_description = new_line[5]
                start_location = new_line[6]
                end_location = new_line[7]
                score = new_line[8]
                dict_value = [seq_length, program, signature_identifier, signature_description, start_location, end_location, score]
                if isoform in domain_dict:
                    domain_dict[isoform].append(dict_value)
                elif isoform not in domain_dict:
                    domain_dict.update({isoform:[dict_value]})
            elif len(new_line) == 15:
                isoform = new_line[0]
                seq_length = new_line[2]
                program = new_line[3]
                signature_identifier = new_line[4]
                signature_description = new_line[5]
                start_location = new_line[6]
                end_location = new_line[7]
                score = new_line[8]
                interpro_accession = new_line[11]
                interpro_description = new_line[12]
                go_term = new_line[13]
                pathway_info = new_line[14]
                dict_value = [seq_length, program, signature_identifier, signature_description, start_location, end_location, score, interpro_accession, interpro_description, go_term, pathway_info]
                if isoform in domain_dict:
                    domain_dict[isoform].append(dict_value)
                elif isoform not in domain_dict:
                    domain_dict.update({isoform:[dict_value]})
    return domain_dict


#examining the signature description
def examine_descriptions():
    domain_dict = read_tsv()
    for isoform in domain_dict:
        single_isoform = domain_dict[isoform]
        for entry in single_isoform:
            print(entry)


examine_descriptions()

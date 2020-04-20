#summarizing interproscan output for all novel isoforms
#basically will be trying to recreate IPRStats output as this program will not work for me
#will use tsv output from interproscan as input
#will separate out results based on program (analysis column 3 in tsv), will write one output with each of the program results then a log file with the programs detected and the number of values from each program
#to run script: python Summarize_All_InterProScan_Results.py <tsv from interproscan> <output> >> <log file>
#Author: Alice Naftaly, April 2020

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
                dict_value = [isoform, seq_length, program, signature_identifier, signature_description, start_location, end_location, score]
                if program in domain_dict:
                    domain_dict[program].append(dict_value)
                elif program not in domain_dict:
                    domain_dict.update({program:[dict_value]})
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
                dict_value = [isoform, seq_length, program, signature_identifier, signature_description, start_location, end_location, score, interpro_accession, interpro_description, go_term, pathway_info]
                if program in domain_dict:
                    domain_dict[program].append(dict_value)
                elif program not in domain_dict:
                    domain_dict.update({program:[dict_value]})
    return domain_dict

#count predictions per program
#
def write_stats():
    domain_dict = read_tsv()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for program in domain_dict:
            single_program = domain_dict[program]
            program_header = "Number of predictions from %s database: %d" % (str(program), len(single_program))
            print(program_header)
            header = "Program\tIsoform.ID\tIdentifier\tDescription\tInterProScan.Accession\tInterProScan.Description\tGO.Term\tScore\n"
            out.write(header)
            for entry in single_program:
                if len(entry) == 8:
                    isoform_id = entry[0]
                    sig_identifier = entry[3]
                    sig_description = entry[4]
                    score = entry[7]
                    final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(program), str(isoform_id), str(sig_identifier), str(sig_description), ".", ".", ".",str(score))
                    out.write(final)
                elif len(entry) == 12:
                    isoform_id = entry[0]
                    sig_identifier = entry[3]
                    sig_description = entry[4]
                    score = entry[7]
                    interproscan_accession = entry[8]
                    interproscan_description = entry[9]
                    go_terms = entry[10]
                    if go_terms.startswith("GO"):
                        continue
                    else:
                        go_terms = "."
                    final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(program), str(isoform_id), str(sig_identifier), str(sig_description), str(interproscan_accession), str(interproscan_description), str(go_terms),str(score))
                    out.write(final)

write_stats()

#write file for R
def write_stats_R():
    domain_dict = read_tsv()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for program in domain_dict:
            single_program = domain_dict[program]
            header = "Program\tIsoform.ID\tIdentifier\tInterProScan.Accession\tScore\n"
            out.write(header)
            for entry in single_program:
                if len(entry) == 8:
                    isoform_id = entry[0]
                    sig_identifier = entry[3]
                    sig_description = entry[4]
                    score = entry[7]
                    final = "%s\t%s\t%s\t%s\t%s\n" % (str(program), str(isoform_id), str(sig_identifier), ".",str(score))
                    out.write(final)
                elif len(entry) == 12:
                    isoform_id = entry[0]
                    sig_identifier = entry[3]
                    sig_description = entry[4]
                    score = entry[7]
                    interproscan_accession = entry[8]
                    interproscan_description = entry[9]
                    go_terms = entry[10]
                    if go_terms.startswith("GO"):
                        continue
                    else:
                        go_terms = "."
                    final = "%s\t%s\t%s\t%s\t%s\n" % (str(program), str(isoform_id), str(sig_identifier), str(interproscan_accession),str(score))
                    out.write(final)

write_stats_R()

#removing duplicates identified from BLAST results in combined sexes isoseq data
#if query start and end matches subject start and end and only matches to itself = perfect matches
#to run script: python3 Remove_duplicates_BLAST.py <blast results for combined sexes> <combined sexes classification file> <output new classification file without duplicates>
#Author: Alice Naftaly, February 2020

import sys

#read in blast output
#returns dictionary with key = isoform id and value = alignment
def read_blast():
    blast_output = sys.argv[1]
    blast_dict = {}
    with open(blast_output, 'r') as blast_results:
        for line in blast_results:
            if line.startswith("#"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                if isoform_id in blast_dict:
                    blast_dict[isoform_id].append(new_line)
                elif isoform_id not in blast_dict:
                    blast_dict.update({isoform_id:[new_line]})
    return blast_dict

#need the length of combined sexes isoforms
#will create a dictionary with key = isoform id and value = length of isoform
def pull_combined_isoform():
    combined_isoform_file = sys.argv[2]
    combined_isoform_dict = {}
    with open(combined_isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                length = new_line[3]
                combined_isoform_dict.update({isoform_id:length})
    return combined_isoform_dict

#create dictionary for classification file so duplicates can be removed
#key = first value (isoform for first line as header) and isoform id for the rest of the lines and the value == the line
def create_class_dict():
    class_file = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_input:
        for line in class_input:
            new_line = line.split()
            key_value = new_line[0]
            class_dict.update({key_value:line})
    return class_dict


#sorting through blast results:

#if blast matches are perfect between query and subject (should be only 1 isoform matching itself best)
def sort_matches():
    blast_results = read_blast()
    isoform_lengths = pull_combined_isoform()
    full_matches = {}
    partial_matches = {}
    duplicates = []
    for isoform in blast_results:
        single_isoform = blast_results[isoform]
        isoform_length = isoform_lengths[isoform]
        if len(single_isoform) == 1:
            full_matches.update({isoform:single_isoform})
        elif len(single_isoform) > 1:
            for single in single_isoform:
                q_isoform_id = single[0]
                s_isoform_id = single[1]
                q_isoform_length = single[3]
                q_start = single[6]
                q_end = single[7]
                s_start = single[8]
                s_end = single[9]
                #identifies full matches where the ids are the same and the length aligned matches the sequence length from the classification file
                if q_isoform_id == s_isoform_id and int(isoform_length) == int(q_isoform_length) and q_start == s_start and q_end == s_end:
                    full_matches.update({isoform:single})
                #partial matches cover most of the sequence (either starting at 1 and finishing at the same length, but skipping some nucleotides in the middle of the sequence)
                #perhaps this is because of repetitive regions?
                elif q_isoform_id == s_isoform_id and int(isoform_length) != int(q_isoform_length) and q_start == s_start and q_end == s_end:
                    if isoform in partial_matches:
                        partial_matches[isoform].append(single)
                    elif isoform not in partial_matches:
                        partial_matches.update({isoform:[single]})
                #identifies duplicates (where different isoforms match perfectly); will need to remove these isoforms from full matches and partial matches
                elif q_isoform_id != s_isoform_id and int(q_isoform_length) == int(isoform_lengths[s_isoform_id]) and q_start == s_start and q_end == s_end:
                    duplicates.append(q_isoform_id)
                    duplicates.append(s_isoform_id)
    for iso in duplicates:
        if iso in full_matches:
            del full_matches[iso]
        if iso in partial_matches:
            del partial_matches[iso]
    set_duplicates = list(set(duplicates))
    missing_isoforms = []
    for i in isoform_lengths:
        if i not in full_matches and i not in partial_matches and i not in set_duplicates and i != "isoform":
            missing_isoforms.append(i)
    return full_matches, partial_matches, set_duplicates, missing_isoforms

#creating classification file with duplicates removed
def sort_class_dict():
    full_match, partial_match, set_dups, missing_isoforms = sort_matches()
    classification_dict = create_class_dict()
    for iso in set_dups:
        if iso in classification_dict:
            del classification_dict[iso]
    return classification_dict

#writing new filtered classification filter
def write_new_classification_file():
    classification_final_dict = sort_class_dict()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for key in classification_final_dict:
            if key == "isoform":
                out.write(classification_final_dict[key])
            else:
                out.write(classification_final_dict[key])

write_new_classification_file()

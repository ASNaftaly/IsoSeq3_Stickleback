#Sorting BLAST results from BLAST_form9_pull_fullmatches.py
#will need to either get a single match for an isoform or decide to examine the isoform later
#should also add in chrom number of individual isoform and combined sexes isoform to further sort out which isoform matches to what
#to run script: python3 Sort_best_matches_BLAST.py <output from BLAST_form9_pull_fullmatches.py> <individual tissue classification file> <combined tissues classification file> <output sample name> <output file>
#Author: Alice Naftaly, February 2020

import sys


#reads best matches output
#creates dictionary with key = query isoform id and value == whole line
def read_best_matches():
    best_matches_input = sys.argv[1]
    best_matches_dict = {}
    with open(best_matches_input, 'r') as best_matches:
        for line in best_matches:
            new_line = line.split()
            isoform_id = new_line[0]
            if isoform_id in best_matches_dict:
                best_matches_dict[isoform_id].append(new_line)
            elif isoform_id not in best_matches_dict:
                best_matches_dict.update({isoform_id:[new_line]})
    return best_matches_dict

#read in individual tissues classification to get isoform chrom num and strand
def read_individual_class():
    individual_class_file = sys.argv[2]
    ind_class_dict = {}
    ind_length = {}
    with open(individual_class_file, 'r') as ind_class:
        for line in ind_class:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                length = int(new_line[3])
                dict_value = [chr_num, strand]
                ind_length.update({isoform:length})
                ind_class_dict.update({isoform:dict_value})
    return ind_class_dict, ind_length


#read in combined tissues classification to get isoform chrom num and strand
def read_combined_class():
    combined_class_file = sys.argv[3]
    combined_class_dict = {}
    combined_length = {}
    with open(combined_class_file, 'r') as com_class:
        for line in com_class:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                length = int(new_line[3])
                dict_value = [chr_num, strand]
                combined_length.update({isoform:length})
                combined_class_dict.update({isoform:dict_value})
    return combined_class_dict, combined_length


#sorting out best matches:
#include matches where query sequenced covers at least 60% of the subject sequence
#anything that doesn't match this criteria will be in excluded_matches
def sort_matches():
    best_matches = read_best_matches()
    individual_isoforms, individual_lengths = read_individual_class()
    combined_isoforms, combined_length = read_combined_class()
    kept_matches = {}
    excluded_matches = {}
    total = 0
    kept_count = 0
    removed_count = 0
    for isoform in best_matches:
        if isoform in individual_isoforms:
            #total += 1
            single_isoform = best_matches[isoform]
            #if there is only 1 isoform match
            if len(single_isoform) == 1:
                single = single_isoform[0]
                single_query_isoform = single[0]
                single_subject_isoform = single[1]
                q_length = int(single[3])
                q_start = int(single[6])
                q_end = int(single[7])
                s_start = int(single[8])
                s_end = int(single[9])
                marker = single[len(single)-1]
                comb_length = combined_length[single_subject_isoform]
                comb_length_60per = int(round(0.60 * comb_length, 0))
                ind_chr = individual_isoforms[single_query_isoform]
                comb_chr = combined_isoforms[single_subject_isoform]
                #total += 1
                #if both individual isoform and combined analysis isoform have the same chromosome and strand
                if ind_chr == comb_chr:
                    #if the alignment is an exact match
                    if q_start == s_start and q_end == s_end:
                        if isoform in kept_matches:
                            kept_matches[isoform].append(single)
                        elif isoform not in kept_matches:
                            kept_matches.update({isoform:[single]})
                        #kept_count += 1
                    #if alignment covers more than 60% of the subject isoform
                    elif q_length >= comb_length_60per:
                        if isoform in kept_matches:
                            kept_matches[isoform].append(single)
                        elif isoform not in kept_matches:
                            kept_matches.update({isoform:[single]})
                        #kept_count += 1
                    #if the alignment is less than 60% of the subject isoform, these will be removed for now
                    else:
                        if isoform in excluded_matches:
                            excluded_matches[isoform].append(single)
                        elif isoform not in excluded_matches:
                            excluded_matches.update({isoform:[single]})
                        #removed_count += 1
                #if the chromosome or strand isn't the same, then the matches are not correct
                else:
                    if isoform in excluded_matches:
                        excluded_matches[isoform].append(single)
                    elif isoform not in excluded_matches:
                        excluded_matches.update({isoform:[single]})
                    #removed_count += 1
            #if there is more than 1 match for a single isoform (matches can be full matches or where the length of alignment is greater than half of the length of the query sequence)
            elif len(single_isoform) > 1:
                for single in single_isoform:
                    single_query_isoform = single[0]
                    single_subject_isoform = single[1]
                    q_length = int(single[3])
                    q_start = int(single[6])
                    q_end = int(single[7])
                    s_start = int(single[8])
                    s_end = int(single[9])
                    marker = single[len(single)-1]
                    comb_length = combined_length[single_subject_isoform]
                    comb_length_60per = int(round(0.60 * comb_length, 0))
                    ind_chr = individual_isoforms[single_query_isoform]
                    comb_chr = combined_isoforms[single_subject_isoform]
                    #total += 1
                    #if both individual isoform and combined analysis isoform have the same chromosome and strand
                    if ind_chr == comb_chr:
                        #if the alignment is an exact match
                        if q_start == s_start and q_end == s_end:
                            if isoform in kept_matches:
                                kept_matches[isoform].append(single)
                            elif isoform not in kept_matches:
                                kept_matches.update({isoform:[single]})
                            #kept_count += 1
                        #if alignment covers more than 60% of the subject isoform
                        elif q_length >= comb_length_60per:
                            if isoform in kept_matches:
                                kept_matches[isoform].append(single)
                            elif isoform not in kept_matches:
                                kept_matches.update({isoform:[single]})
                            #kept_count += 1
                        #if the alignment is less than 60% of the subject isoform, these will be removed for now
                        else:
                            if isoform in excluded_matches:
                                excluded_matches[isoform].append(single)
                            elif isoform not in excluded_matches:
                                excluded_matches.update({isoform:[single]})
                            #print((q_length/comb_length)*100)
                            #removed_count += 1
                    #if the chromosome isn't the same, then the matches are not correct
                    else:
                        if isoform in excluded_matches:
                            excluded_matches[isoform].append(single)
                        elif isoform not in excluded_matches:
                            excluded_matches.update({isoform:[single]})
                        #removed_count += 1
    print("Number of Isoforms that met filter criteria 1")
    print(len(kept_matches))
    print("Number of Isoforms that were removed from filter criteria 1")
    print(len(excluded_matches))
    return kept_matches, excluded_matches



#need to go through the kept_matches from the sort function
#if an individual isoform matches more than 1 isoform from the combined data set, need to condense this to 1 isoform (preferably best match)
#   going to go with the isoform that matches best to the subject isoform (highest alignment length)
#if the matches are identical, will likely just remove isoform from analysis
def sort_kept_matches():
    matches, non_matches = sort_matches()
    individual_isoforms, individual_lengths = read_individual_class()
    combined_isoforms, combined_length = read_combined_class()
    single_matches = {}
    for isoform in matches:
        single_isoform = matches[isoform]
        if len(single_isoform) == 1:
            single = single_isoform[0]
            single_matches.update({isoform:single})
        else:
            q_length = individual_lengths[isoform]
            alignment_lengths = []
            for single in single_isoform:
                alignment_lengths.append(int(single[3]))
            max_alignment_length = max(alignment_lengths)
            alignment_index = alignment_lengths.index(max_alignment_length)
            single_matches.update({isoform:single_isoform[alignment_index]})
    return single_matches


#made excluded matches as dictionary so I can go back through it if I decide to (i.e. if more transcripts should be kept for criteria that are not present in the above functions)

#need to now sort out if single tissue matches to multiple combined sexes isoforms
#chooses the best alignment length for subject isoform so only one single tissue isoform matches only on combined tissues isoform
def collapse_matches():
    all_matches = sort_kept_matches()
    subject_dict = {}
    final_dict = {}
    for isoform in all_matches:
        single_isoform = all_matches[isoform]
        s_isoform = single_isoform[1]
        if s_isoform in subject_dict:
            subject_dict[s_isoform].append(single_isoform)
        elif s_isoform not in subject_dict:
            subject_dict.update({s_isoform:[single_isoform]})
    for key in subject_dict:
        single_key = subject_dict[key]
        if len(single_key) == 1:
            single = single_key[0]
            final_dict.update({key:single})
        else:
            alignment_lengths = []
            for single in single_key:
                alignment_lengths.append(int(single[3]))
            max_length = max(alignment_lengths)
            max_index = alignment_lengths.index(max_length)
            final_dict.update({key:single_key[max_index]})
    print("Number of Isoforms that matched only one combined sexes isoform")
    print(len(final_dict))
    return final_dict


#write BLAST matches to new output file
def write():
    final_matches = collapse_matches()
    sample_name = sys.argv[4]
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = sample_name + "\tCombined.Sexes.Isoform\n"
        out.write(header)
        for key in final_matches:
            single_key = final_matches[key]
            final = "%s\t%s\n" % (str(single_key[0]), str(single_key[1]))
            out.write(final)


write()

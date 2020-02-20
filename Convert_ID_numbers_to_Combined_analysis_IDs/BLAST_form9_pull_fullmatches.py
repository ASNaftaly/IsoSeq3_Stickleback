#Pulling BLAST matches (output format #9 - tab delimited file) for getting isoform id from combined sexes analysis with individual analyses
#BLAST parameters: Database = combined sexes nucleotides; query = individual tissue analysis
#will first pull out matches that cover most of the isoform sequence
#to run script: python3 BLAST_form9_pull_fullmatches.py <individual tissue classification file> <blast results> <output file>
#Author: Alice Naftaly, February 2020

import sys


#first need the length of individual tissue analysis isoforms
#will create a dictionary with key = isoform id and value = [chrom num, strand, length of isoform]
def pull_ind_isoform():
    ind_isoform_file = sys.argv[1]
    ind_isoform_dict = {}
    with open(ind_isoform_file, 'r') as isoforms:
        for line in isoforms:
            new_line = line.split()
            isoform_id = new_line[0]
            length = new_line[3]
            ind_isoform_dict.update({isoform_id:length})
    print("Number of Isoforms in Single Tissue")
    print(len(ind_isoform_dict))
    return ind_isoform_dict

#need to pull BLAST information
#used output format 9 = where there are 4 comment lines before each query where fourth line specifies labels for query lines
#order = query id, subject id (id in combined sexes), % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
#need to pull query id, subject id, alignment length, q. start, q. end, s. start, s.end
#will first create dictionary with key = query id and value = each line that matches to query id
def pull_BLAST_results():
    BLAST_output = sys.argv[2]
    BLAST_dict = {}
    with open(BLAST_output, 'r') as blast_results:
        for line in blast_results:
            if line.startswith("#"):
                continue
            else:
                new_line = line.split()
                query_id = new_line[0]
                if query_id in BLAST_dict:
                    BLAST_dict[query_id].append(new_line)
                elif query_id not in BLAST_dict:
                    BLAST_dict.update({query_id:[new_line]})
    print("Total BLAST matches")
    print(len(BLAST_dict))
    return BLAST_dict

#next, will sort through BLAST results and pull out those that are close to ind_query sequence length
#keeps all matches that match at least 1/2 of the individual query sequence
#returns dictionary with key = isoform and value == [blast output for each line] + "*" for full match and "-" for partial matches
def sort_BLAST_bestmatch():
    blast_dict = pull_BLAST_results()
    individual_isoforms = pull_ind_isoform()
    best_match_dict = {}
    for isoform in blast_dict:
        length_list = []
        full_matches = []
        partial_matches = []
        single_blast_results = blast_dict[isoform]
        single_individual_isoform = individual_isoforms[isoform]
        for index, single in enumerate(single_blast_results):
            length_of_query_match = int(single[3])
            value = [index, length_of_query_match]
            length_list.append(value)
        #sort length list based on position 1 in list (length of query match)
        sorted_length_list = sorted(length_list, key=lambda x: max(x[1:]), reverse=True)
        for match in sorted_length_list:
            match_index = match[0]
            match_length = match[1]
            if int(match_length) == int(single_individual_isoform):
                dict_value = single_blast_results[match_index]
                dict_value += "*"
                full_matches.append(dict_value)
            #this pulls matches that match at least have the length of the query sequence
            elif int(match_length) > int(int(single_individual_isoform)/2):
                dict_value = single_blast_results[match_index]
                dict_value += "-"
                partial_matches.append(dict_value)
        for fm in full_matches:
            if isoform in best_match_dict:
                best_match_dict[isoform].append(fm)
            elif isoform not in best_match_dict:
                best_match_dict.update({isoform:[fm]})
        for pm in partial_matches:
            if isoform in best_match_dict:
                best_match_dict[isoform].append(pm)
            elif isoform not in best_match_dict:
                best_match_dict.update({isoform:[pm]})
    print("Number of Isoforms that matched greater than 50% of Subject")
    print(len(best_match_dict))
    return best_match_dict

#write best matches to an output file
#writes blast output to a file to be further processed to only 1 match per isoform
def write_best_matches():
    best_match_dict = sort_BLAST_bestmatch()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for isoform in best_match_dict:
            single_isoform = best_match_dict[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                final_value = "\t".join(single)
                out.write(final_value + "\n")
            else:
                for single in single_isoform:
                    final_value = "\t".join(single)
                    out.write(final_value + "\n")


write_best_matches()

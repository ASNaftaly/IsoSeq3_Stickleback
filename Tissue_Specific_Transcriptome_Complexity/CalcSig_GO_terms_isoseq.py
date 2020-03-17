#Determining if GO terms are enriched in different sets of tissue comparisons or tissue specific genes
#1. each random distribution file has the format: iteration# \t GO.term \t counts of that go term in that iteration
#2. Make dictionary of random GO terms with key = GO term and dictionary value = number of occurrences for each GO term and iteration of distribution
#3. add "0" to each key until the total length is 10,000 values
#4. read in observed go term occurrences for specfic observed file and save as dictionary with key = GO term and value = number of observed occurrences
#5. Go key by key and ask how many times the observed value is greater then the random value; get total and use 1- (total/10000) to get pvalue
#6. Write final file: GO term \t pvalue
#to run script: python3 CalcSig_GO_terms_isoseq.py <permutations 10k file> <observed GO terms file> <output file; format = GO.Term \t P>
#Author: Alice Naftaly, March 2020

import sys

#pull permutation counts
#each line = iteration \t GO term \t counts of GO term in that iteration
#creates dictionary with key == GO term and value == list of counts
def pull_permutations():
    permutation_file = sys.argv[1]
    permutation_dict = {}
    with open(permutation_file, "r") as permutations:
        for line in permutations:
            new_line = line.split()
            go_term = new_line[1]
            counts = new_line[2]
            if go_term in permutation_dict:
                permutation_dict[go_term].append(int(counts))
            elif go_term not in permutation_dict:
                permutation_dict.update({go_term:[int(counts)]})
    return permutation_dict

#not all GO terms were seen in every permutation (in fact most are not)
#need to add "0" values get the total number of counts per go term to 10,000
#returns dictionary with key == GO term and value == list of counts (total 10k per go term)
def add_counts():
    permutations = pull_permutations()
    for go_term in permutations:
        num_counts = len(permutations[go_term])
        while num_counts < 10000:
            permutations[go_term].append(0)
            num_counts += 1
    return permutations


#pulls observed GO terms counts
#return dictionary with key == go term and value = number of observed occurrances of that go term in actual set of genes/GO terms seen
def pull_observed_terms():
    observed_file = sys.argv[2]
    observed_dict = {}
    final_dict = {}
    with open(observed_file, 'r') as observed:
        for line in observed:
            new_line = line.split()
            go_terms = new_line[1].split(",")
            for go in go_terms:
                if go in observed_dict:
                    observed_dict[go].append("1")
                elif go not in observed_dict:
                    observed_dict.update({go:["1"]})
    for key in observed_dict:
        num_counts = len(observed_dict[key])
        final_dict.update({key:num_counts})
    return final_dict

#calculating P value for over enriched GO terms
#will definitely need to correct for multiple testing
#returns dictionary with key == go term and value = p value (overenriched)
def calc_P():
    observed_go_terms = pull_observed_terms()
    permutations = add_counts()
    p_values_dict = {}
    for go_term in observed_go_terms:
        if go_term in permutations:
            single_observed = observed_go_terms[go_term]
            single_permutation = permutations[go_term]
            single_overenriched = 0
            single_underenriched = 0
            for single_count in single_permutation:
                #if observed value is over enriched; it should be greater than the counts seen in the permutations
                if single_count <= int(single_observed):
                    single_overenriched += 1
            #over-enriched p value
            p_overenriched = round(1 - (single_overenriched/10000),5)
            p_values_dict.update({go_term:p_overenriched})
    return p_values_dict

#write p values to file with GO terms:
#final output = GO.Term \t P
def write():
    p_values = calc_P()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "GO.Term\tP\n"
        for go_term in p_values:
            single_p_value = p_values[go_term]
            final = "%s\t%s\n" % (str(go_term), str(single_p_value))
            out.write(final)

write()

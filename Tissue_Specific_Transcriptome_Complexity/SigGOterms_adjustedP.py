#Script to pull significant GO terms after adjusting P value
#Need to pull GO terms with P values from pvalues file, read in adjusted p value,
#then need to write new output file with GO terms that are significant based on the adjusted P value
#to run script: SigGOterms_adjustedP.py <p values file with all observed GO terms and unadjusted p values> <adjusted p value in decimal form> <output file, format: GO.Term\tP>
#Author: Alice Naftaly, March 2020

import sys

#read in p values file
#returns dictionary with key == go term and value == unadjusted p value
def read_significance_file():
    sig_file = sys.argv[1]
    sig_dict = {}
    with open(sig_file, 'r') as sig:
        for line in sig:
            new_line = line.split()
            go_term = new_line[0]
            p_value = new_line[1]
            sig_dict.update({go_term:p_value})
    return sig_dict


#pull significant p values and GO terms based on bonferroni adjusted pvalues
#returns p values that are significant given the adjusted p value in a dictionary
#dictionary has key == go term and value == p value
def adjusted_p_values():
    sig_dict = read_significance_file()
    adjusted_p = float(sys.argv[2])
    adjusted_dict = {}
    for go_term in sig_dict:
        single_p = sig_dict[go_term]
        if float(single_p) <= adjusted_p:
            adjusted_dict.update({go_term:single_p})
    return adjusted_dict


#write significant go terms to new file
def write():
    adjusted_dict = adjusted_p_values()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header="GO.Term\tP\n"
        out.write(header)
        for go_term in adjusted_dict:
            single_p = adjusted_dict[go_term]
            final = "%s\t%s\n" % (str(go_term), str(single_p))
            out.write(final)

write()

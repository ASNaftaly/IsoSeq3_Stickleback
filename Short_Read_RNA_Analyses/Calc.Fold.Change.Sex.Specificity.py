#calculating fold change in TPM between female/male isoforms compared to female and male sequences (using Kallisto)
#will read in kallisto files for female to female or male to male to compare to female to male or male to female
#to run script: python3 Calc.Fold.Change.Sex.Specificity.py <kallisto tsv for female/female or male/male file> <kallisto tsv for female to male or male to female file> <output file>
#Author: Alice Naftaly, June 2020

import sys

#read kallisto tsv for isoforms to same fastqs
#returns dictionary with key == isoform and value == [transcript length, tpm value]
def read_same_tsv():
    tsv_file = sys.argv[1]
    same_tsv_dict = {}
    with open(tsv_file, 'r') as tsv:
        for line in tsv:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                transcript_length = int(new_line[1])
                tpm = float(new_line[4])
                dict_value = [transcript_length, tpm]
                same_tsv_dict.update({isoform:dict_value})
    return same_tsv_dict


#read kallisto tsv for isoforms to opposite fastqs
#returns dictionary with key == isoform and value == [transcript length, tpm value]
def read_diff_tsv():
    tsv_file = sys.argv[2]
    diff_tsv_dict = {}
    with open(tsv_file, 'r') as tsv:
        for line in tsv:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                transcript_length = int(new_line[1])
                tpm = float(new_line[4])
                dict_value = [transcript_length, tpm]
                diff_tsv_dict.update({isoform:dict_value})
    return diff_tsv_dict


#calculate fold Change
def calc_fold_change():
    same_tsv = read_same_tsv()
    diff_tsv = read_diff_tsv()
    fold_change_dict = {}
    for isoform in same_tsv:
        single_same_isoform = same_tsv[isoform]
        single_diff_isoform = diff_tsv[isoform]
        single_same_transcript_length = single_same_isoform[0]
        single_same_tpm = single_same_isoform[1]
        single_diff_tpm = single_diff_isoform[1]
        if single_same_tpm == 0 or single_diff_tpm == 0:
            fold_change == "0"
        #positive fold change = diff is greater than same
        elif single_same_tpm < single_diff_tpm:
            fold_change = single_diff_tpm/single_same_tpm
        #no Change
        elif single_same_tpm == single_diff_tpm:
            fold_change = single_diff_tpm/single_same_tpm
        #negative fold change = same is greater than diff
        elif single_same_tpm > single_diff_tpm:
            fold_change = single_same_tpm/single_diff_tpm * - 1
        dict_value = [single_same_transcript_length, fold_change]
        fold_change_dict.update({isoform:dict_value})
    return fold_change_dict

#write fold Change
def write():
    fold_change = calc_fold_change()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Isoform.ID\tTranscript.Length\tFold.Change\n"
        for isoform in fold_change:
            single = fold_change[isoform]
            final = "%s\t%s\t%s\n" % (str(isoform), str(single[0]), str(single[1]))
            out.write(final)


write()

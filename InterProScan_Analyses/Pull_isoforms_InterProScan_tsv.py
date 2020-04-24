#pulling specific isoforms from InterProScan tsv output file
#to run script: python3 Pull_isoforms_InterProScan_tsv.py <full tsv file> <novel gene or novel transcript classification file> <output filtered tsv file>
#Author: Alice Naftaly, April 2020

import sys


#read tsv file
#returns dictionary with key == isoform id and value == line for each isoform for tsv file
def read_tsv():
    tsv_file = sys.argv[1]
    tsv_dict = {}
    with open(tsv_file, 'r') as tsv:
        for line in tsv:
            new_line = line.split("\t")
            isoform = new_line[0]
            if isoform in tsv_dict:
                tsv_dict[isoform].append(line)
            elif isoform not in tsv_dict:
                tsv_dict.update({isoform:[line]})
    return tsv_dict


#read in classification file with novel transcripts or novel genes
#returns list of isoform ids
def read_class():
    class_file = sys.argv[2]
    isoform_list = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                isoform_list.append(isoform_id)
    return isoform_list


#pull tsv lines for specific isoforms:
#returns list of lines to go into final output
def filter_tsv():
    tsv_dict = read_tsv()
    novel_isoforms = read_class()
    final_lines = []
    for isoform in novel_isoforms:
        if isoform in tsv_dict:
            single_tsv = tsv_dict[isoform]
            for value in single_tsv:
                final_lines.append(value)
    return final_lines


#write output
def write():
    final_lines = filter_tsv()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for line in final_lines:
            out.write(line)

write()

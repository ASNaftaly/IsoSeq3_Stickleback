#sort classification file based on final BLAST comparison
#read in classification file and exon counts file
#to run script: python3 Sort_Class_file_by_exon_counts.py <classification file> <exon counts file> <output filtered class file>
#Author: Alice Naftaly, April 2020

import sys

#read classification file
#returns dictionary with key == isoform id and value == line
def read_class_file():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                class_dict.update({isoform:new_line})
    return class_dict


#read exon counts file
#return list of novel_isoforms
def read_exon_counts():
    exon_file = sys.argv[2]
    all_isoforms = []
    with open(exon_file, 'r') as exon_counts:
        for line in exon_counts:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                all_isoforms.append(isoform)
    return all_isoforms


#filter class dict:
#returns filtered class_dict
def filter_class():
    class_dict = read_class_file()
    isoforms = read_exon_counts()
    filtered_dict = {}
    for isoform in isoforms:
        if isoform in class_dict:
            filtered_dict.update({isoform:class_dict[isoform]})
    return filtered_dict

#write filtered class dict
def write():
    filtered_dict = filter_class()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for isoform in filtered_dict:
            single_isoform = filtered_dict[isoform]
            final = "\t".join(single_isoform)
            out.write(final + "\n")

write()

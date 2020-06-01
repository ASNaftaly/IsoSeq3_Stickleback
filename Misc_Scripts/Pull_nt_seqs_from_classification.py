#pull fasta sequences based on the classification file
#similar to Pull_seqs_exons_filtered.py except only pull nucleotide sequences
#to run script: python3 Pull_nt_seqs_from_classification.py <classification file> <fasta file> <output new fasta file with only ids in classification file>
#Author: Alice Naftaly, May 2020

import sys


#read in classification file
#just need the isoform ids
#returns list of isoform ids
def pull_isoform_ids():
    class_file = sys.argv[1]
    isoforms = []
    with open(class_file, 'r') as class_input:
        for line in class_input:
            if line.startswith("PB"):
                new_line = line.split()
                isoforms.append(new_line[0])
    return isoforms


#read in fasta file
#returns fasta dictionary with key = isoform id and value == header line for sequence, followed by sequence with \n kept
def read_fasta():
    fasta_file = sys.argv[2]
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                isoform_id = new_line[0].strip(">")
                fasta_dict.update({isoform_id:[line]})
            else:
                fasta_dict[isoform_id].append(line)
    return fasta_dict


#filter out ids not in classification file
def filter_fasta():
    isoforms = pull_isoform_ids()
    fasta_dict = read_fasta()
    isos_to_remove = []
    for iso in fasta_dict:
        if iso not in isoforms:
            isos_to_remove.append(iso)
    for i in isos_to_remove:
        if i in fasta_dict:
            del fasta_dict[i]
    return fasta_dict


#writing duplicate removed fasta, faa, and gtf files

def write_fasta():
    filtered_fasta = filter_fasta()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for iso in filtered_fasta:
            single_iso = filtered_fasta[iso]
            for single in single_iso:
                out.write(single)

write_fasta()

#pulling all non coding sequences from combined sexes files with duplicates removed
#need to read in faa headers, fasta file, gtf file, and classification file
#to run script: Pull_noncoding_isoforms.py <combined sexes classification file> <faa file> <fasta file> <gtf file> <outout noncoding isoforms classification file> <output noncoding isoforms fasta file> <output noncoding isoforms gtf file>
#Author: Alice Naftaly, April 2020

import sys

#read in classification file:
#returns dictionary with key == isoform and value == classification line
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("isoform"):
                class_dict.update({"header":line})
            elif line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                class_dict.update({isoform:line})
    return class_dict


#read in faa file
#just need the isoform ids from the headers
#returns isoform list
def read_faa():
    faa_file = sys.argv[2]
    isoform_list = []
    with open(faa_file, 'r') as amino_acids:
        for line in amino_acids:
            if line.startswith(">"):
                new_line = line.split('\t')
                isoform = new_line[0].strip(">")
                isoform_list.append(isoform)
    return isoform_list


#read in fasta file
#first dictionary gets all of the sequence for each isoform together
#final dictionary returned has key == isoform and value = nucleotide sequence
def read_fasta():
    fasta_file = sys.argv[3]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as nucleotides:
        for line in nucleotides:
            if line.startswith(">"):
                new_line = line.split()
                isoform_id = new_line[0].strip(">")
            else:
                if isoform_id in fasta_dict:
                    fasta_dict[isoform_id].append(line.strip("\n"))
                elif isoform_id not in fasta_dict:
                    fasta_dict.update({isoform_id:[line.strip("\n")]})
    for key in fasta_dict:
        single_key = fasta_dict[key]
        final_seq = []
        for seq in single_key:
            list_seq = list(seq)
            final_seq += list_seq
        final_fasta_dict.update({key:final_seq})
    return final_fasta_dict


#read in gtf file
#returns dictionary with key == isoform and value == every exon in the isoform
def read_gtf():
    gtf_file = sys.argv[4]
    gtf_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            isoform_info = new_line[8].split(";")
            isoform_full = isoform_info[0].split(" ")[1]
            final_isoform = isoform_full.strip("\"")
            if final_isoform in gtf_dict:
                gtf_dict[final_isoform].append(line)
            elif final_isoform not in gtf_dict:
                gtf_dict.update({final_isoform:[line]})
    return gtf_dict


#Now to filter for only non protein coding isoforms
#will do this by comparing the class_dict and faa_dict where if the isoform is not in the faa_dict, then it was predicted to be noncoding
#returns list of noncoding isoforms
def pull_noncoding():
    class_dict = read_class()
    faa_list = read_faa()
    noncoding_isoforms = []
    for isoform in class_dict:
        if isoform not in faa_list and isoform.startswith("PB"):
            noncoding_isoforms.append(isoform)
    return noncoding_isoforms


#now need to write noncoding isoforms classification, fasta, and gtf files
def write_noncoding_class():
    class_dict = read_class()
    noncoding_isoforms = pull_noncoding()
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = class_dict["header"]
        out.write(header)
        for isoform in noncoding_isoforms:
            single_class_entry = class_dict[isoform]
            out.write(single_class_entry)


def write_noncoding_fasta():
    fasta_dict = read_fasta()
    noncoding_isoforms = pull_noncoding()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for isoform in noncoding_isoforms:
            single_fasta = fasta_dict[isoform]
            sequence = "".join(single_fasta)
            final_header = ">%s\n" % str(isoform)
            final_seq = sequence + "\n"
            out.write(final_header)
            out.write(final_seq)


def write_noncoding_gtf():
    gtf_dict = read_gtf()
    noncoding_isoforms = pull_noncoding()
    output = sys.argv[7]
    with open(output, 'a') as out:
        for isoform in noncoding_isoforms:
            single_gtf = gtf_dict[isoform]
            for exon in single_gtf:
                out.write(exon)


#call all functions
def call():
    class_out = write_noncoding_class()
    fasta_out = write_noncoding_fasta()
    gtf_out = write_noncoding_gtf()

call()

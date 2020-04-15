#Sort non coding isoforms based on size
#anything equal to or under 200bp long will be considered a short ncRNA
#anything over 200 bp long will be considered a long ncRNA
#will go ahead and split classification, fasta and gtf for short vs long ncRNAs (will likely want to upload these when we submit the paper)
#to run script: Sort_noncoding_RNAs_by_size.py <non coding isoforms classification file> <non coding isoforms fasta file> <non coding isoforms gtf file> <output short noncoding isoforms classification> <output long noncoding isoforms classification> <output short noncoding isoforms fasta> <output long noncoding isoforms fasta> <output short noncoding isoforms gtf> <output long noncoding isoforms gtf>
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

#read in fasta file
#first dictionary gets all of the sequence for each isoform together
#final dictionary returned has key == isoform and value = nucleotide sequence
def read_fasta():
    fasta_file = sys.argv[2]
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
    gtf_file = sys.argv[3]
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


#separating short vs long ncRNAs
#returns dictionary where key = "short" or "long" and value == list of isoforms that fall under one category
def sort_by_size():
    nc_seqs = read_fasta()
    size_dict = {}
    for isoform in nc_seqs:
        single_isoform_seq = nc_seqs[isoform]
        length_seq = len(single_isoform_seq)
        if length_seq <= 200:
            if "short" in size_dict:
                size_dict["short"].append(isoform)
            elif "short" not in size_dict:
                size_dict.update({"short":[isoform]})
        elif length_seq > 200:
            if "long" in size_dict:
                size_dict["long"].append(isoform)
            elif "long" not in size_dict:
                size_dict.update({"long":[isoform]})
    return size_dict


#write new split classification, fasta, and gtf files:
def write_short_class():
    class_dict = read_class()
    size_dict = sort_by_size()
    short_isoforms = size_dict["short"]
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = class_dict["header"]
        out.write(header)
        for isoform in short_isoforms:
            single_class_entry = class_dict[isoform]
            out.write(single_class_entry)

def write_long_class():
    class_dict = read_class()
    size_dict = sort_by_size()
    long_isoforms = size_dict["long"]
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = class_dict["header"]
        out.write(header)
        for isoform in long_isoforms:
            single_class_entry = class_dict[isoform]
            out.write(single_class_entry)

def write_short_fasta():
    fasta_dict = read_fasta()
    size_dict = sort_by_size()
    short_isoforms = size_dict["short"]
    output = sys.argv[6]
    with open(output, 'a') as out:
        for isoform in short_isoforms:
            single_fasta = fasta_dict[isoform]
            sequence = "".join(single_fasta)
            final_header = ">%s\n" % str(isoform)
            final_seq = sequence + "\n"
            out.write(final_header)
            out.write(final_seq)

def write_long_fasta():
    fasta_dict = read_fasta()
    size_dict = sort_by_size()
    long_isoforms = size_dict["long"]
    output = sys.argv[7]
    with open(output, 'a') as out:
        for isoform in long_isoforms:
            single_fasta = fasta_dict[isoform]
            sequence = "".join(single_fasta)
            final_header = ">%s\n" % str(isoform)
            final_seq = sequence + "\n"
            out.write(final_header)
            out.write(final_seq)

def write_short_gtf():
    gtf_dict = read_gtf()
    size_dict = sort_by_size()
    short_isoforms = size_dict["short"]
    output = sys.argv[8]
    with open(output, 'a') as out:
        for isoform in short_isoforms:
            single_gtf = gtf_dict[isoform]
            for exon in single_gtf:
                out.write(exon)

def write_long_gtf():
    gtf_dict = read_gtf()
    size_dict = sort_by_size()
    long_isoforms = size_dict["long"]
    output = sys.argv[9]
    with open(output, 'a') as out:
        for isoform in long_isoforms:
            single_gtf = gtf_dict[isoform]
            for exon in single_gtf:
                out.write(exon)


#call all functions:
def call():
    short_class = write_short_class()
    long_class = write_long_class()
    short_fasta = write_short_fasta()
    long_fasta = write_long_fasta()
    short_gtf = write_short_gtf()
    long_gtf = write_long_gtf()

call()

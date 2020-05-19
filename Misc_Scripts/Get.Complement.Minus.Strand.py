#converting sequences on - strand to complementary strand
#in transcriptome (.fasta file) all sequences are in the 5' to 3' direction
#simply need to get complement of isoforms on the - strand
#will need to read in classification file to get strand information
#to run script: python3 Get.Complement.Minus.Strand.py <classification file> <fasta file> <output fasta file>
#Author: Alice Naftaly, May 2020

import sys

#read classification file and get isoform id with strand information
#returns dictionary with key == strand (+ or -) and value == list of isoform ids
def read_class():
    class_file = sys.argv[1]
    strand_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                strand = new_line[2]
                if strand in strand_dict:
                    strand_dict[strand].append(isoform)
                elif strand not in strand_dict:
                    strand_dict.update({strand:[isoform]})
    return strand_dict


#read fasta file
#returns dictionary with key == isoform id and value == sequence
def read_fasta():
    fasta_file = sys.argv[2]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                isoform_id = new_line[0].strip(">")
                fasta_dict.update({isoform_id:[line]})
            else:
                fasta_dict[isoform_id].append(line)
    for isoform in fasta_dict:
        single_isoform = fasta_dict[isoform]
        full_seq = []
        for value in single_isoform:
            if value.startswith(">"):
                final_fasta_dict.update({isoform:[value]})
            else:
                split_value = value.strip("\n")
                full_seq += list(split_value)
        final_fasta_dict[isoform].append(full_seq)
    return final_fasta_dict


#filter fasta dict into 2 dictionaries for + and - strands
#returns 2 dictionaries with one having all + strand isoforms and the other having all - strand isoforms
#each dictionary is formatted as: key == isoform and value == sequence
def filter_fasta():
    strand_dict = read_class()
    sequences = read_fasta()
    plus_strand_sequences = {}
    minus_strand_sequences = {}
    for strand in strand_dict:
        single_strand = strand_dict[strand]
        for isoform in single_strand:
            if strand == "+":
                plus_strand_sequences.update({isoform:sequences[isoform]})
            elif strand == "-":
                minus_strand_sequences.update({isoform:sequences[isoform]})
    return plus_strand_sequences, minus_strand_sequences



#get complement for minus strand
def get_complement():
    plus_strand_sequences, minus_strand_sequences = filter_fasta()
    complement_minus_strand = {}
    for isoform in minus_strand_sequences:
        single_seq = minus_strand_sequences[isoform]
        header = single_seq[0]
        complement_minus_strand.update({isoform:[header]})
        sequence = single_seq[1]
        new_seq = []
        for value in sequence:
            for nt in value:
                if nt == "A" or nt == "a":
                    new_nt = "T"
                elif nt == "C" or nt == "c":
                    new_nt = "G"
                elif nt == "G" or nt == "g":
                    new_nt = "C"
                elif nt == "T" or nt == "t":
                    new_nt = "A"
                new_seq += new_nt
        new_seq.reverse()
        final_seq = ''.join(new_seq)
        final_sequence = final_seq + "\n"
        complement_minus_strand[isoform].append(final_sequence)
    return complement_minus_strand


#write new fasta file
def write():
    plus_strand_sequences, old_minus_strand_sequences = filter_fasta()
    correct_minus_strand_sequences = get_complement()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for plus_isoform in plus_strand_sequences:
            single_sequence = plus_strand_sequences[plus_isoform]
            header = single_sequence[0]
            sequence = single_sequence[1]
            joined_seq = "".join(sequence)
            out.write(header)
            final = joined_seq + "\n"
            out.write(final)
        for minus_isoform in correct_minus_strand_sequences:
            single_sequence = correct_minus_strand_sequences[minus_isoform]
            header = single_sequence[0]
            sequence = single_sequence[1]
            out.write(header)
            out.write(sequence)


write()

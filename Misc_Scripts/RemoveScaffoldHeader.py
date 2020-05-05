#remove parts of header lines for scaffolds in fasta file
#to run script: python3 RemoveScaffoldHeader.py <ensembl unmasked.fa > <output file>
#output file is in fasta format
#Author: Alice Naftaly, May 2020


import sys


#read in fasta file:
#returns dictionary with key = header (>) and value == sequence
def read_fasta():
    fasta_file = sys.argv[1]
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                header_line = line.strip("\n")
            else:
                if header_line in fasta_dict:
                    fasta_dict[header_line].append(line)
                elif header_line not in fasta_dict:
                    fasta_dict.update({header_line:[line]})
    return fasta_dict


#cut header lines for scaffolds
#all saffolds have headers like: >scaffold_XXXX dna:scaffold scaffold:BROADS1:scaffold_XXX:1:XX:1 REF
#want to remove everything after >scaffold_XXXX
#returns dictionary with key == position (group or scaffold) and value == sequence
def format_headers():
    fasta_dict = read_fasta()
    final_dict = {}
    for position in fasta_dict:
        if  position.startswith(">group"):
            single_position = fasta_dict[position]
            final_dict.update({position:single_position})
        elif position.startswith(">scaffold"):
            single_position = fasta_dict[position]
            header = position.split(" ")
            final_position = header[0]
            final_dict.update({final_position:single_position})
    return final_dict

#write final fasta
def write():
    final_dict = format_headers()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for header in final_dict:
            print(header)
            single_position = final_dict[header]
            out.write(header + "\n")
            for value in single_position:
                out.write(value)

write()

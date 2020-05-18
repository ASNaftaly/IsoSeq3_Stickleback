#calculate normalization factor for PSI calculations
#normalization factor = # of maximum mappable positions / # of unique mappable positions
# # of maximum mappable positions = 135 (k-15)
# # of unique mappable positions = in file Unique.Kmers.Alignment.for.Normalization.txt
#need to read in unique kmers alignments and calculate normalization factor for each exon junction pair
#output format: Junction.Identifier \t Normalization.Factor
#to run script: python3 Calc.PSI.Normalization.Factor.py <135> <Unique,Kmers.Alignment.for.Normalization.txt> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read in unique mappable kmers
#returns dictionary with key == junction identifier and value == number of unique mappable kmers
def read_unique_positions():
    kmers_file = sys.argv[2]
    junction_dict = {}
    with open(kmers_file, 'r') as kmers:
        for line in kmers:
            new_line = line.split()
            junction_iden = new_line[0]
            num_unique_positions = int(new_line[1])
            junction_dict.update({junction_iden:num_unique_positions})
    return junction_dict


#calculate normalization factors
#returns dictionary with key = junction identifier and value == normalization factor rounded to 2 decimal points
def calc_norm_factor():
    junctions = read_unique_positions()
    max_mappable_positions = int(sys.argv[1])
    norm_factor_dict = {}
    for junction in junctions:
        unique_pos = junctions[junction]
        if unique_pos == 0:
            continue
        else:
            norm_factor = round(max_mappable_positions/unique_pos,2)
            norm_factor_dict.update({junction:norm_factor})
    return norm_factor_dict


#write normalization factors to file
def write():
    norm_factors = calc_norm_factor()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for junction in norm_factors:
            final = "%s\t%s\n" % (str(junction),str(norm_factors[junction]))
            out.write(final)

write()

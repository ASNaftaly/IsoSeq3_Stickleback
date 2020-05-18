#Counting number of unique kmers that aligned to any given junction pair
#reads in bed file created from mapq 30 filtered and sorted bam file and outputs a file with the following format:
#Junction.Name \t Unique.kmer.Alignments
#to run script: python3 Count.Kmers.Alignments.py <bed file> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read in bed file
#returns dictionary with key == junction identifier (chr num + strand with gene id and exon number as well as no_as1 (exon 1-2 junction), no_as2 (exon 2-3 junction), or as1 (exon 1-3 junction))
def read_bed():
    bed_file = sys.argv[1]
    bed_dict = {}
    with open(bed_file, 'r') as bed_info:
        for line in bed_info:
            new_line = line.split()
            junction_iden = new_line[0]
            start_pos = int(new_line[1])
            end_pos = int(new_line[2])
            kmer_iden = new_line[3]
            dict_value = [start_pos, end_pos, kmer_iden]
            if junction_iden in bed_dict:
                bed_dict[junction_iden].append(dict_value)
            elif junction_iden not in bed_dict:
                bed_dict.update({junction_iden:[dict_value]})
    return bed_dict


#filter out non unique Alignments
#if two or more kmers align to the same junction, all will be thrown out
def filter_unique_kmers():
    bed_dict = read_bed()
    filtered_dict = {}
    for junction in bed_dict:
        single_junction = bed_dict[junction]
        kmers_start = []
        kmers_end = []
        kmers_iden = []
        for kmer in single_junction:
            kmers_start.append(kmer[0])
            kmers_end.append(kmer[1])
            kmers_iden.append(kmer[2])
        #identify which start values are duplicates
        unique_kmers = set()
        duplicate_kmers = []
        for idx, val in enumerate(kmers_start):
            if val not in unique_kmers:
                unique_kmers.add(val)
            else:
                duplicate_kmers.append(val)
        overlapping_kmer_alignments = []
        final_unique_kmers = []
        #get indeces for duplicate values
        for overlap in duplicate_kmers:
            overlapping_kmer_alignments += [i for i, x in enumerate(kmers_start) if x == overlap]
        #remove indeces from kmer_idens
        for index, k in enumerate(kmers_iden):
            if index in overlapping_kmer_alignments:
                continue
            elif index not in overlapping_kmer_alignments:
                final_unique_kmers.append(k)
        filtered_dict.update({junction:len(final_unique_kmers)})
    return filtered_dict


#write unique mappable positions for each junction
def write():
    unique_kmers = filter_unique_kmers()
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "Junction.Name\tUnique.Kmer.Alignment\n"
        for junction in unique_kmers:
            final = "%s\t%s\n" % (str(junction), str(unique_kmers[junction]))
            out.write(final)

write()

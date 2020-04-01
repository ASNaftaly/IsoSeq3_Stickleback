#Counting the number of transcripts covered by CCS reads to test for sequencing saturation
#will read in blast output (format 9), fasta file from downsampled bam file for each sample, clasification file with transcripts for each sample
#to run script: python3 Counting_Transcripts_Seq_Saturation.py <classification file for sample> <fasta file from downsampled bam file> <blast output> <output file, format = Isoform.ID \t Number.of.CCS.Reads>
#Author: Alice Naftaly, April 2020

import sys
import re

#read classification file
#returns dictionary with key == isoform id and value == [gene id, isoform length]
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                gene_id = new_line[6]
                length = new_line[3]
                final = [gene_id, length]
                class_dict.update({isoform_id:final})
    print("Number of Total Isoforms")
    print(len(class_dict))
    return class_dict


#read fasta file and pull length of each ccs read
#creates first dictionary where key = read identifier and value == list of sequence
#returns final dictionary where key = read identifier and value == length of sequence
def read_ccs_fasta():
    fasta_file = sys.argv[2]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.strip("\n")
                final_identifier = new_line.strip(">")
            else:
                new_line = line.strip('\n')
                list_line = list(new_line)
                fasta_dict.update({final_identifier:list_line})
    for read in fasta_dict:
        single_read = fasta_dict[read]
        read_length = len(single_read)
        final_fasta_dict.update({read:read_length})
    print("Number of ccs reads")
    print(len(final_fasta_dict))
    return final_fasta_dict


#read in blast output
#returns dictionary with key == query id (ccs read identifier) and value = [subject_id (isoform id), alignment length, q_start, q_end, s_start, s_end]
def read_blast():
    blast_file = sys.argv[3]
    blast_dict = {}
    with open(blast_file, 'r') as blast:
        for line in blast:
            if line.startswith("#"):
                continue
            else:
                new_line = line.split("\t")
                query_id = new_line[0]
                subject_id = new_line[1]
                alignment_length = new_line[3]
                q_start = new_line[6]
                q_end = new_line[7]
                s_start = new_line[8]
                s_end = new_line[9]
                final = [subject_id, alignment_length, q_start, q_end, s_start, s_end]
                if query_id in blast_dict:
                    blast_dict[query_id].append(final)
                elif query_id not in blast_dict:
                    blast_dict.update({query_id:[final]})
    return blast_dict


#compare alignment length and length of ccs read
#pulls ccs reads that have blast hits of greater than 50%
#returns dictionary where key = ccs identifier and value == [isoform id, alignment_length, q_start, q_end, s_start, s_end]
#there will be entries that have more than one alignment
def compare_length_ccs_read():
    blast_dict = read_blast()
    ccs_reads = read_ccs_fasta()
    filtered_blast_dict = {}
    total_reads = 0
    for identifier in ccs_reads:
        single_ccs_read = int(ccs_reads[identifier])
        partial_length_50per = int(round(single_ccs_read * 0.50,0))
        if identifier in blast_dict:
            single_blast_dict = blast_dict[identifier]
            for single in single_blast_dict:
                alignment_length = int(single[1])
                if alignment_length >= partial_length_50per:
                    total_reads += 1
                    if identifier in filtered_blast_dict:
                        filtered_blast_dict[identifier].append(single)
                    elif identifier not in filtered_blast_dict:
                        filtered_blast_dict.update({identifier:[single]})
    print("Number of CCS read that pass 50% alignment filter")
    print(total_reads)
    return filtered_blast_dict


#compare alignment length to length of isoform
#pulls isoforms that have alignments of greater than 50%
#returns dictionary where key == isoform id and value == list of "1"s for each time a CCS read matches the requirements of read being >=50% of alignment and alignment being >= 50% of the isoform length
#some isoforms will have more than one ccs read that matches them
def compare_length_isoform():
    blast_dict = compare_length_ccs_read()
    isoform_dict = read_class()
    final_isoforms = {}
    total_counts = 0
    for ccs_identifier in blast_dict:
        single_blast = blast_dict[ccs_identifier]
        if len(single_blast) == 1:
            single = single_blast[0]
            isoform_id = single[0]
            alignment_length = int(single[1])
            if isoform_id in isoform_dict:
                single_isoform = isoform_dict[isoform_id]
                isoform_length = int(single_isoform[1])
                #50% of isoform length
                isoform_50per = int(round(isoform_length * 0.50,0))
                if alignment_length >= isoform_50per:
                    total_counts += 1
                    if isoform_id in final_isoforms:
                        final_isoforms[isoform_id].append("1")
                    elif isoform_id not in final_isoforms:
                        final_isoforms.update({isoform_id:["1"]})
        elif len(single_blast) > 1:
            for single in single_blast:
                isoform_id = single[0]
                alignment_length = int(single[1])
                if isoform_id in isoform_dict:
                    single_isoform = isoform_dict[isoform_id]
                    #50% of isoform length
                    isoform_50per = int(round(isoform_length * 0.50, 0))
                    if alignment_length >= isoform_50per:
                        total_counts += 1
                        if isoform_id in final_isoforms:
                            final_isoforms[isoform_id].append("1")
                        elif isoform_id not in final_isoforms:
                            final_isoforms.update({isoform_id:["1"]})
    print("Number of isoforms that pass isoform alignment filter")
    print(total_counts)
    return final_isoforms


#write final isoforms to a file
#file format: Isoform.ID\tNumber.of.CCS.Reads
def write():
    isoforms = compare_length_isoform()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Isoform.ID\tNumber.CCS.Reads\n"
        out.write(header)
        for isoform in isoforms:
            number_reads = len(isoforms[isoform])
            final = "%s\t%s\n" % (str(isoform), str(number_reads))
            out.write(final)

write()

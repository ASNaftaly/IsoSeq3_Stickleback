#creating final gtf for isoseq data
#gtf format will have 9 columns = chr num, source, feature, start, end, score, strand, frame, attributes
#source will be PacBio; feature can be gene, isoform (transcript), exon, CDS, start_codon, stop_codon, five_prime_utr, three_prime_utr; frame is only applicable to protein coding genes in the CDS
#non protein coding genes will only have gene, transcript, and exon features
#protein coding genes can have all features (though some may be missing)
#attributes will include (when applicable): gene id, transcript id, gene source, transcript source, gene biotype, transcript biotype, exon number, gene name
#if a value is not available, this will show up as "." if only value in column or be left blank if in the attribute section
#will need to read in several files to get all of the above positions and information
#determining start and stop codons, UTRs and CDS regions for 19 isoforms that have 5 or less amino acids in the predicted protein coding region
#to run script: python3 Create.Full.GTF.IsoformIssuesHandling.2.py <source as string; in this case all sources will be PacBio> <classification file for combined analyses> <19 isoforms only bed file> <full Ensembl GTF> <isoseq gtf with exon positions> <isoseq fasta file> <isoseq faa file> <start.stop positions file for 19 isoforms> <output file>
#Author: Alice Naftaly, Aug 2020

import sys
import time

#read classification file
#returns dictionary with key = isoform and value = [isoform, chr_num, strand, gene id, transcript id, protein coding]
def read_class():
    class_file = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                protein_coding = new_line[27]
                dict_value = [isoform, chr_num, strand, gene_id, transcript_id, protein_coding]
                class_dict.update({isoform:dict_value})
    return class_dict

#read bed file with all transcript start and stop positions
#returns dictionary with key == isoform and value == [isoform, chr_num, start_pos, end_pos, strand]
def read_shortproteins_bed():
    bed_file = sys.argv[3]
    bed_dict = {}
    with open(bed_file, 'r') as bed_info:
        for line in bed_info:
            new_line = line.split()
            chr_num = new_line[0]
            start_pos = new_line[1]
            end_pos = new_line[2]
            isoform = new_line[3]
            strand = new_line[5]
            dict_value = [isoform, chr_num, start_pos, end_pos, strand]
            bed_dict.update({isoform:dict_value})
    return bed_dict


#read in ensembl GTF
#will read in the ensembl gtf in a few ways
#this returns chr num, gene id, and gene name if available
def read_ensembl_gtf_genes():
    gtf_file = sys.argv[4]
    gene_gtf_dict = {}
    gene_name = ""
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            chr_num = new_line[0]
            feature = new_line[2]
            if feature == "gene":
                gene_info = new_line[8].split(";")
                for value in gene_info:
                    if value.startswith("gene_id"):
                        gene_id_full = value.split(" ")
                        gene_id = gene_id_full[1].strip("\"")
                    elif value.startswith(" gene_name"):
                        gene_name_full = value.split(" ")
                        gene_name = gene_name_full[2].strip("\"")
                if len(gene_name) == 0:
                    dict_value = [chr_num, gene_id, "."]
                elif len(gene_name) > 0:
                    dict_value = [chr_num, gene_id, gene_name]
                if gene_id in gene_gtf_dict:
                    gene_gtf_dict[gene_id].append(dict_value)
                elif gene_id not in gene_gtf_dict:
                    gene_gtf_dict.update({gene_id:[dict_value]})
    return gene_gtf_dict

#this returns chr num, gene id, transcript id, gene name if available
def read_ensembl_gtf_transcripts():
    gtf_file = sys.argv[4]
    transcripts_gtf_dict = {}
    gene_name = ""
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            chr_num = new_line[0]
            feature = new_line[2]
            if feature == "transcript":
                gene_info = new_line[8].split(";")
                for value in gene_info:
                    if value.startswith("gene_id"):
                        gene_id_full = value.split(" ")
                        gene_id = gene_id_full[1].strip("\"")
                    elif value.startswith(" gene_name"):
                        gene_name_full = value.split(" ")
                        gene_name = gene_name_full[2].strip("\"")
                    elif value.startswith(" transcript_id"):
                        transcript_id_full = value.split(" ")
                        transcript_id = transcript_id_full[2].strip("\"")
                if len(gene_name) == 0:
                    dict_value = [chr_num, gene_id, transcript_id, "."]
                elif len(gene_name) > 0:
                    dict_value = [chr_num, gene_id, transcript_id, gene_name]
                if transcript_id in transcripts_gtf_dict:
                    transcripts_gtf_dict[transcript_id].append(dict_value)
                elif transcript_id not in transcripts_gtf_dict:
                    transcripts_gtf_dict.update({transcript_id:[dict_value]})
    return transcripts_gtf_dict

#read in isoseq gtf with exon positions
#returns dictionary with key == isoform id and value == list of exons with [strand, exon start, exon end]
#flipped the order of exons for - strand
#first exon is 1bp off compared to bed file (need to subtract 1bp from the first exon start position)
def read_isoseq_gtf():
    gtf_file = sys.argv[5]
    exon_dict = {}
    final_exon_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            exon_start = new_line[3]
            exon_end = new_line[4]
            strand = new_line[6]
            isoform_id_full = new_line[9].strip(";")
            isoform_id = isoform_id_full.strip("\"")
            dict_value = [strand, exon_start, exon_end]
            if isoform_id in exon_dict:
                exon_dict[isoform_id].append(dict_value)
            elif isoform_id not in exon_dict:
                exon_dict.update({isoform_id:[dict_value]})
    for isoform in exon_dict:
        single_isoform = exon_dict[isoform]
        if len(single_isoform) == 1:
            single = single_isoform[0]
            final_single = [single[0], str(int(single[1])-1), single[2]]
            final_exon_dict.update({isoform:[final_single]})
        elif len(single_isoform) > 1:
            strand = single_isoform[0][0]
            if strand == "+":
                first_exon = single_isoform[0]
                final_first = [first_exon[0], str(int(first_exon[1])-1), first_exon[2]]
                single_isoform[0] = final_first
                final_exon_dict.update({isoform:single_isoform})
            elif strand == "-":
                single_isoform.reverse()
                first_exon = single_isoform[0]
                final_first = [first_exon[0], str(int(first_exon[1])-1), first_exon[2]]
                last_exon = single_isoform[len(single_isoform)-1]
                final_last = [last_exon[0], str(int(last_exon[1])-1), last_exon[2]]
                single_isoform[0] = final_first
                single_isoform[len(single_isoform)-1] = final_last
                final_exon_dict.update({isoform:single_isoform})
    return final_exon_dict

#read in nucleotides fasta file for isoforms
#returns dictionary with key == isoform id and value == sequence in a list for each isoform
def read_isoseq_nt_fasta():
    fasta_file = sys.argv[6]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                full_isoform_id = new_line[0].strip("\'")
                isoform_id = full_isoform_id.strip(">")
            else:
                new_line = line.strip("\n")
                if isoform_id in fasta_dict:
                    fasta_dict[isoform_id].append(new_line)
                elif isoform_id not in fasta_dict:
                    fasta_dict.update({isoform_id:[new_line]})
        for isoform in fasta_dict:
            final_seq = []
            single_isoform = fasta_dict[isoform]
            for seq in single_isoform:
                final_seq += seq
            final_fasta_dict.update({isoform:final_seq})
    return final_fasta_dict


#read in amino acid faa file for isoforms
#returns dictionary with key == isoform id and value == sequence in a list for each isoform
def read_isoseq_aa_faa():
    faa_file = sys.argv[7]
    faa_dict = {}
    final_faa_dict = {}
    with open(faa_file, 'r') as faa:
        for line in faa:
            if line.startswith(">"):
                new_line = line.split("\t")
                isoform_id = new_line[0].strip(">")
            else:
                new_line = line.strip("\n")
                if isoform_id in faa_dict:
                    faa_dict[isoform_id].append(new_line)
                elif isoform_id not in faa_dict:
                    faa_dict.update({isoform_id:[new_line]})
        for isoform in faa_dict:
            final_seq = []
            single_isoform = faa_dict[isoform]
            for seq in single_isoform:
                final_seq += seq
            final_faa_dict.update({isoform:final_seq})
    return final_faa_dict

#only keep the 19 weird isoforms
def filtered_class():
    weird_isoforms = read_shortproteins_bed()
    class_dict = read_class()
    final_class_dict = {}
    for isoform in class_dict:
        if isoform in weird_isoforms:
            final_class_dict.update({isoform:class_dict[isoform]})
    return final_class_dict

def filtered_isoseq_gtf():
    weird_isoforms = read_shortproteins_bed()
    gtf_dict = read_isoseq_gtf()
    final_gtf_dict = {}
    for isoform in gtf_dict:
        if isoform in weird_isoforms:
            final_gtf_dict.update({isoform:gtf_dict[isoform]})
    return final_gtf_dict

def filtered_isoseq_fasta():
    weird_isoforms = read_shortproteins_bed()
    fasta_dict = read_isoseq_nt_fasta()
    final_fasta_dict = {}
    for isoform in fasta_dict:
        if isoform in weird_isoforms:
            final_fasta_dict.update({isoform:fasta_dict[isoform]})
    return final_fasta_dict

def filtered_isoseq_faa():
    weird_isoforms = read_shortproteins_bed()
    faa_dict = read_isoseq_aa_faa()
    final_faa_dict = {}
    for isoform in faa_dict:
        if isoform in weird_isoforms:
            final_faa_dict.update({isoform:faa_dict[isoform]})
    return final_faa_dict

#amino acid conversion
#converting amino acid (single letter amino acid code) to codons (DNA sequence) in 5'-3' direction
#returns a dictionary with key == single letter amino acid code and value == list of codons in DNA in 5'-3'
def amino_acid_to_codons():
    aa_codon_dict = {}
    glycine = ["GGT", "GGC", "GGA", "GGG"]
    alanine = ["GCT", "GCC", "GCA", "GCG"]
    leucine = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
    methionine = ["ATG"]
    phenylalanine = ["TTT", "TTC"]
    tryptophan = ["TGG"]
    lysine = ["AAA", "AAG"]
    glutamine = ["CAA", "CAG"]
    glutamic_acid = ["GAA", "GAG"]
    serine = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
    proline = ["CCT", "CCC", "CCA", "CCG"]
    valine = ["GTT", "GTC", "GTA", "GTG"]
    isoleucine = ["ATT", "ATC", "ATA"]
    cysteine = ["TGT", "TGC"]
    tyrosine = ["TAT", "TAC"]
    histidine = ["CAT", "CAC"]
    arginine = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
    asparagine = ["AAT", "AAC"]
    aspartic_acid = ["GAT", "GAC"]
    threonine = ["ACT", "ACC", "ACA", "ACG"]
    stop = ["TAA", "TAG", "TGA"]
    aa_codon_dict.update({"G":glycine})
    aa_codon_dict.update({"A":alanine})
    aa_codon_dict.update({"L":leucine})
    aa_codon_dict.update({"M":methionine})
    aa_codon_dict.update({"F":phenylalanine})
    aa_codon_dict.update({"W":tryptophan})
    aa_codon_dict.update({"K":lysine})
    aa_codon_dict.update({"Q":glutamine})
    aa_codon_dict.update({"E":glutamic_acid})
    aa_codon_dict.update({"S":serine})
    aa_codon_dict.update({"P":proline})
    aa_codon_dict.update({"I":isoleucine})
    aa_codon_dict.update({"C":cysteine})
    aa_codon_dict.update({"Y":tyrosine})
    aa_codon_dict.update({"H":histidine})
    aa_codon_dict.update({"R":arginine})
    aa_codon_dict.update({"N":asparagine})
    aa_codon_dict.update({"D":aspartic_acid})
    aa_codon_dict.update({"T":threonine})
    aa_codon_dict.update({"V":valine})
    aa_codon_dict.update({"Stop":stop})
    return aa_codon_dict


#read in start and stop positions that I did manually
def read_start_stop_codons():
    codon_file = sys.argv[8]
    start_codon_dict = {}
    stop_codon_dict = {}
    with open(codon_file, 'r') as codons:
        for line in codons:
            new_line = line.split()
            isoform = new_line[0]
            start_codon_pos = [new_line[1], new_line[2]]
            stop_codon_pos = [new_line[3], new_line[4]]
            start_codon_dict.update({isoform:start_codon_pos})
            stop_codon_dict.update({isoform:stop_codon_pos})
    return start_codon_dict, stop_codon_dict

#convert start positions based on strand
#+ strand is fine
#- strand is in 5'-3' direction for nucleotide sequence, but need to flip position because poistions in gtf are 5'-3' for + strand
#there are 186 isoforms that the start codon spans the exon-intron boundary; these will be handled in a separate script/manually
def convert_start_positions():
    start_codons,stop_codons = read_start_stop_codons()
    exon_positions = filtered_isoseq_gtf()
    converted_start_codons = {}
    count = 0
    #key_list = ["PB.2804.1", "PB.8298.4","PB.5392.2"]
    for isoform in start_codons:
        single_start_codon = start_codons[isoform]
        single_codon_start = int(single_start_codon[0])
        single_codon_end = int(single_start_codon[1])
        frame = single_start_codon[1]
        exon_list = exon_positions[isoform]
        strand = exon_list[0][0]
        if strand == "+":
            if len(exon_list) == 1:
                exon = exon_list[0]
                exon_start = int(exon[1])
                exon_end = int(exon[2])
                converted_start_codon_start = exon_start + single_codon_start
                converted_start_codon_end = converted_start_codon_start + 2
                new_start_codon_pos = [str(converted_start_codon_start), str(converted_start_codon_end)]
                converted_start_codons.update({isoform:new_start_codon_pos})
            elif len(exon_list) > 1:
                exon_1 = exon_list[0]
                exon_1_start = int(exon_1[1])
                exon_1_end = int(exon_1[2])
                exon_1_length = exon_1_end - exon_1_start
                if exon_1_length > single_codon_start:
                    #print("start codon in first exon")
                    converted_start_codon_start = exon_1_start + single_codon_start
                    converted_start_codon_end = converted_start_codon_start + 2
                    new_start_codon_pos = [str(converted_start_codon_start), str(converted_start_codon_end)]
                    converted_start_codons.update({isoform:new_start_codon_pos})
                #this is a strange scenario where the start codon actually spans exon-intron boundary
                #will probably want to specify this in the final gtf file
                elif exon_1_length == single_codon_start:
                    continue
                elif exon_1_length < single_codon_start:
                    #print("start codon not in first exon")
                    exon_length_list = []
                    for exon in exon_list:
                        exon_start = int(exon[1])
                        exon_end = int(exon[2])
                        exon_length = exon_end - exon_start
                        exon_length_list.append(exon_length)
                    x = 0
                    moving_start_codon = single_codon_start
                    while x < len(exon_length_list):
                        if moving_start_codon - exon_length_list[x] > 0:
                            moving_start_codon = moving_start_codon - exon_length_list[x]
                            x += 1
                        elif moving_start_codon - exon_length_list[x] <= 0:
                            converted_start = int(exon_list[x][1]) + moving_start_codon
                            #there is some strange counting/math issue here that you need to subtract the exon number from the converted start
                            #if I figure out what this is; I'll explain it here or try to code for it better
                            final_converted_start = converted_start - x
                            final_converted_end = final_converted_start + 2
                            new_start_codon_pos = [str(final_converted_start), str(final_converted_end)]
                            converted_start_codons.update({isoform:new_start_codon_pos})
                            break
        elif strand == "-":
            if len(exon_list) == 1:
                exon = exon_list[0]
                exon_start = int(exon[2])
                exon_end = int(exon[1])
                converted_start_codon_start = exon_start - single_codon_start
                converted_start_codon_end = converted_start_codon_start - 2
                new_start_codon_pos = [str(converted_start_codon_start), str(converted_start_codon_end)]
                converted_start_codons.update({isoform:new_start_codon_pos})
            elif len(exon_list) > 1:
                exon_1 = exon_list[0]
                exon_1_start = int(exon_1[2])
                exon_1_end = int(exon_1[1])
                exon_1_length = exon_1_start - exon_1_end
                if exon_1_length > single_codon_start:
                    converted_start_codon_start = exon_1_start - single_codon_start + 1
                    converted_start_codon_end = converted_start_codon_start - 2
                    new_start_codon_pos = [str(converted_start_codon_start), str(converted_start_codon_end)]
                    converted_start_codons.update({isoform:new_start_codon_pos})
                #this is a strange scenario where the start codon actually spans exon-intron boundary
                #will probably want to specify this in the final gtf file
                elif exon_1_length == single_codon_start:
                    continue
                elif exon_1_length < single_codon_start:
                    exon_length_list = []
                    for exon in exon_list:
                        exon_start = int(exon[2])
                        exon_end = int(exon[1])
                        exon_length = exon_start - exon_end
                        exon_length_list.append(exon_length)
                    x = 0
                    moving_start_codon = single_codon_start
                    while x < len(exon_length_list):
                        if moving_start_codon - exon_length_list[x] > 0:
                            moving_start_codon = moving_start_codon - exon_length_list[x]
                            x += 1
                        elif moving_start_codon - exon_length_list[x] <= 0:
                            converted_start = int(exon_list[x][2]) - moving_start_codon
                            #there is some strange counting/math issue here that you need to add the exon number from the converted start
                            #if I figure out what this is; I'll explain it here or try to code for it better
                            final_converted_start = converted_start + x - 1
                            final_converted_end = final_converted_start - 2
                            new_start_codon_pos = [str(final_converted_start), str(final_converted_end)]
                            converted_start_codons.update({isoform:new_start_codon_pos})
                            break
    return converted_start_codons

#convert stop positions based on strand
#+ strand is fine
#- strand is in 5'-3' direction for nucleotide sequence, but need to flip position because poistions in gtf are 5'-3' for + strand
def convert_stop_positions():
    start_codons,stop_codons = read_start_stop_codons()
    converted_start_positions = convert_start_positions()
    exon_positions = filtered_isoseq_gtf()
    aa_sequences = filtered_isoseq_faa()
    converted_stop_codons = {}
    for isoform in stop_codons:
        single_start_position = converted_start_positions[isoform]
        single_protein_seq = aa_sequences[isoform]
        single_stop_codon = stop_codons[isoform]
        single_codon_start = int(single_stop_codon[0])
        single_codon_end = int(single_stop_codon[1])
        exon_list = exon_positions[isoform]
        strand = exon_list[0][0]
        if strand == "+":
            if len(exon_list) == 1:
                exon = exon_list[0]
                exon_start = int(exon[1])
                exon_end = int(exon[2])
                converted_stop_codon_start = exon_start + single_codon_start
                converted_stop_codon_end = converted_stop_codon_start + 2
                new_stop_codon_pos = [str(converted_stop_codon_start), str(converted_stop_codon_end)]
                converted_stop_codons.update({isoform:new_stop_codon_pos})
            elif len(exon_list) > 1:
                exon_1 = exon_list[0]
                exon_1_start = int(exon_1[1])
                exon_1_end = int(exon_1[2])
                exon_1_length = exon_1_end - exon_1_start
                if exon_1_length > single_codon_start:
                    converted_stop_codon_start = exon_1_start + single_codon_start
                    converted_stop_codon_end = converted_stop_codon_start + 2
                    new_stop_codon_pos = [str(converted_stop_codon_start), str(converted_stop_codon_end)]
                    converted_stop_codons.update({isoform:new_stop_codon_pos})
                elif exon_1_length < single_codon_start:
                    exon_length_list = []
                    for exon in exon_list:
                        exon_start = int(exon[1])
                        exon_end = int(exon[2])
                        exon_length = exon_end - exon_start
                        exon_length_list.append(exon_length)
                    x = 0
                    moving_stop_codon = single_codon_start
                    while x < len(exon_length_list):
                        if moving_stop_codon - exon_length_list[x] > 0:
                            moving_stop_codon = moving_stop_codon - exon_length_list[x]
                            x += 1
                        elif moving_stop_codon - exon_length_list[x] <= 0:
                            converted_stop = int(exon_list[x][1]) + moving_stop_codon
                            protein_length = len(single_protein_seq)
                            seq_length = (converted_stop - int(single_start_position[0]))
                            if str(int(seq_length/3)).endswith(".0"):
                                if protein_length > 1:
                                    final_converted_stop = converted_stop
                                    final_converted_end = final_converted_stop + 2
                                elif protein_length == 1:
                                    final_converted_stop = converted_stop - x
                                    final_converted_end = final_converted_stop + 2
                            elif str(int(seq_length-1/3)).endswith(".0"):
                                if protein_length > 1:
                                    final_converted_stop = converted_stop
                                    final_converted_end = final_converted_stop + 2
                                elif protein_length == 1:
                                    final_converted_stop = converted_stop - x
                                    final_converted_end = final_converted_stop + 2
                            if str(int(seq_length-2/3)).endswith(".0"):
                                if protein_length > 1:
                                    final_converted_stop = converted_stop
                                    final_converted_end = final_converted_stop + 2
                                elif protein_length == 1:
                                    final_converted_stop = converted_stop - x
                                    final_converted_end = final_converted_stop + 2
                            else:
                                #there is some strange counting/math issue here that you need to subtract the exon number from the converted start
                                #if I figure out what this is; I'll explain it here or try to code for it better
                                final_converted_stop = converted_stop - x
                                final_converted_end = final_converted_stop + 2
                            new_stop_codon_pos = [str(final_converted_stop), str(final_converted_end)]
                            converted_stop_codons.update({isoform:new_stop_codon_pos})
                            break
        elif strand == "-":
            if len(exon_list) == 1:
                exon = exon_list[0]
                exon_start = int(exon[2])
                exon_end = int(exon[1])
                converted_stop_codon_start = exon_start - single_codon_start
                converted_stop_codon_end = converted_stop_codon_start - 2
                seq_length =  int(single_start_position[0]) - converted_stop_codon_start
                new_stop_codon_pos = [str(converted_stop_codon_start), str(converted_stop_codon_end)]
                converted_stop_codons.update({isoform:new_stop_codon_pos})
            elif len(exon_list) > 1:
                exon_1 = exon_list[0]
                exon_1_start = int(exon_1[2])
                exon_1_end = int(exon_1[1])
                exon_1_length = exon_1_start - exon_1_end
                if exon_1_length > single_codon_start:
                    converted_stop_codon_start = exon_1_start - single_codon_start
                    converted_stop_codon_end = converted_stop_codon_start - 2
                    new_stop_codon_pos = [str(converted_stop_codon_start), str(converted_stop_codon_end)]
                    converted_stop_codons.update({isoform:new_stop_codon_pos})
                elif exon_1_length < single_codon_start:
                    exon_length_list = []
                    for exon in exon_list:
                        exon_start = int(exon[2])
                        exon_end = int(exon[1])
                        exon_length = exon_start - exon_end
                        exon_length_list.append(exon_length)
                    x = 0
                    moving_stop_codon = single_codon_start
                    while x < len(exon_length_list):
                        if moving_stop_codon - exon_length_list[x] > 0:
                            moving_stop_codon = moving_stop_codon - exon_length_list[x]
                            x += 1
                        elif moving_stop_codon - exon_length_list[x] <= 0:
                            converted_stop = int(exon_list[x][2]) - moving_stop_codon
                            protein_length = len(single_protein_seq)
                            seq_length = str(( int(single_start_position[0]) - converted_stop - 2)/3)
                            if seq_length.endswith(".0"):
                                final_converted_stop = converted_stop
                                final_converted_end = final_converted_stop + 2
                            else:
                            #there is some strange counting/math issue here that you need to add the exon number from the converted start
                            #if I figure out what this is; I'll explain it here or try to code for it better
                                final_converted_stop = converted_stop + x - 1
                                final_converted_end = final_converted_stop - 2
                            new_stop_codon_pos = [str(final_converted_stop), str(final_converted_end)]
                            converted_stop_codons.update({isoform:new_stop_codon_pos})
                            break
    return converted_stop_codons


#determine CDS positions, 5' UTRs, and 3'UTRs
def determine_cds_utrs_positions():
    exon_dict = filtered_isoseq_gtf()
    start_codons = convert_start_positions()
    stop_codons = convert_stop_positions()
    five_prime_utr_dict = {}
    cds_dict = {}
    three_prime_utr_dict = {}
    for isoform in exon_dict:
        if isoform in start_codons and isoform in stop_codons:
            single_exon_list = exon_dict[isoform]
            single_start_position = start_codons[isoform]
            single_stop_position = stop_codons[isoform]
            strand = single_exon_list[0][0]
            if strand == "+":
                if len(single_exon_list) == 1:
                    exon_start = int(single_exon_list[0][1])
                    exon_end = int(single_exon_list[0][2])
                    start_codon_start_pos = int(single_start_position[0])
                    start_codon_end_pos = int(single_start_position[1])
                    stop_codon_start_pos = int(single_stop_position[0])
                    stop_codon_end_pos = int(single_stop_position[1])
                    five_prime_utr_pos = [str(exon_start), str(start_codon_start_pos - 1)]
                    cds_pos = [str(start_codon_start_pos), str(stop_codon_end_pos), "1"]
                    three_prime_utr_pos = [str(stop_codon_end_pos+1), str(exon_end)]
                    five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                    cds_dict.update({isoform:[cds_pos]})
                    three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                elif len(single_exon_list) > 1:
                    for index, exon in enumerate(single_exon_list):
                        exon_start = int(exon[1])
                        exon_end = int(exon[2])
                        exon_pos = [exon_start, exon_end]
                        start_codon_start_pos = int(single_start_position[0])
                        start_codon_end_pos = int(single_start_position[1])
                        stop_codon_start_pos = int(single_stop_position[0])
                        stop_codon_end_pos = int(single_stop_position[1])
                        if start_codon_start_pos > exon_end:
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(exon_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[exon_pos]})
                        elif exon_start < start_codon_start_pos < exon_end and stop_codon_start_pos > exon_end:
                            end_five_prime_utr = start_codon_start_pos - 1
                            five_prime_utr_pos = [exon_start, end_five_prime_utr]
                            cds_start_pos = start_codon_start_pos
                            cds_pos = [cds_start_pos, exon_end, str(index+1)]
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(five_prime_utr_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                        elif exon_start < start_codon_start_pos < exon_end and stop_codon_start_pos < exon_end:
                            end_five_prime_utr = start_codon_start_pos - 1
                            five_prime_utr_pos = [exon_start, end_five_prime_utr]
                            cds_start_pos = start_codon_start_pos
                            cds_pos = [cds_start_pos, stop_codon_end_pos, str(index+1)]
                            start_three_prime_utr = stop_codon_end_pos + 1
                            three_prime_utr_pos = [start_three_prime_utr, exon_end]
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(five_prime_utr_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(three_prime_utr_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                        elif stop_codon_start_pos > exon_end and start_codon_start_pos < exon_end:
                            final = [exon_start, exon_end, str(index+1)]
                            if isoform in cds_dict:
                                cds_dict[isoform].append(final)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[final]})
                        elif exon_start < stop_codon_start_pos < exon_end:
                            cds_stop_pos = stop_codon_end_pos
                            cds_pos = [exon_start, cds_stop_pos, str(index+1)]
                            three_prime_utr_start_pos = stop_codon_end_pos + 1
                            three_prime_utr_pos = [three_prime_utr_start_pos, exon_end]
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(three_prime_utr_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                        elif stop_codon_start_pos < exon_end and stop_codon_start_pos < exon_start:
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(exon_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[exon_pos]})
            elif strand == "-":
                if len(single_exon_list) == 1:
                    exon_start = int(single_exon_list[0][2])
                    exon_end = int(single_exon_list[0][1])
                    start_codon_start_pos = int(single_start_position[0])
                    start_codon_end_pos = int(single_start_position[1])
                    stop_codon_start_pos = int(single_stop_position[0])
                    stop_codon_end_pos = int(single_stop_position[1])
                    five_prime_utr_pos = [str(exon_start), str(start_codon_start_pos + 1)]
                    cds_pos = [str(start_codon_start_pos), str(stop_codon_end_pos), "1"]
                    three_prime_utr_pos = [str(stop_codon_end_pos-1), str(exon_end)]
                    five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                    cds_dict.update({isoform:[cds_pos]})
                    three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                elif len(single_exon_list) > 1:
                    for ind, exon in enumerate(single_exon_list):
                        exon_start = int(exon[2])
                        exon_end = int(exon[1])
                        exon_pos = [exon_start, exon_end]
                        start_codon_start_pos = int(single_start_position[0])
                        start_codon_end_pos = int(single_start_position[1])
                        stop_codon_start_pos = int(single_stop_position[0])
                        stop_codon_end_pos = int(single_stop_position[1])
                        if start_codon_start_pos < exon_end:
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(exon_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[exon_pos]})
                        elif exon_start > start_codon_start_pos > exon_end and stop_codon_start_pos < exon_end:
                            end_five_prime_utr = start_codon_start_pos + 1
                            five_prime_utr_pos = [exon_start, end_five_prime_utr]
                            cds_start_pos = start_codon_start_pos
                            cds_pos = [cds_start_pos, exon_end, str(ind+1)]
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(five_prime_utr_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                        elif exon_start > start_codon_start_pos > exon_end and stop_codon_start_pos > exon_end:
                            end_five_prime_utr = start_codon_start_pos + 1
                            five_prime_utr_pos = [exon_start, end_five_prime_utr]
                            cds_start_pos = start_codon_start_pos
                            cds_end_pos = stop_codon_end_pos
                            cds_pos = [cds_start_pos, cds_end_pos, str(ind+1)]
                            three_prime_utr_pos = [str(stop_codon_end_pos-1), exon_end]
                            if isoform in five_prime_utr_dict:
                                five_prime_utr_dict[isoform].append(five_prime_utr_pos)
                            elif isoform not in five_prime_utr_dict:
                                five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(three_prime_utr_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                        elif stop_codon_start_pos < exon_end and start_codon_start_pos > exon_end:
                            final = [exon_start, exon_end, str(ind+1)]
                            if isoform in cds_dict:
                                cds_dict[isoform].append(final)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[final]})
                        elif exon_start > stop_codon_start_pos > exon_end:
                            cds_stop_pos = stop_codon_end_pos
                            cds_pos = [exon_start, cds_stop_pos, str(ind+1)]
                            three_prime_utr_start_pos = stop_codon_end_pos -1
                            three_prime_utr_pos = [three_prime_utr_start_pos, exon_end]
                            if isoform in cds_dict:
                                cds_dict[isoform].append(cds_pos)
                            elif isoform not in cds_dict:
                                cds_dict.update({isoform:[cds_pos]})
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(three_prime_utr_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
                        elif stop_codon_end_pos > exon_start:
                            three_prime_utr_pos = [exon_start, exon_end]
                            if isoform in three_prime_utr_dict:
                                three_prime_utr_dict[isoform].append(three_prime_utr_pos)
                            elif isoform not in three_prime_utr_dict:
                                three_prime_utr_dict.update({isoform:[three_prime_utr_pos]})
        #pulling isoform data without a stop codon
        if isoform in start_codons and isoform not in stop_codons:
            single_exon_list = exon_dict[isoform]
            single_start_position = start_codons[isoform]
            strand = single_exon_list[0][0]
            #strand is + and there is only 1 exon for this isoform
            five_prime_utr_start = int(single_exon_list[0][1])
            five_prime_utr_end = int(single_start_position[0]) - 1
            five_prime_utr_pos = [five_prime_utr_start, five_prime_utr_end]
            cds_pos = [single_start_position[0], single_exon_list[0][2], "1"]
            five_prime_utr_dict.update({isoform:[five_prime_utr_pos]})
            cds_dict.update({isoform:[cds_pos]})
    return five_prime_utr_dict, cds_dict, three_prime_utr_dict

#create transcript feature
#this will include:
#chr num, source=PacBio, feature=transcript, start, end, score=".", strand, frame = ".", attributes = gene id, gene source, gene biotype, gene name if available, transcript id, transcript source, transcript biotype; isoform id
#difference between transcript and isoform id: transcript id will be either ENSGACT or novel and isoform id is PB.XXXX.X
def create_transcript_feature():
    class_dict = filtered_class()
    bed_dict = read_shortproteins_bed()
    ensembl_dict = read_ensembl_gtf_transcripts()
    source = sys.argv[1]
    transcript_feature_dict = {}
    for isoform in class_dict:
        single_isoform_class = class_dict[isoform]
        single_bed_isoform = bed_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        coding_potential = single_isoform_class[5]
        if coding_potential == "coding":
            final_potential = "protein_coding"
        elif coding_potential == "non_coding":
            final_potential = "nonprotein_coding"
        isoform_start = single_bed_isoform[2]
        isoform_end = single_bed_isoform[3]
        if transcript_id in ensembl_dict:
            single_ensembl = ensembl_dict[transcript_id][0]
            if single_ensembl[3] == ".":
                transcript_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            else:
                gene_name = single_ensembl[3]
                transcript_attributes = "gene_id %s; gene_source %s; gene_biotype %s; gene_name %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s;" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
            transcript_feature = [str(chr_num), str(source), "transcript", str(isoform_start), str(isoform_end), ".", str(strand), ".", transcript_attributes]
            if isoform in transcript_feature_dict:
                transcript_feature_dict[isoform].append(transcript_feature)
            elif isoform not in transcript_feature_dict:
                transcript_feature_dict.update({isoform:[transcript_feature]})
        elif transcript_id not in ensembl_dict:
            transcript_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            transcript_feature = [str(chr_num), str(source), "transcript", str(isoform_start), str(isoform_end), ".", str(strand), ".", transcript_attributes]
            if isoform in transcript_feature_dict:
                transcript_feature_dict[isoform].append(transcript_feature)
            elif isoform not in transcript_feature_dict:
                transcript_feature_dict.update({isoform:[transcript_feature]})
    return transcript_feature_dict

#create exon feature
#exon feature will contain
#chr num, source=PacBio, feature=exon, start, end, score=".", strand, frame = ".", attributes = gene id, gene source, gene biotype, gene name if available, transcript id, transcript source, transcript biotype; isoform id; exon number
def create_exon_feature():
    class_dict = filtered_class()
    ensembl_dict = read_ensembl_gtf_genes()
    exon_dict = filtered_isoseq_gtf()
    source = sys.argv[1]
    exon_feature_dict = {}
    for isoform in class_dict:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        coding_potential = single_isoform_class[5]
        if coding_potential == "coding":
            final_potential = "protein_coding"
        elif coding_potential == "non_coding":
            final_potential = "nonprotein_coding"
        exon_list = exon_dict[isoform]
        exon_count = 1
        if len(exon_list) == 1:
            single_exon = exon_list[0]
            exon_start = single_exon[1]
            exon_end = single_exon[2]
            if gene_id in ensembl_dict:
                single_ensembl = ensembl_dict[gene_id][0]
                if single_ensembl[2] == ".":
                    exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                else:
                    gene_name = single_ensembl[2]
                    exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; gene_name %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                exon_feature = [str(chr_num), str(source), "exon", str(exon_start), str(exon_end), ".", str(strand), ".", exon_attributes]
                if isoform in exon_feature_dict:
                    exon_feature_dict[isoform].append(exon_feature)
                elif isoform not in exon_feature_dict:
                    exon_feature_dict.update({isoform:[exon_feature]})
            elif gene_id not in ensembl_dict:
                exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                exon_feature = [str(chr_num), str(source), "exon", str(exon_start), str(exon_end), ".", str(strand), ".", exon_attributes]
                if isoform in exon_feature_dict:
                    exon_feature_dict[isoform].append(exon_feature)
                elif isoform not in exon_feature_dict:
                    exon_feature_dict.update({isoform:[exon_feature]})
        elif len(exon_list) > 1:
            for exon in exon_list:
                exon_start = exon[1]
                exon_end = exon[2]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                    else:
                        gene_name = single_ensembl[2]
                        exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; gene_name %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                    exon_feature = [str(chr_num), str(source), "exon", str(exon_start), str(exon_end), ".", str(strand), ".", exon_attributes]
                    if isoform in exon_feature_dict:
                        exon_feature_dict[isoform].append(exon_feature)
                    elif isoform not in exon_feature_dict:
                        exon_feature_dict.update({isoform:[exon_feature]})
                elif gene_id not in ensembl_dict:
                    exon_attributes = "gene_id %s; gene_source %s; gene_biotype %s; transcript_id %s; transcript_source %s; transcript_biotype %s; isoform_id %s; exon_number %s;" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_count))
                    exon_feature = [str(chr_num), str(source), "exon", str(exon_start), str(exon_end), ".", str(strand), ".", exon_attributes]
                    if isoform in exon_feature_dict:
                        exon_feature_dict[isoform].append(exon_feature)
                    elif isoform not in exon_feature_dict:
                        exon_feature_dict.update({isoform:[exon_feature]})
                exon_count += 1
    return exon_feature_dict

#create start codon feature
def create_start_codon_feature():
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    start_codons = convert_start_positions()
    source = sys.argv[1]
    start_codon_feature_dict = {}
    for isoform in start_codons:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        final_potential = "protein_coding"
        single_start_codon = start_codons[isoform]
        if strand == "+":
            codon_start = single_start_codon[0]
            codon_end = single_start_codon[1]
        elif strand == "-":
            codon_start = single_start_codon[1]
            codon_end = single_start_codon[0]
        if gene_id in ensembl_dict:
            single_ensembl = ensembl_dict[gene_id][0]
            if single_ensembl[2] == ".":
                start_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            else:
                gene_name = single_ensembl[2]
                start_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
            start_codon_feature = [str(chr_num), str(source), "start_codon", str(codon_start), str(codon_end), ".", str(strand), ".", start_codon_attributes]
            if isoform in start_codon_feature_dict:
                start_codon_feature_dict[isoform].append(start_codon_feature)
            elif isoform not in start_codon_feature_dict:
                start_codon_feature_dict.update({isoform:[start_codon_feature]})
        elif gene_id not in ensembl_dict:
            start_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            start_codon_feature = [str(chr_num), str(source), "start_codon", str(codon_start), str(codon_end), ".", str(strand), ".", start_codon_attributes]
            if isoform in start_codon_feature_dict:
                start_codon_feature_dict[isoform].append(start_codon_feature)
            elif isoform not in start_codon_feature_dict:
                start_codon_feature_dict.update({isoform:[start_codon_feature]})
    return start_codon_feature_dict


#create stop codon feature
def create_stop_codon_feature():
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    stop_codons = convert_stop_positions()
    source = sys.argv[1]
    stop_codon_feature_dict = {}
    for isoform in stop_codons:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        coding_potential = single_isoform_class[5]
        final_potential = "protein_coding"
        single_stop_codon = stop_codons[isoform]
        if strand == "+":
            codon_start = single_stop_codon[0]
            codon_end = single_stop_codon[1]
        elif strand == "-":
            codon_start = single_stop_codon[1]
            codon_end = single_stop_codon[0]
        if gene_id in ensembl_dict:
            single_ensembl = ensembl_dict[gene_id][0]
            if single_ensembl[2] == ".":
                stop_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            else:
                gene_name = single_ensembl[2]
                stop_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
            stop_codon_feature = [str(chr_num), str(source), "stop_codon", str(codon_start), str(codon_end), ".", str(strand), ".", stop_codon_attributes]
            if isoform in stop_codon_feature_dict:
                stop_codon_feature_dict[isoform].append(stop_codon_feature)
            elif isoform not in stop_codon_feature_dict:
                stop_codon_feature_dict.update({isoform:[stop_codon_feature]})
        elif gene_id not in ensembl_dict:
            stop_codon_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
            stop_codon_feature = [str(chr_num), str(source), "stop_codon", str(codon_start), str(codon_end), ".", str(strand), ".", stop_codon_attributes]
            if isoform in stop_codon_feature_dict:
                stop_codon_feature_dict[isoform].append(stop_codon_feature)
            elif isoform not in stop_codon_feature_dict:
                stop_codon_feature_dict.update({isoform:[stop_codon_feature]})
    return stop_codon_feature_dict

#create five prime utr feature
def create_five_prime_UTR_feature():
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    five_prime_utrs, cds_regions, three_prime_utrs = determine_cds_utrs_positions()
    source = sys.argv[1]
    five_prime_utr_feature_dict = {}
    for isoform in five_prime_utrs:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        final_potential = "protein_coding"
        single_five_prime_utr = five_prime_utrs[isoform]
        if strand == "+":
            if len(single_five_prime_utr) == 1:
                five_prime_utr_start = single_five_prime_utr[0][0]
                five_prime_utr_end = single_five_prime_utr[0][1]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    else:
                        gene_name = single_ensembl[2]
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                    five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                    if isoform in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                    elif isoform not in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
                elif gene_id not in ensembl_dict:
                    five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                    if isoform in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                    elif isoform not in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
            elif len(single_five_prime_utr) > 1:
                for single in single_five_prime_utr:
                    five_prime_utr_start = single[0]
                    five_prime_utr_end = single[1]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        else:
                            gene_name = single_ensembl[2]
                            five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                        five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                        if isoform in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                        elif isoform not in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
                    elif gene_id not in ensembl_dict:
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                        if isoform in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                        elif isoform not in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
        elif strand == "-":
            if len(single_five_prime_utr) == 1:
                five_prime_utr_start = single_five_prime_utr[0][1]
                five_prime_utr_end = single_five_prime_utr[0][0]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    else:
                        gene_name = single_ensembl[2]
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                    five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                    if isoform in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                    elif isoform not in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
                elif gene_id not in ensembl_dict:
                    five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                    if isoform in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                    elif isoform not in five_prime_utr_feature_dict:
                        five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
            elif len(single_five_prime_utr) > 1:
                for single in single_five_prime_utr:
                    five_prime_utr_start = single[1]
                    five_prime_utr_end = single[0]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        else:
                            gene_name = single_ensembl[2]
                            five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                        five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                        if isoform in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                        elif isoform not in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
                    elif gene_id not in ensembl_dict:
                        five_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        five_prime_utr_feature = [str(chr_num), str(source), "five_prime_utr", str(five_prime_utr_start), str(five_prime_utr_end), ".", str(strand), ".", five_prime_utr_attributes]
                        if isoform in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict[isoform].append(five_prime_utr_feature)
                        elif isoform not in five_prime_utr_feature_dict:
                            five_prime_utr_feature_dict.update({isoform:[five_prime_utr_feature]})
    return five_prime_utr_feature_dict

#create three prime utr feature
def create_three_prime_UTR_feature():
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    three_prime_utrs, cds_regions, three_prime_utrs = determine_cds_utrs_positions()
    source = sys.argv[1]
    three_prime_utr_feature_dict = {}
    for isoform in three_prime_utrs:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        final_potential = "protein_coding"
        single_three_prime_utr = three_prime_utrs[isoform]
        if strand == "+":
            if len(single_three_prime_utr) == 1:
                three_prime_utr_start = single_three_prime_utr[0][0]
                three_prime_utr_end = single_three_prime_utr[0][1]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    else:
                        gene_name = single_ensembl[2]
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                    three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                    if isoform in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                    elif isoform not in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
                elif gene_id not in ensembl_dict:
                    three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                    if isoform in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                    elif isoform not in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
            elif len(single_three_prime_utr) > 1:
                for single in single_three_prime_utr:
                    three_prime_utr_start = single[0]
                    three_prime_utr_end = single[1]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        else:
                            gene_name = single_ensembl[2]
                            three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                        three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                        if isoform in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                        elif isoform not in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
                    elif gene_id not in ensembl_dict:
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                        if isoform in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                        elif isoform not in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
        elif strand == "-":
            if len(single_three_prime_utr) == 1:
                three_prime_utr_start = single_three_prime_utr[0][1]
                three_prime_utr_end = single_three_prime_utr[0][0]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    else:
                        gene_name = single_ensembl[2]
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                    three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                    if isoform in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                    elif isoform not in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
                elif gene_id not in ensembl_dict:
                    three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                    three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                    if isoform in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                    elif isoform not in three_prime_utr_feature_dict:
                        three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
            elif len(single_three_prime_utr) > 1:
                for single in single_three_prime_utr:
                    three_prime_utr_start = single[1]
                    three_prime_utr_end = single[0]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        else:
                            gene_name = single_ensembl[2]
                            three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform))
                        three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                        if isoform in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                        elif isoform not in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
                    elif gene_id not in ensembl_dict:
                        three_prime_utr_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform))
                        three_prime_utr_feature = [str(chr_num), str(source), "three_prime_utr", str(three_prime_utr_start), str(three_prime_utr_end), ".", str(strand), ".", three_prime_utr_attributes]
                        if isoform in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict[isoform].append(three_prime_utr_feature)
                        elif isoform not in three_prime_utr_feature_dict:
                            three_prime_utr_feature_dict.update({isoform:[three_prime_utr_feature]})
    return three_prime_utr_feature_dict

#create cds feature
def create_CDS_feature():
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    five_prime_utrs, cds_regions, three_prime_utrs = determine_cds_utrs_positions()
    source = sys.argv[1]
    cds_feature_dict = {}
    for isoform in cds_regions:
        single_isoform_class = class_dict[isoform]
        chr_num = single_isoform_class[1]
        strand = single_isoform_class[2]
        gene_id = single_isoform_class[3]
        transcript_id = single_isoform_class[4]
        final_potential = "protein_coding"
        single_cds_region = cds_regions[isoform]
        if strand == "+":
            if len(single_cds_region) == 1:
                cds_region_start = single_cds_region[0][0]
                cds_region_end = single_cds_region[0][1]
                exon_number = single_cds_region[0][2]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    else:
                        gene_name = single_ensembl[2]
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                    if isoform in cds_feature_dict:
                        cds_feature_dict[isoform].append(cds_region_feature)
                    elif isoform not in cds_feature_dict:
                        cds_feature_dict.update({isoform:[cds_region_feature]})
                elif gene_id not in ensembl_dict:
                    cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                    if isoform in cds_feature_dict:
                        cds_feature_dict[isoform].append(cds_region_feature)
                    elif isoform not in cds_feature_dict:
                        cds_feature_dict.update({isoform:[cds_region_feature]})
            elif len(single_cds_region) > 1:
                for single in single_cds_region:
                    cds_region_start = single[0]
                    cds_region_end = single[1]
                    exon_number = single[2]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        else:
                            gene_name = single_ensembl[2]
                            cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                        if isoform in cds_feature_dict:
                            cds_feature_dict[isoform].append(cds_region_feature)
                        elif isoform not in cds_feature_dict:
                            cds_feature_dict.update({isoform:[cds_region_feature]})
                    elif gene_id not in ensembl_dict:
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                        if isoform in cds_feature_dict:
                            cds_feature_dict[isoform].append(cds_region_feature)
                        elif isoform not in cds_feature_dict:
                            cds_feature_dict.update({isoform:[cds_region_feature]})
        elif strand == "-":
            if len(single_cds_region) == 1:
                cds_region_start = single_cds_region[0][1]
                cds_region_end = single_cds_region[0][0]
                exon_number = single_cds_region[0][2]
                if gene_id in ensembl_dict:
                    single_ensembl = ensembl_dict[gene_id][0]
                    if single_ensembl[2] == ".":
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    else:
                        gene_name = single_ensembl[2]
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                    if isoform in cds_feature_dict:
                        cds_feature_dict[isoform].append(cds_region_feature)
                    elif isoform not in cds_feature_dict:
                        cds_feature_dict.update({isoform:[cds_region_feature]})
                elif gene_id not in ensembl_dict:
                    cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                    cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                    if isoform in cds_feature_dict:
                        cds_feature_dict[isoform].append(cds_region_feature)
                    elif isoform not in cds_feature_dict:
                        cds_feature_dict.update({isoform:[cds_region_feature]})
            elif len(single_cds_region) > 1:
                for single in single_cds_region:
                    cds_region_start = single[0]
                    cds_region_end = single[1]
                    exon_number = single[2]
                    if gene_id in ensembl_dict:
                        single_ensembl = ensembl_dict[gene_id][0]
                        if single_ensembl[2] == ".":
                            cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        else:
                            gene_name = single_ensembl[2]
                            cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; gene_name '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(gene_name), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                        if isoform in cds_feature_dict:
                            cds_feature_dict[isoform].append(cds_region_feature)
                        elif isoform not in cds_feature_dict:
                            cds_feature_dict.update({isoform:[cds_region_feature]})
                    elif gene_id not in ensembl_dict:
                        cds_region_attributes = "gene_id '%s'; gene_source '%s'; gene_biotype '%s'; transcript_id '%s'; transcript_source '%s'; transcript_biotype '%s'; isoform_id '%s'; exon_number '%s';" % (str(gene_id), str(source), str(final_potential), str(transcript_id), str(source), str(final_potential), str(isoform), str(exon_number))
                        cds_region_feature = [str(chr_num), str(source), "CDS", str(cds_region_start), str(cds_region_end), ".", str(strand), ".", cds_region_attributes]
                        if isoform in cds_feature_dict:
                            cds_feature_dict[isoform].append(cds_region_feature)
                        elif isoform not in cds_feature_dict:
                            cds_feature_dict.update({isoform:[cds_region_feature]})
    return cds_feature_dict


#writing output to files
#Will write all features to 1 file that will be sorted later
def write():
    transcript_features = create_transcript_feature()
    exon_features = create_exon_feature()
    start_codon_features = create_start_codon_feature()
    stop_codon_features = create_stop_codon_feature()
    five_prime_utr_features = create_five_prime_UTR_feature()
    three_prime_utr_features = create_three_prime_UTR_feature()
    cds_features = create_CDS_feature()
    output = sys.argv[9]
    with open(output, 'a') as out:
        for iso in transcript_features:
            for i in transcript_features[iso]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso1 in exon_features:
            for i in exon_features[iso1]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso2 in start_codon_features:
            for i in start_codon_features[iso2]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso3 in stop_codon_features:
            for i in stop_codon_features[iso3]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso4 in five_prime_utr_features:
            for i in five_prime_utr_features[iso4]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso5 in three_prime_utr_features:
            for i in three_prime_utr_features[iso5]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso6 in cds_features:
            for i in cds_features[iso6]:
                final = "\t".join(i)
                out.write(final + "\n")

write()

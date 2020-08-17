#creating final gtf for isoseq data
#gtf format will have 9 columns = chr num, source, feature, start, end, score, strand, frame, attributes
#source will be PacBio; feature can be gene, isoform (transcript), exon, CDS, start_codon, stop_codon, five_prime_utr, three_prime_utr; frame is only applicable to protein coding genes in the CDS
#non protein coding genes will only have gene, transcript, and exon features
#protein coding genes can have all features (though some may be missing)
#attributes will include (when applicable): gene id, transcript id, gene source, transcript source, gene biotype, transcript biotype, exon number, gene name
#if a value is not available, this will show up as "." if only value in column or be left blank if in the attribute section
#this script handles non protein coding isoforms (only need transcript and exon features)
#to run script: python3 Create.Full.GTF.IsoformIssuesHandling.9.py <source as string; in this case all sources will be PacBio> <classification file for combined analyses> <full bed file> <full Ensembl GTF> <isoseq gtf with exon positions> <isoforms to examine; file contains just isoform ids> <output file>
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
def read_bed():
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


def pull_isoforms_to_analyze():
    isoforms_file = sys.argv[6]
    correct_isoforms = []
    with open(isoforms_file, 'r') as isoforms:
        for line in isoforms:
            new_line = line.strip("\n")
            correct_isoforms.append(new_line)
    return correct_isoforms


#create transcript feature
#this will include:
#chr num, source=PacBio, feature=transcript, start, end, score=".", strand, frame = ".", attributes = gene id, gene source, gene biotype, gene name if available, transcript id, transcript source, transcript biotype; isoform id
#difference between transcript and isoform id: transcript id will be either ENSGACT or novel and isoform id is PB.XXXX.X
def create_transcript_feature():
    class_dict = read_class()
    bed_dict = read_bed()
    correct_isoforms = pull_isoforms_to_analyze()
    ensembl_dict = read_ensembl_gtf_transcripts()
    source = sys.argv[1]
    transcript_feature_dict = {}
    for isoform in correct_isoforms:
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
    class_dict = read_class()
    ensembl_dict = read_ensembl_gtf_genes()
    exon_dict = read_isoseq_gtf()
    correct_isoforms = pull_isoforms_to_analyze()
    source = sys.argv[1]
    exon_feature_dict = {}
    for isoform in correct_isoforms:
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


#writing output to files
#Will write all features to 1 file that will be sorted later
def write():
    transcript_features = create_transcript_feature()
    exon_features = create_exon_feature()
    output = sys.argv[7]
    with open(output, 'a') as out:
        for iso1 in transcript_features:
            for i in transcript_features[iso1]:
                final = "\t".join(i)
                out.write(final + "\n")
        for iso2 in exon_features:
            for i in exon_features[iso2]:
                final = "\t".join(i)
                out.write(final + "\n")
write()

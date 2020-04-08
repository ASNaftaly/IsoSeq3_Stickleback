#creating filtered input for sqanti_qc.py
#have filtered classification file (with artifacts marked) and curated transcriptome (isoforms that are not labeled as isoforms)
#need gtf and abundance files for sqanti_qc.py
#gtf format:
#no headers
#chr# PacBio exon exon.start exon.end . strand . transcript_id "PB.X.X"; gene_id "PB.Y.Y";
#abundance format:
#several headers with #
#Pbid count_fl count_nfl count_nfl_amb norm_fl norm_nfl norm_nfl_amb
#can pull isoform ids (PB id) from curated transcriptome and remove these from gtf and abundance file then run sqanti_qc.py again
#1. pull curated transcriptome isoform ids, unfiltered gtf file, unfiltered abundance file
#2. sort gtf and abundance file to only include "true" isoforms
#3. write new gtf and abundance files
#to run script: python3 FilterSQANTI_input.py <curated file from Sqanti_filter.py> <unfiltered gtf file used as input for Sqanti_qc.py> <unfiltered abundance file used as input for Sqanti_qc.py> <output gtf filtered> <output abundance filtered>
#Author: Alice Naftaly, September 3, 2019

import sys

#pulls transcripts from curated transcripts file
#returns list with each value being a single isoform ID
def pull_filtered_isoforms():
    curated_transcripts_file = sys.argv[1]
    filtered_transcripts = []
    with open(curated_transcripts_file, 'r') as transcripts:
        for line in transcripts:
            new_line = line.strip("\n")
            filtered_transcripts.append(new_line)
    return filtered_transcripts


#pulls gtf file into a dictionary
#key = transcript id (PB.X.X), value = whole gtf line
def pull_unfiltered_gtf():
    gtf_file = sys.argv[2]
    gtf_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            raw_transcript_id = new_line[9].strip(";")
            transcript_id = raw_transcript_id.strip('\"')
            if transcript_id in gtf_dict:
                gtf_dict[transcript_id].append(line)
            elif transcript_id not in gtf_dict:
                gtf_dict.update({transcript_id:[line]})
    return gtf_dict

#pulls headers from abundance file to be added to final abundance file
#returns list with each value being 1 line of the headers (in order)
#should be 14 headers + 1 line for header for data
def pull_abundance_headers():
    abundance_file = sys.argv[3]
    abundance_headers = []
    with open(abundance_file, 'r') as abundance:
        for line in abundance:
            if line.startswith("#"):
                abundance_headers.append(line)
            elif line.startswith("pbid"):
                abundance_headers.append(line)
    return abundance_headers

#pulls abundance data
#returns dictionary with key = transcript Id and value = whole line
def pull_abundance():
    abundance_file = sys.argv[3]
    abundance_dict = {}
    with open(abundance_file, 'r') as abundance:
        for line in abundance:
            if line.startswith("#"):
                continue
            elif line.startswith("pbid"):
                continue
            else:
                new_line = line.split()
                transcript_id = new_line[0]
                if transcript_id in abundance_dict:
                    abundance_dict[transcript_id].append(line)
                elif transcript_id not in abundance_dict:
                    abundance_dict.update({transcript_id:[line]})
    return abundance_dict


#sort out artifact isoforms from gtf
#returns dictionary with key = transcript Id and value = all lines from gtf for a given transcript id
#only returns transcript IDs that are not artifacts
def sort_gtf():
    filtered_transcripts = pull_filtered_isoforms()
    unfiltered_gtf = pull_unfiltered_gtf()
    filtered_gtf = {}
    for id in filtered_transcripts:
        if id in unfiltered_gtf:
            true_isoform = unfiltered_gtf[id]
            if id in filtered_gtf:
                filtered_gtf[id].append(true_isoform)
            elif id not in filtered_gtf:
                filtered_gtf.update({id:true_isoform})
    return filtered_gtf

#sorts abundance file to remove artifacts
#returns dictionary with key = transcript ID and value = lines from abundance file for that transcript ID
def sort_abundance():
    filtered_transcripts = pull_filtered_isoforms()
    unfiltered_abundance = pull_abundance()
    filtered_abundance = {}
    for id in filtered_transcripts:
        if id in unfiltered_abundance:
            true_isoform = unfiltered_abundance[id]
            if id in filtered_abundance:
                filtered_abundance[id].append(true_isoform)
            elif id not in filtered_abundance:
                filtered_abundance.update({id:true_isoform})
    return filtered_abundance


#need to order gtf and abundance file
#returns list of id numbers sorted in numerical order to write the gtf file in numerical order
#can use the same list for both GTF and abudance file
def order_ids():
    gtf_lines = sort_gtf()
    all_ids = []
    final_id_list = []
    for key in gtf_lines:
        id_number = key.strip("PB.")
        split_ids = id_number.split(".")
        ids = [split_ids[0], split_ids[1]]
        all_ids.append(ids)
    sorted_ids = sorted(all_ids,key=lambda x: (int(x[0]), int(x[1])))
    for value in sorted_ids:
        final_id = value[0] + "." + value[1]
        final_id_list.append(final_id)
    return final_id_list


#writing filtered gtf file and abundance file
#writes gtf file in appropriate order
def write_gtf():
    gtf_lines = sort_gtf()
    ordered_ids = order_ids()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for id in ordered_ids:
            full_id = "PB." + str(id)
            gtf_id = gtf_lines[full_id]
            for value in gtf_id:
                out.write(value)

#writes abundance file in appropriate order with headers
def write_abundance():
    abundance_headers = pull_abundance_headers()
    abundance_lines = sort_abundance()
    ordered_ids = order_ids()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for value in abundance_headers:
            out.write(value)
        for id in ordered_ids:
            full_id = "PB." + str(id)
            abundance_id = abundance_lines[full_id]
            for val in abundance_id:
                out.write(val)

#calling all functions:
def call():
    gtf= write_gtf()
    abundance = write_abundance()


call()

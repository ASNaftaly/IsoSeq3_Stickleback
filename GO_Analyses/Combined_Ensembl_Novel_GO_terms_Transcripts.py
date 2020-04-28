#combined annotated ensembl gene go terms and novel gene GO terms
#this script will be used by Randomize_GO_terms_isoseq.py to create random distribution
#will also output just novel isoforms and ensembl transcripts as separate files (all novel isoforms and all ensembl transcripts) in the same format as the combined file to run specific analyses
#to run script: python3 Combined_Ensembl_Novel_GO_terms_Transcripts.py <ensembl table as csv> <novel isoforms go terms file for combined sexes analysis from Pull_GO_terms_Interproscan_Output.py> <output file;combined novel and ensembl format: Transcript.ID \t List of GO terms> <output file; ensembl transcripts> <output file; novel isoforms>
#Author: Alice Naftaly, April 2020

import sys

#first pull genes and go terms from ensembl table
#only need transcript id and go terms
#returns dictionary with key = transcript id and value == list of go terms
def read_ensembl_table():
    ensembl_table_file = sys.argv[1]
    go_term_dict = {}
    final_go_terms = {}
    with open(ensembl_table_file, 'r') as ensembl_table:
        for line in ensembl_table:
            if line.startswith("ENSGACG"):
                new_line = line.split(",")
                transcript_id = new_line[1]
                go_term = new_line[2]
                if transcript_id in go_term_dict:
                    go_term_dict[transcript_id].append(go_term)
                elif transcript_id not in go_term_dict:
                    go_term_dict.update({transcript_id:[go_term]})
    #remove genes with no go terms or empty go term values
    for key in go_term_dict:
        single_key = go_term_dict[key]
        list_of_go_terms = []
        for value in single_key:
            if value.startswith("GO:"):
                list_of_go_terms.append(value)
        if len(list_of_go_terms) > 0:
            final_go_terms.update({key:list_of_go_terms})
    return final_go_terms


#need to add the novel genes detected to the ensembl annotated genes
#read GO terms for novel genes
#returns dictionary with key == isoform id and value == [list of go terms]
def read_novel_GO_terms():
    go_terms_file = sys.argv[2]
    go_term_dict = {}
    with open(go_terms_file, 'r') as go_terms:
        for line in go_terms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                single_go_term = new_line[1]
                if isoform_id in go_term_dict:
                    go_term_dict[isoform_id].append(single_go_term)
                elif isoform_id not in go_term_dict:
                    go_term_dict.update({isoform_id:[single_go_term]})
    return go_term_dict


#need to combine go terms from ensembl annotations and novel genes for combined sexes
#returns one dictionary with key == gene and value == list of go terms for that gene
def combine_ensembl_novel_transcripts():
    ensembl_transcripts = read_ensembl_table()
    novel_transcripts = read_novel_GO_terms()
    combined_final_dict = {}
    combined_final_dict.update(ensembl_transcripts)
    combined_final_dict.update(novel_transcripts)
    return combined_final_dict

#write combined dictionary to output file
#final format should be gene id \t list of go terms as list separated by commas
def write_combined():
    final_dict = combine_ensembl_novel_transcripts()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for transcript in final_dict:
            single_transcript = final_dict[transcript]
            go_term_list = ",".join(single_transcript)
            final = "%s\t%s\n" % (str(transcript), go_term_list)
            out.write(final)

def write_ensembl():
    final_dict = read_ensembl_table()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for transcript in final_dict:
            single_transcript = final_dict[transcript]
            go_term_list = ",".join(single_transcript)
            final = "%s\t%s\n" % (str(transcript), go_term_list)
            out.write(final)

def write_novel():
    final_dict = read_novel_GO_terms()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for transcript in final_dict:
            single_transcript = final_dict[transcript]
            go_term_list = ",".join(single_transcript)
            final = "%s\t%s\n" % (str(transcript), go_term_list)
            out.write(final)

#call all functions()
def call():
    combined = write_combined()
    ensembl = write_ensembl()
    novel = write_novel()

call()

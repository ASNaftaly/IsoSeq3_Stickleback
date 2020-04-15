#breakdown of annotated vs novel genes for ncRNAs
#will examine short and long ncRNAs separately
#to run script: python3 Annotated_vs_Novel_ncRNAs.py <classification file for noncoding RNAs> <output gene ids for annotated genes> <output gene ids for novel genes> <output isoform ids for annotated transcripts> <outptu isoform ids for novel transcripts> >> <summary.log>
#Author: Alice Naftaly, April 2020

import sys

#read classification file
#returns dictionary with key == isoform and value == [gene id, transcript id]
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                dict_value = [gene_id, transcript_id]
                class_dict.update({isoform:dict_value})
    return class_dict

#counting annotated genes vs novel genes:
#returns lists the annotated genes and novel genes
def classify_genes():
    class_dict = read_class()
    annotated_genes = []
    novel_genes = []
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        gene_id = single_isoform[0]
        if gene_id.startswith("ENSGACG"):
            annotated_genes.append(gene_id)
        elif gene_id.startswith("novel"):
            novel_genes.append(gene_id)
    set_annotated_genes = set(annotated_genes)
    set_novel_genes = set(novel_genes)
    return set_annotated_genes, set_novel_genes


#counting annotated transcripts vs novel transcripts
#returns lists of isoform ids for annotated transcripts and novel transcripts
def classify_transcripts():
    class_dict = read_class()
    annotated_transcripts = []
    novel_transcripts = []
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        transcript_id = single_isoform[1]
        if transcript_id.startswith("ENSGACT"):
            annotated_transcripts.append(isoform)
        elif transcript_id.startswith("novel"):
            novel_transcripts.append(isoform)
    return annotated_transcripts, novel_transcripts


#writing output files and summary statements

def write_annotated_genes():
    annotated_genes, novel_genes = classify_genes()
    output = sys.argv[2]
    num_annotated_genes = len(annotated_genes)
    summary = "Total number of annotated ncRNA genes: %s" % str(num_annotated_genes)
    print(summary)
    with open(output ,'a') as out:
        for isoform in annotated_genes:
            final = "%s\n" % str(isoform)
            out.write(final)

def write_novel_genes():
    annotated_genes, novel_genes = classify_genes()
    output = sys.argv[3]
    num_novel_genes = len(novel_genes)
    summary = "Total number of novel ncRNA genes: %s" % str(num_novel_genes)
    print(summary)
    with open(output ,'a') as out:
        for isoform in novel_genes:
            final = "%s\n" % str(isoform)
            out.write(final)


def write_annotated_transcripts():
    annotated_transcripts, novel_transcripts = classify_transcripts()
    output = sys.argv[4]
    num_annotated_transcripts = len(annotated_transcripts)
    summary = "Total number of annotated ncRNA transcripts: %s" % str(num_annotated_transcripts)
    print(summary)
    with open(output ,'a') as out:
        for isoform in annotated_transcripts:
            final = "%s\n" % str(isoform)
            out.write(final)

def write_novel_transcripts():
    annotated_transcripts, novel_transcripts = classify_transcripts()
    output = sys.argv[5]
    num_novel_transcripts = len(novel_transcripts)
    summary = "Total number of novel ncRNA transcripts: %s" % str(num_novel_transcripts)
    print(summary)
    with open(output ,'a') as out:
        for isoform in novel_transcripts:
            final = "%s\n" % str(isoform)
            out.write(final)

#call all functions:
def call():
    annotated_genes = write_annotated_genes()
    novel_genes = write_novel_genes()
    annotated_transcripts = write_annotated_transcripts()
    novel_transcripts = write_novel_transcripts()
call()

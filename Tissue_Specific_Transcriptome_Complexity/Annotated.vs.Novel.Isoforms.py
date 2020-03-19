#Pull annotated vs novel isoforms
#First have to look at genes, to pull annotated genes and novel genes separately
#Don't have to look at novel genes to count novel transcripts as they all will be novel
#The gene ids have already been converted to the combined sexes gene ids
#to run script: python3 Annotated.vs.Novel.Isoforms.py <individual tissue exon counts which has Isoform ids of correctly matched isoforms> <combined sexes classification file> <output, annotated genes Isoform.ID\tGene.ID> <output, novel genes Isoform.ID\tGene.ID> <output, annotated genes, annotated transcripts Isoform.ID\tGene.ID\tTranscript.ID> <output, annotated genes, novel transcripts Isoform.ID\tGene.ID\tTranscript.ID>

import sys

#read individual exon file
#returns a list of isoform ids
def read_isoform_ids():
    exon_file = sys.argv[1]
    isoform_ids = []
    with open(exon_file, 'r') as exons:
        for line in exons:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_ids.append(new_line[0])
    return isoform_ids

#read combined classification file
#need the gene id and the transcript id
#returns dictionary with key = isoform id and value == [gene id, transcript id]
def read_combined_class():
    class_file = sys.argv[2]
    isoform_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                final = [gene_id, transcript_id]
                isoform_dict.update({isoform_id:final})
    return isoform_dict

#pull genes and transcripts to look at
#there are some isoforms that have mismatches to the combined sexes isoforms (these are removed in the exon file)
#returns dictionary with key = isoform id and value == [gene id, transcript id]
def remove_mismatches():
    isoforms = read_isoform_ids()
    class_dict = read_combined_class()
    final_dict = {}
    for isoform in isoforms:
        single_class = class_dict[isoform]
        final_dict.update({isoform:single_class})
    return final_dict

#look at annotated genes
#separates out annotated genes and determine how many of these have annotated transcripts vs novel transcripts
#all novel genes are novel Transcripts
#prints out summary statements
def transcript_summary():
    final_dict = remove_mismatches()
    annotated_transcripts = []
    annotated_novel_transcripts = []
    novel_transcripts = []
    for isoform in final_dict:
        single_isoform = final_dict[isoform]
        gene_id = single_isoform[0]
        transcript_id = single_isoform[1]
        if gene_id.startswith("EN"):
            if transcript_id.startswith("EN"):
                annotated_transcripts.append(transcript_id)
            elif transcript_id.startswith("novel"):
                annotated_novel_transcripts.append(transcript_id)
        elif gene_id.startswith("novel"):
            novel_transcripts.append(transcript_id)
    summary_annotated_genes_annotated_transcripts = "Number of Annotated Genes, Annotated Transcripts: %s\n" % str(len(annotated_transcripts))
    summary_annotated_genes_novel_transcripts = "Number of Annotated Genes, Novel Transcripts: %s\n" % str(len(annotated_novel_transcripts))
    summary_novel_genes_novel_transcripts = "Number of Novel Genes, Novel Transcripts: %s\n" % str(len(novel_transcripts))
    print(summary_annotated_genes_annotated_transcripts)
    print(summary_annotated_genes_novel_transcripts)
    print(summary_novel_genes_novel_transcripts)

#now write annotated and novel genes and annotated and novel transcripts to files
#will write isoform id and gene id/transcript id
#writes Isoform id \t gene Id
def write_annotated_genes():
    final_dict = remove_mismatches()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Isoform.ID\tGene.ID\n"
        out.write(header)
        for isoform in final_dict:
            single_isoform = final_dict[isoform]
            gene_id = single_isoform[0]
            if gene_id.startswith("EN"):
                final = "%s\t%s\n" % (str(isoform), str(gene_id))

def write_novel_genes():
    final_dict = remove_mismatches()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Isoform.ID\tGene.ID\n"
        out.write(header)
        for isoform in final_dict:
            single_isoform = final_dict[isoform]
            gene_id = single_isoform[0]
            if gene_id.startswith("novel"):
                final = "%s\t%s\n" % (str(isoform), str(gene_id))

def write_annotated_genes_annotated_transcripts():
    final_dict = remove_mismatches()
    output = sys.argv[5]
    with open(output, 'a') as out:
        header = "Isoform.ID\tGene.ID\tTranscript.ID\n"
        out.write(header)
        for isoform in final_dict:
            single_isoform = final_dict[isoform]
            gene_id = single_isoform[0]
            transcript_id = single_isoform[1]
            if gene_id.startswith("EN"):
                if transcript_id.startswith("EN"):
                    final = "%s\t%s\t%s\n" % (str(isoform), str(gene_id),str(transcript_id))

def write_annotated_genes_novel_transcripts():
    final_dict = remove_mismatches()
    output = sys.argv[6]
    with open(output, 'a') as out:
        header = "Isoform.ID\tGene.ID\tTranscript.ID\n"
        out.write(header)
        for isoform in final_dict:
            single_isoform = final_dict[isoform]
            gene_id = single_isoform[0]
            transcript_id = single_isoform[1]
            if gene_id.startswith("EN"):
                if transcript_id.startswith("novel"):
                    final = "%s\t%s\t%s\n" % (str(isoform), str(gene_id),str(transcript_id))

#call all functions:
def call():
    summary = transcript_summary()
    annotated_genes = write_annotated_genes()
    novel_genes = write_annotated_genes()
    annotated_transcripts = write_annotated_genes_annotated_transcripts()
    novel_transcripts = write_annotated_genes_novel_transcripts()

call()

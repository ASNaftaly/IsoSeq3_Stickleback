#examining differences among novel transcripts
#need to read in novel transcripts from classification and gtf files and ensembl genes to compare differences
#to run script: python3 NovelIsoform_Comparison_to_Ensembl.py

import sys


#read in novel transcripts classification file
#returns dictionary with key == isoform id and value == [gene_id, splice type, coding]
def read_novel_class():
    novel_file = sys.argv[1]
    novel_dict = {}
    with open(novel_file, 'r') as novel:
        for line in novel:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                splice_type = new_line[5]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                coding = new_line[27]
                if transcript_id == "novel" and gene_id.startswith("ENSGACG"):
                    final = [gene_id, splice_type, coding]
                    novel_dict.update({isoform_id:final})
    return novel_dict


#read novel gtf file
#returns dictionary with key == isoform id and value == [exon_start, exon_end] for each exon
def read_gtf():
    novel_gtf = sys.argv[2]
    novel_gtf_dict = {}
    with open(novel_gtf, 'r') as novel:
        for line in novel:
            new_line = line.split("\t")
            exon_pos_1 = new_line[3]
            exon_pos_2 = new_line[4]
            strand = new_line[6]
            transcript_info = new_line[8].split(";")
            isoform_id_full = transcript_info[0].split(" ")
            isoform_id = isoform_id_full[1].strip("\"")
            if strand == "+":
                exon_start = exon_pos_1
                exon_end = exon_pos_2
            elif strand == "-":
                exon_start = exon_pos_2
                exon_end = exon_pos_1
            final = [exon_start, exon_end]
            if isoform_id in novel_gtf_dict:
                novel_gtf_dict[isoform_id].append(final)
            elif isoform_id not in novel_gtf_dict:
                novel_gtf_dict.update({isoform_id:[final]})

    return novel_gtf_dict


#filter novel gtf to only novel transcripts
#returns a final dictionary with only novel transcripts
def filter_novel_gtf():
    full_gtf = read_gtf()
    novel_transcripts = read_novel_class()
    novel_isoforms = novel_transcripts.keys()
    final_gtf_dict = {}
    for isoform in novel_isoforms:
        if isoform in full_gtf:
            final_gtf_dict.update({isoform:full_gtf[isoform]})
    return final_gtf_dict


#read in ensembl gtf file
#will want gene id and exons
#returns dictionary with key == gene id and value == [transcript id, [exons]]
def read_ensembl_gtf():
    ensembl_file = sys.argv[3]
    transcript_dict = {}
    exon_dict = {}
    final_dict = {}
    with open(ensembl_file, 'r') as ensembl:
        for line in ensembl:
            new_line = line.split("\t")
            region = new_line[2]
            if region == "transcript":
                full_transcript = new_line[8].split(";")
                gene_id_full = full_transcript[0].split(" ")
                gene_id = gene_id_full[1].strip("\"")
                transcript_id_full = full_transcript[2].split(" ")
                transcript_id = transcript_id_full[2].strip("\"")
                transcript_dict.update({transcript_id:gene_id})
            elif region == "exon":
                strand = new_line[6]
                exon_pos_1 = new_line[3]
                exon_pos_2 = new_line[4]
                full_transcript = new_line[8].split(";")
                transcript_id_full = full_transcript[2].split(" ")
                transcript_id = transcript_id_full[2].strip("\"")
                if strand == "_":
                    exon_start = exon_pos_1
                    exon_end = exon_pos_2
                elif strand == "-":
                    exon_start = exon_pos_2
                    exon_end = exon_pos_1
                final = [exon_start, exon_end]
                if transcript_id in exon_dict:
                    exon_dict[transcript_id].append(final)
                elif transcript_id not in exon_dict:
                    exon_dict.update({transcript_id:[final]})
        for transcript in transcript_dict:
            associated_gene = transcript_dict[transcript]
            exons = exon_dict[transcript]
            final_value = [transcript, exons]
            if associated_gene in final_dict:
                final_dict[associated_gene].append(final_value)
            elif associated_gene not in final_dict:
                final_dict.update({associated_gene:[final_value]})
    return final_dict



#comparing novel transcripts to ensembl genes
def compare():
    novel_transcripts = read_novel_class()
    novel_gtf = filter_novel_gtf()
    ensembl_dict = read_ensembl_gtf()
    for isoform in novel_transcripts:
        single_isoform_info = novel_transcripts[isoform]
        gene_id = single_isoform_info[0]
        single_isoform_exons = novel_gtf[isoform]
        isoseq_number_exons = len(single_isoform_exons)
        if gene_id in ensembl_dict:
            single_ensembl_gene = ensembl_dict[gene_id]
            for transcript in single_ensembl_gene:
                ensembl_exons = transcript[1]
                ensembl_num_exons = len(ensembl_exons)
                #if the number of exons is the same between isoseq and ensembl:
                if isoseq_number_exons == ensembl_num_exons:
                    print(isoform)
                    print(single_isoform_exons)
                    print(ensembl_exons)
compare()

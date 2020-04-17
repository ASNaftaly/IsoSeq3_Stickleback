#determine if ncRNAs are intergenic (not overlapping any Ensembl or Isoseq gene/transcript); intronic (completely within the intron of an Ensembl gene or an Isoseq gene), or antisense (in opposite direction/strand of a protein coding gene)
#will read in ensembl GTF (with tss positions) and isoseq classification (combined sexes) and collapsed gff from cDNAcupcake files for gene/transcript positions; gtf files for exon positions; directionality will be determined by positions in classification and gtf files
#to run script: python3 Locating_ncRNAS_GenomeWide_Search.py <ensembl gtf full> <isoseq classification file with duplicates removed> <isoseq gff file from cDNAcupcake> <isoseq gtf from SQANTI> <noncoding RNAs classification file> <output intergenic ncRNAs isoform ids> <output intronic ncRNAs isoform ids> <output antisense ncRNAs isoform ids> <output remaining ncRNAs that do not fit the previous three categories isoform ids> >> <log file>
#author: Alice Naftaly, April 2020


import sys

#read in gtf file from Ensembl with exons, tss, separately
#will want to read in transcript for start and end positions with strand
#will want to read in exons with start, end positions, strand, and exon number
#returns dictionary with key == transcript id and value == [gene_id, transcript_id, chr_num, transcript_start_pos, transcript_end_pos, strand]
def read_ensembl_gtf_transcript_positions():
    ensembl_gtf = sys.argv[1]
    ensembl_transcript_dict = {}
    with open(ensembl_gtf, 'r') as ensembl:
        for line in ensembl:
            new_line = line.split("\t")
            chr_num = new_line[0]
            type = new_line[2]
            type_start_pos = new_line[3]
            type_end_pos = new_line[4]
            strand = new_line[6]
            gene_info = new_line[8].split(";")
            for value in gene_info:
                if value.startswith("gene_id"):
                    gene_id_full = value.split(" ")
                    gene_id = gene_id_full[1].strip("\"")
                elif value.startswith(" transcript_id"):
                    transcript_id_full = value.split(" ")
                    transcript_id = transcript_id_full[2].strip("\"")
            if type == "transcript":
                if strand == "+":
                    start_pos = new_line[3]
                    end_pos = new_line[4]
                elif strand == "-":
                    start_pos = new_line[4]
                    end_pos = new_line[3]
                dict_value = [gene_id, transcript_id, chr_num, start_pos, end_pos, strand]
                ensembl_transcript_dict.update({transcript_id:dict_value})
    return ensembl_transcript_dict


#returns dictionary with key == transcript id and value == for all lexons *[gene_id, transcript_id, chr_num, exon_start_pos, exon_end_pos, strand, exon number]
def read_ensembl_gtf_exon_positions():
    ensembl_gtf = sys.argv[1]
    ensembl_exons_dict = {}
    with open(ensembl_gtf, 'r') as ensembl:
        for line in ensembl:
            new_line = line.split("\t")
            chr_num = new_line[0]
            type = new_line[2]
            strand = new_line[6]
            gene_info = new_line[8].split(";")
            for value in gene_info:
                if value.startswith("gene_id"):
                    gene_id_full = value.split(" ")
                    gene_id = gene_id_full[1].strip("\"")
                elif value.startswith(" transcript_id"):
                    transcript_id_full = value.split(" ")
                    transcript_id = transcript_id_full[2].strip("\"")
                elif value.startswith(" exon_number"):
                    exon_number_full = value.split(" ")
                    exon_number = exon_number_full[2].strip("\"")
            if type == "exon":
                if strand == "+":
                    start_pos = new_line[3]
                    end_pos = new_line[4]
                elif strand == "-":
                    start_pos = new_line[4]
                    end_pos = new_line[3]
                dict_value = [gene_id, transcript_id, chr_num, start_pos, end_pos, strand, exon_number]
                if transcript_id in ensembl_exons_dict:
                    ensembl_exons_dict[transcript_id].append(dict_value)
                elif transcript_id not in ensembl_exons_dict:
                    ensembl_exons_dict.update({transcript_id:[dict_value]})
    return ensembl_exons_dict


#the collapsed gff file created from cDNAcupcake for the collapsing identical isoforms has the start and end positions for all of the transcripts.
#will have to read in the classification file to get the isoforms after the SQANTI filtering and removal of duplicates then read in gff and pull transcript positions for all transcripts (including novel transcripts)

#just need isoform ids from this file
#returns a list of isoform ids
def read_nodups_classification():
    class_file = sys.argv[2]
    isoform_ids_list = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                isoform_ids_list.append(isoform)
    return isoform_ids_list

#read in gff file
#returns dictionary with key == isoform id and value == [chr_num, isoform_id, start_pos, end_pos, strand]
def read_isoseq_collapsed_gff():
    gff_file = sys.argv[3]
    gff_dict = {}
    with open(gff_file, 'r') as gff:
        for line in gff:
            new_line = line.split()
            chr_num = new_line[0]
            type = new_line[2]
            strand = new_line[6]
            isoform_id_full = new_line[11].strip(";")
            isoform_id = isoform_id_full.strip("\"")
            if type == "transcript":
                if strand == "+":
                    start_pos = new_line[3]
                    end_pos = new_line[4]
                elif strand == "-":
                    start_pos = new_line[4]
                    end_pos = new_line[3]
                dict_value = [chr_num, isoform_id, start_pos, end_pos, strand]
                gff_dict.update({isoform_id:dict_value})
    return gff_dict

#removing isoforms from the gff dictionary as the gff has all of the isoforms prior to SQANTI filtering and duplicate removal
#returns dictionary with key = isoform id and value == [chr_num, isoform id, start pos, end pos, strand]
def filter_isoseq():
    filtered_isoforms = read_nodups_classification()
    all_isoforms = read_isoseq_collapsed_gff()
    final_isoseq_transcript_dict = {}
    for isoform in filtered_isoforms:
        if isoform in all_isoforms:
            single_isoform = all_isoforms[isoform]
            final_isoseq_transcript_dict.update({isoform:single_isoform})
    return final_isoseq_transcript_dict


#now reading in exons for isoseq isoforms (can use the gtf file created from SQANTI)
#takes into account strand to flip exons around
#returns dictionary with key == isoform and value == [chr_num, isoform id, strand, exon number, exon start, exon end]
def read_isoseq_gtf():
    gtf_file = sys.argv[4]
    isoseq_exon_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            chr_num = new_line[0]
            strand = new_line[6]
            isoform_id_full = new_line[9].strip(";")
            isoform_id = isoform_id_full.strip("\"")
            if strand == "+":
                exon_start = new_line[3]
                exon_end = new_line[4]
            elif strand == "-":
                exon_start = new_line[4]
                exon_end = new_line[3]
            if isoform_id in isoseq_exon_dict:
                 exon_number += 1
                 dict_value = [chr_num, isoform_id, strand, str(exon_number), exon_start, exon_end]
                 isoseq_exon_dict[isoform_id].append(dict_value)
            elif isoform_id not in isoseq_exon_dict:
                exon_number = 1
                dict_value = [chr_num, isoform_id, strand, str(exon_number), exon_start, exon_end]
                isoseq_exon_dict.update({isoform_id:[dict_value]})
    return isoseq_exon_dict

#need to read in ncRNAs isoform ids and remove these from the isoseq dictionaries and ensembl as needed.
#returns dictionary with key == isoform id and value == transcript id (needed to remove genes from ensembl data sets)
def read_ncRNAs():
    ncRNAs_file = sys.argv[5]
    ncRNAs_isoforms = {}
    with open(ncRNAs_file, 'r') as ncRNAs:
        for line in ncRNAs:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                transcript_id = new_line[7]
                ncRNAs_isoforms.update({isoform:transcript_id})
    return ncRNAs_isoforms

#need to get the ncRNAs positions (i.e. transcript)
def get_ncRNAs_transcript_positions():
    ncRNAs = read_ncRNAs()
    isoseq_positions = filter_isoseq()
    final_ncRNAs = {}
    for isoform in ncRNAs:
        single_position = isoseq_positions[isoform]
        final_ncRNAs.update({isoform:single_position})
    return final_ncRNAs

#remove ncRNAs from isoseq and ensembl dictionaries
def remove_ncRNAs_ensembl_transcripts():
    ncRNAs_isoforms = read_ncRNAs()
    ensembl_transcripts = read_ensembl_gtf_transcript_positions()
    for isoform in ncRNAs_isoforms:
        transcript = ncRNAs_isoforms[isoform]
        if transcript in ensembl_transcripts:
            del ensembl_transcripts[transcript]
    return ensembl_transcripts

def remove_ncRNAs_isoseq_transcripts():
    ncRNAs_isoforms = read_ncRNAs()
    isoseq_transcripts = filter_isoseq()
    for isoform in ncRNAs_isoforms:
        if isoform in isoseq_transcripts:
            del isoseq_transcripts[isoform]
    return isoseq_transcripts

def remove_ncRNAs_ensembl_exons():
    ncRNAs_isoforms = read_ncRNAs()
    ensembl_exons = read_ensembl_gtf_exon_positions()
    for isoform in ncRNAs_isoforms:
        transcript = ncRNAs_isoforms[isoform]
        if transcript in ensembl_exons:
            del ensembl_exons[transcript]
    return ensembl_exons

def remove_ncRNAs_isoseq_exons():
    ncRNAs_isoforms = read_ncRNAs()
    isoseq_exons = read_isoseq_gtf()
    for isoform in ncRNAs_isoforms:
        if isoform in isoseq_exons:
            del isoseq_exons[isoform]
    return isoseq_exons


########Analysis portion of script

#how many ncRNAs are intergenic (do not overlap with any genes)

#ensembl examined first
def compare_ncRNAs_ensembl_transcripts():
    ensembl_transcripts = remove_ncRNAs_ensembl_transcripts()
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    intergenic_ensembl_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        counts = 0
        #iterate through ensembl transcripts
        for transcript in ensembl_transcripts:
            single_ensembl_transcript = ensembl_transcripts[transcript]
            ensembl_gene_id = single_ensembl_transcript[0]
            ensembl_chr_num = single_ensembl_transcript[2]
            ensembl_start_pos = int(single_ensembl_transcript[3])
            ensembl_end_pos = int(single_ensembl_transcript[4])
            ensembl_strand = single_ensembl_transcript[5]
            #if ncRNA overlaps at all with the ensembl transcript counts will increase by 1
            if ncRNA_chr_num == ensembl_chr_num and ncRNA_strand == ensembl_strand:
                if ensembl_strand == "+":
                    #if the ncRNA is completely within an ensembl transcript
                    if ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the beginning of the transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the end of the transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                elif ensembl_strand == "-":
                    #if ncRNA is completely within an ensembl transcripts
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
        if counts == 0:
            intergenic_ensembl_ncRNAs.update({isoform:single_ncRNA_isoform})
    return intergenic_ensembl_ncRNAs

#examining isoseq transcripts
def compare_ncRNAs_isoseq_transcripts():
    isoseq_transcripts = remove_ncRNAs_isoseq_transcripts()
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    intergenic_isoseq_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        counts = 0
        #iterate through ensembl transcripts
        for transcript in isoseq_transcripts:
            single_isoseq_transcript = isoseq_transcripts[transcript]
            isoseq_gene_id = single_isoseq_transcript[0]
            isoseq_chr_num = single_isoseq_transcript[1]
            isoseq_start_pos = int(single_isoseq_transcript[2])
            isoseq_end_pos = int(single_isoseq_transcript[3])
            isoseq_strand = single_isoseq_transcript[4]
            #if ncRNA overlaps at all with the isoseq transcript counts will increase by 1
            if ncRNA_chr_num == isoseq_chr_num and ncRNA_strand == isoseq_strand:
                if isoseq_strand == "+":
                    #if the ncRNA is completely within an isoseq transcript
                    if ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the beginning of the transcript
                    elif ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the end of the transcript
                    elif ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
                elif isoseq_strand == "-":
                    #if ncRNA is completely within an ensembl transcripts
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
        if counts == 0:
            intergenic_isoseq_ncRNAs.update({isoform:single_ncRNA_isoform})
    return intergenic_isoseq_ncRNAs


#pulling the intergenic ncRNAs
#if a ncRNA overlaps and ensembl or isoseq transcript, the transcript will not be in the list for the other transcriptome
#need to only keep the ncRNAs which do not overlap an transcript in both the ensembl and the isoseq database
#returns dictionary with ncRNAs that do not overlap with ensembl or isoseq transcripts
def combine_intergenic_ncRNAs():
    ensembl_intergenic_ncRNAs = compare_ncRNAs_ensembl_transcripts()
    isoseq_intergenic_ncRNAs = compare_ncRNAs_isoseq_transcripts()
    final_intergenic_ncRNAs = {}
    print("Total number of ncRNAs that do not overlap with ensembl transcripts")
    print(len(ensembl_intergenic_ncRNAs))
    print("Total number of ncRNAs that do not overlap with isoseq transcripts")
    print(len(isoseq_intergenic_ncRNAs))
    for isoform in ensembl_intergenic_ncRNAs:
        if isoform in isoseq_intergenic_ncRNAs:
            final_intergenic_ncRNAs.update({isoform:ensembl_intergenic_ncRNAs[isoform]})
    print("\n")
    return final_intergenic_ncRNAs


#how many ncRNAs are fully located in an intron?

#examining ensembl exons
#returns dictionary with key == isoform id and value == [chr num, isoform id, start pos, end pos, strand] *for ncRNA
def compare_ncRNAs_ensembl_exons():
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    ensembl_transcripts = remove_ncRNAs_ensembl_transcripts()
    ensembl_exons = remove_ncRNAs_ensembl_exons()
    ncRNAs_that_overlap_ensembl_transcripts = {}
    final_intronic_ensembl_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        for transcript in ensembl_transcripts:
            single_ensembl_transcript = ensembl_transcripts[transcript]
            ensembl_chr_num = single_ensembl_transcript[2]
            ensembl_start_pos = int(single_ensembl_transcript[3])
            ensembl_end_pos = int(single_ensembl_transcript[4])
            ensembl_strand = single_ensembl_transcript[5]
            #make sure chromosome number and strand are the same
            counts = 0
            if ncRNA_chr_num == ensembl_chr_num and ncRNA_strand == ensembl_strand:
                if ensembl_strand == "+":
                    #if the ncRNA is completely within an ensembl transcript
                    if ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the beginning of the transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the end of the transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                elif ensembl_strand == "-":
                    #if ncRNA is completely within an ensembl transcripts
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
            if counts > 0:
                ncRNAs_that_overlap_ensembl_transcripts.update({transcript:single_ncRNA_isoform})
    for transcript_overlap in ncRNAs_that_overlap_ensembl_transcripts:
        single_isoform = ncRNAs_that_overlap_ensembl_transcripts[transcript_overlap]
        ncRNA_start = int(single_isoform[2])
        ncRNA_end = int(single_isoform[3])
        single_ensembl_exons = ensembl_exons[transcript_overlap]
        #ncRNAs can only be in an intron if there is more than 1 exon
        exon_counts = 0
        if len(single_ensembl_exons) > 1:
            x = 1
            while x < len(single_ensembl_exons)-1:
                single_exon = single_ensembl_exons[x]
                single_exon_start = int(single_exon[3])
                single_exon_end = int(single_exon[4])
                single_strand = single_exon[5]
                if single_strand == "+":
                    #if the ncRNA is in exon
                    if ncRNA_start <= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end <= single_exon_end:
                        exon_counts ++ 1
                    #if the ncRNA overlaps the beginning of the exon
                    elif ncRNA_start >= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end <= single_exon_end:
                        exon_counts += 1
                    #if the ncRNA overlaps the end of the exon
                    elif ncRNA_start >= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end >= single_exon_end:
                        exon_counts += 1
                elif single_strand == "-":
                    #if ncRNA is completely within an ensembl transcripts
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                x += 1
            if counts == 0:
                #just need to check if ncRNA overlaps first or last exon at all
                last_counts = 0
                first_exon = single_ensembl_exons[0]
                last_exon = single_ensembl_exons[x]
                first_exon_start = int(first_exon[3])
                first_exon_end = int(first_exon[4])
                last_exon_start = int(last_exon[3])
                last_exon_end = int(last_exon[4])
                #checking for overlap with first_exon
                #+ strand
                if first_exon_start <= ncRNA_start and first_exon_end <= ncRNA_end and first_exon_start <= ncRNA_end and first_exon_end <= ncRNA_start:
                    last_counts += 1
                #- strand
                elif first_exon_start >= ncRNA_start and first_exon_end >= ncRNA_end and first_exon_start >= ncRNA_end and first_exon_end >= ncRNA_start:
                    last_counts += 1
                #checking for overlap with teh last exon
                #+ strand
                if last_exon_start >= ncRNA_start and last_exon_end >= ncRNA_end and last_exon_start >= ncRNA_end and last_exon_end >= ncRNA_start:
                    last_counts += 1
                #- strand
                elif last_exon_start <= ncRNA_start and last_exon_end <= ncRNA_end and last_exon_start <= ncRNA_end and last_exon_end <= ncRNA_start:
                    last_counts += 1
                if last_counts == 2:
                    isoform_id = single_isoform[1]
                    final_intronic_ensembl_ncRNAs.update({isoform_id:single_isoform})
    return final_intronic_ensembl_ncRNAs

#examining isoseq exons
def compare_ncRNAs_isoseq_exons():
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    isoseq_transcripts = remove_ncRNAs_isoseq_transcripts()
    isoseq_exons = remove_ncRNAs_isoseq_exons()
    ncRNAs_that_overlap_isoseq_transcripts = {}
    final_intronic_isoseq_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        for transcript in isoseq_transcripts:
            single_isoseq_transcript = isoseq_transcripts[transcript]
            isoseq_chr_num = single_isoseq_transcript[1]
            isoseq_start_pos = int(single_isoseq_transcript[2])
            isoseq_end_pos = int(single_isoseq_transcript[3])
            isoseq_strand = single_isoseq_transcript[4]
            #make sure chromosome number and strand are the same
            counts = 0
            if ncRNA_chr_num == isoseq_chr_num and ncRNA_strand == isoseq_strand:
                if isoseq_strand == "+":
                    #if the ncRNA is completely within an isoseq transcript
                    if ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the beginning of the transcript
                    elif ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    #if the ncRNA overlaps the end of the transcript
                    elif ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
                elif isoseq_strand == "-":
                    #if ncRNA is completely within an isoseq transcripts
                    if ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                            counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
            if counts > 0:
                ncRNAs_that_overlap_isoseq_transcripts.update({isoform:single_ncRNA_isoform})
    for isoform_overlap in ncRNAs_that_overlap_isoseq_transcripts:
        single_isoform = ncRNAs_that_overlap_isoseq_transcripts[isoform_overlap]
        ncRNA_start = int(single_isoform[2])
        ncRNA_end = int(single_isoform[3])
        single_isoseq_exons = isoseq_exons[isoform_overlap]
        #ncRNAs can only be in an intron if there is more than 1 exon
        exon_counts = 0
        if len(single_isoseq_exons) > 1:
            x = 1
            while x < len(single_isoseq_exons)-1:
                single_exon = single_isoseq_exons[x]
                single_exon_start = int(single_exon[3])
                single_exon_end = int(single_exon[4])
                single_strand = single_exon[5]
                if single_strand == "+":
                    #if the ncRNA is in exon
                    if ncRNA_start <= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end <= single_exon_end:
                        exon_counts ++ 1
                    #if the ncRNA overlaps the beginning of the exon
                    elif ncRNA_start >= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end <= single_exon_end:
                        exon_counts += 1
                    #if the ncRNA overlaps the end of the exon
                    elif ncRNA_start >= single_exon_start and ncRNA_start <= single_exon_end and ncRNA_end >= single_exon_start and ncRNA_end >= single_exon_end:
                        exon_counts += 1
                elif single_strand == "-":
                    #if ncRNA is completely within an ensembl transcripts
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the beginning of a transcript
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    #if ncRNA overlaps the end of a transcript
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                x += 1
            if counts == 0:
                #just need to check if ncRNA overlaps first or last exon at all
                last_counts = 0
                first_exon = single_isoseq_exons[0]
                last_exon = single_isoseq_exons[x]
                first_exon_start = int(first_exon[3])
                first_exon_end = int(first_exon[4])
                last_exon_start = int(last_exon[3])
                last_exon_end = int(last_exon[4])
                #checking for overlap with first_exon
                #+ strand
                if first_exon_start <= ncRNA_start and first_exon_end <= ncRNA_end and first_exon_start <= ncRNA_end and first_exon_end <= ncRNA_start:
                    last_counts += 1
                #- strand
                elif first_exon_start >= ncRNA_start and first_exon_end >= ncRNA_end and first_exon_start >= ncRNA_end and first_exon_end >= ncRNA_start:
                    last_counts += 1
                #checking for overlap with teh last exon
                #+ strand
                if last_exon_start >= ncRNA_start and last_exon_end >= ncRNA_end and last_exon_start >= ncRNA_end and last_exon_end >= ncRNA_start:
                    last_counts += 1
                #- strand
                elif last_exon_start <= ncRNA_start and last_exon_end <= ncRNA_end and last_exon_start <= ncRNA_end and last_exon_end <= ncRNA_start:
                    last_counts += 1
                if last_counts == 2:
                    isoform_id = single_isoform[1]
                    final_intronic_isoseq_ncRNAs.update({isoform_id:single_isoform})
    return final_intronic_isoseq_ncRNAs


#pulling the intronic ncRNAs
#the ncRNA will be within an intron of an ensembl or isoseq transcript
#returns dictionary with ncRNAs that are located in an intron
def combine_intronic_ncRNAs():
    ensembl_intronic_ncRNAs = compare_ncRNAs_ensembl_exons()
    isoseq_intronic_ncRNAs = compare_ncRNAs_isoseq_exons()
    all_intronic_ncRNAs = {}
    final_intronic_ncRNAs = {}
    print("Total number of ncRNAs that are located in an intron of an ensembl transcript")
    print(len(ensembl_intronic_ncRNAs))
    print("Total number of ncRNAs that are located in an intron of an isoseq transcripts")
    print(len(isoseq_intronic_ncRNAs))
    all_intronic_ncRNAs.update(ensembl_intronic_ncRNAs)
    all_intronic_ncRNAs.update(isoseq_intronic_ncRNAs)
    for key in all_intronic_ncRNAs:
        single_key = all_intronic_ncRNAs[key]
        if len(single_key) == 5 and isinstance(single_key[0], list) == False:
            final_intronic_ncRNAs.update({key:single_key})
    print("\n")
    return final_intronic_ncRNAs



#how many ncRNAs are antisense?
#this means they overlap a transcript, but on the opposite strand

def compare_antisense_ensembl_transcripts():
    ensembl_transcripts = remove_ncRNAs_ensembl_transcripts()
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    antisense_ensembl_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        for transcript in ensembl_transcripts:
            single_ensembl_transcript = ensembl_transcripts[transcript]
            ensembl_gene_id = single_ensembl_transcript[0]
            ensembl_chr_num = single_ensembl_transcript[2]
            ensembl_start_pos = int(single_ensembl_transcript[3])
            ensembl_end_pos = int(single_ensembl_transcript[4])
            ensembl_strand = single_ensembl_transcript[5]
            counts = 0
            #if ncRNA overlaps at all with the ensembl transcript counts will increase by 1
            if ncRNA_chr_num == ensembl_chr_num and ncRNA_strand != ensembl_strand:
                if ensembl_strand == "+" and ncRNA_strand == "-":
                    if ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                    elif ncRNA_start_pos >= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos <= ensembl_end_pos:
                        counts += 1
                elif ensembl_strand == "-" and ncRNA_strand == "+":
                    if ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos <= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos <= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
                    elif ncRNA_start_pos <= ensembl_start_pos and ncRNA_start_pos >= ensembl_end_pos and ncRNA_end_pos >= ensembl_start_pos and ncRNA_end_pos >= ensembl_end_pos:
                        counts += 1
            if counts > 0:
                antisense_ensembl_ncRNAs.update({isoform:single_ncRNA_isoform})
    return antisense_ensembl_ncRNAs


#examining isoseq transcripts
def compare_antisense_isoseq_transcripts():
    isoseq_transcripts = remove_ncRNAs_isoseq_transcripts()
    ncRNAs_transcripts = get_ncRNAs_transcript_positions()
    antisense_isoseq_ncRNAs = {}
    for isoform in ncRNAs_transcripts:
        single_ncRNA_isoform = ncRNAs_transcripts[isoform]
        ncRNA_chr_num = single_ncRNA_isoform[0]
        ncRNA_start_pos = int(single_ncRNA_isoform[2])
        ncRNA_end_pos = int(single_ncRNA_isoform[3])
        ncRNA_strand = single_ncRNA_isoform[4]
        for transcript in isoseq_transcripts:
            single_isoseq_transcript = isoseq_transcripts[transcript]
            isoseq_gene_id = single_isoseq_transcript[0]
            isoseq_chr_num = single_isoseq_transcript[1]
            isoseq_start_pos = int(single_isoseq_transcript[2])
            isoseq_end_pos = int(single_isoseq_transcript[3])
            isoseq_strand = single_isoseq_transcript[4]
            counts = 0
            #if ncRNA overlaps at all with the ensembl transcript counts will increase by 1
            if ncRNA_chr_num == isoseq_chr_num and ncRNA_strand != isoseq_strand:
                if isoseq_strand == "+" and ncRNA_strand == "-":
                    if ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    elif ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                    elif ncRNA_start_pos >= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos <= isoseq_end_pos:
                        counts += 1
                elif isoseq_strand == "-" and ncRNA_strand == "+":
                    if ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos <= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
                    elif ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos <= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
                    elif ncRNA_start_pos <= isoseq_start_pos and ncRNA_start_pos >= isoseq_end_pos and ncRNA_end_pos >= isoseq_start_pos and ncRNA_end_pos >= isoseq_end_pos:
                        counts += 1
            if counts > 0:
                antisense_isoseq_ncRNAs.update({isoform:single_ncRNA_isoform})
    return antisense_isoseq_ncRNAs


#pulling the antisense ncRNAs
#the ncRNA is antisense to the ensembl or isoseq transcript
#returns dictionary with ncRNAs that are antisense
def combine_antisense_ncRNAs():
    ensembl_antisense_ncRNAs = compare_antisense_ensembl_transcripts()
    isoseq_antisense_ncRNAs = compare_antisense_isoseq_transcripts()
    all_antisense_ncRNAs = {}
    final_antisense_ncRNAs = {}
    print("Total number of ncRNAs that are antisense to a ensembl transcript")
    print(len(ensembl_antisense_ncRNAs))
    print("Total number of ncRNAs that are antisense to a isoseq transcripts")
    print(len(isoseq_antisense_ncRNAs))
    all_antisense_ncRNAs.update(ensembl_antisense_ncRNAs)
    all_antisense_ncRNAs.update(isoseq_antisense_ncRNAs)
    for key in all_antisense_ncRNAs:
        single_key = all_antisense_ncRNAs[key]
        if len(single_key) == 5 and isinstance(single_key[0], list) == False:
            final_antisense_ncRNAs.update({key:single_key})
    print("\n")
    return final_antisense_ncRNAs


#writing all categories:
def write_all():
    all_ncRNAs = read_ncRNAs()
    intergenic = combine_intergenic_ncRNAs()
    intronic = combine_intronic_ncRNAs()
    antisense = combine_antisense_ncRNAs()
    remaining_isoforms = []
    intergenic_output = sys.argv[6]
    intronic_output = sys.argv[7]
    antisense_output = sys.argv[8]
    remaining_output = sys.argv[9]
    shared_intergenic_intronic_antisense = []
    shared_intergenic_intronic = []
    shared_intergenic_antisense = []
    shared_intronic_antisense = []
    #because of requiring the intergenic and intronic categories to have the same strand, there may be some isoforms in antisense that are also in intergenic or intronic categories (these need to be removed from the intergenic and intronic categories)
    for isoform_a in antisense:
        if isoform_a in intergenic:
            del intergenic[isoform_a]
        elif isoform_a in intronic:
            del intronic[isoform_a]
    for isoform in all_ncRNAs:
        if isoform not in intergenic and isoform not in intronic and isoform not in antisense:
            remaining_isoforms.append(isoform)
    print("Final number of ncRNAs that are intergenic")
    print(len(intergenic))
    print("Final number of ncRNAs that are intronic")
    print(len(intronic))
    print("Final number of ncRNAs that are antisense")
    print(len(antisense))
    print("Number of ncRNAs that either overlap a transcript in the same strand, either with or without an exon:")
    print(len(remaining_isoforms))
    with open(intergenic_output, 'a') as intergenic_out, open(intronic_output, 'a') as intronic_out, open(antisense_output, 'a') as antisense_out, open(remaining_output, 'a') as remaining_out:
        for key in intergenic:
            final = "%s\n" % str(key)
            intergenic_out.write(final)
        for key2 in intronic:
            final = "%s\n" % str(key2)
            intronic_out.write(final)
        for key3 in antisense:
            final = "%s\n" % str(key3)
            antisense_out.write(final)
        for key4 in remaining_isoforms:
            final = "%s\n" % str(key4)
            remaining_out.write(final)

write_all()

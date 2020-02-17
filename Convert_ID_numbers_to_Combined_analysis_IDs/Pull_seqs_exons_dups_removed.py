#need to pull fasta sequences, protein sequences, and exons from combined sexes classification with duplicates removed
#to run script: python3 Pull_seqs_exons_dups_removed.py <no dups classification file> <fasta file> <protein faa file> <exons gtf> <output filtered fasta file> <output filtered protein faa file> <output filtered exons gtf file>
#Author: Alice Naftaly, February 2020



import sys

#read in classification file
#just need the isoform ids
#returns list of isoform ids
def pull_isoform_ids():
    class_file = sys.argv[1]
    isoforms = []
    with open(class_file, 'r') as class_input:
        for line in class_input:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoforms.append(new_line[0])
    return isoforms


#read in fasta file
#returns fasta dictionary with key = isoform id and value == header line for sequence, followed by sequence with \n kept
def read_fasta():
    fasta_file = sys.argv[2]
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                isoform_id = new_line[0].strip(">")
                fasta_dict.update({isoform_id:[line]})
            else:
                fasta_dict[isoform_id].append(line)
    return fasta_dict


#read in protein faa file
#returns dictionary with key =isoform id and value == header line for sequence, followed by sequece with \n kept
def read_protein_faa():
    protein_file = sys.argv[3]
    protein_dict = {}
    with open(protein_file, 'r') as protein:
        for line in protein:
            if line.startswith(">"):
                new_line = line.split("\t")
                isoform_id = new_line[0].strip(">")
                protein_dict.update({isoform_id:[line]})
            else:
                protein_dict[isoform_id].append(line)
    return protein_dict

#read in exon gtf
#returns dictionary with key = isoform id and value == each exon present for each isoform
def read_exon_gtf():
    gtf_file = sys.argv[4]
    gtf_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            transcript_info = new_line[8].split(" ")
            transcript_id_partial = transcript_info[1].strip(";")
            transcript_id_final = transcript_id_partial.strip("\"")
            if transcript_id_final in gtf_dict:
                gtf_dict[transcript_id_final].append(line)
            elif transcript_id_final not in gtf_dict:
                gtf_dict.update({transcript_id_final:[line]})
    return gtf_dict


#filter fasta, faa, and exon gtf files

def filter_fasta():
    isoforms = pull_isoform_ids()
    fasta_dict = read_fasta()
    isos_to_remove = []
    for iso in fasta_dict:
        if iso not in isoforms:
            isos_to_remove.append(iso)
    for i in isos_to_remove:
        if i in fasta_dict:
            del fasta_dict[i]
    return fasta_dict


def filter_protein():
    isoforms = pull_isoform_ids()
    proteins_dict = read_protein_faa()
    isos_to_remove = []
    for iso in proteins_dict:
        if iso not in isoforms:
            isos_to_remove.append(iso)
    for i in isos_to_remove:
        if i in proteins_dict:
            del proteins_dict[i]
    return proteins_dict

def filter_gtf():
    isoforms = pull_isoform_ids()
    gtf_dict = read_exon_gtf()
    isos_to_remove = []
    for iso in gtf_dict:
        if iso not in isoforms:
            isos_to_remove.append(iso)
    for i in isos_to_remove:
        if i in gtf_dict:
            del gtf_dict[i]
    return gtf_dict


#writing duplicate removed fasta, faa, and gtf files

def write_fasta():
    filtered_fasta = filter_fasta()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for iso in filtered_fasta:
            single_iso = filtered_fasta[iso]
            for single in single_iso:
                out.write(single)


def write_proteins():
    filtered_proteins = filter_protein()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for iso in filtered_proteins:
            single_iso = filtered_proteins[iso]
            for single in single_iso:
                out.write(single)


def write_gtf():
    filtered_gtf = filter_gtf()
    output = sys.argv[7]
    with open(output, 'a') as out:
        for iso in filtered_gtf:
            single_iso = filtered_gtf[iso]
            if len(single_iso) == 1:
                out.write(single_iso[0])
            elif len(single_iso) > 1:
                for single in single_iso:
                    out.write(single)


#write all new files
def write_all():
    fasta = write_fasta()
    proteins = write_proteins()
    gtf = write_gtf()

write_all()

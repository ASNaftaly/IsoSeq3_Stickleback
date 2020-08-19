#sort through concatenated features file to make sure all isoforms are accounted for correctly
#read in concatenated file and then sort by feature
#to run script: python3 Check.Corrected.Features.py <concat file> <classification file to get relationship between isoform and gene ids> <bed file with all isoforms> <weird isoforms corrected bed file> <isoseq gtf file> <output gtf file>
#Author: Alice Naftaly, Aug 2020

import sys

#read in concat file
#returns dictionary with key == feature and value == line
def read_concat_file():
    input_file = sys.argv[1]
    concat_dict = {}
    with open(input_file, 'r') as concat:
        for line in concat:
            new_line = line.split("\t")
            feature = new_line[2]
            if feature in concat_dict:
                concat_dict[feature].append(line)
            elif feature not in concat_dict:
                concat_dict.update({feature:[line]})
    return concat_dict


#read classification file
#returns dictionary with key == isoform id and value == gene id
def read_class_genes():
    class_file = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                gene_id = new_line[6]
                class_dict.update({isoform_id:gene_id})
    return class_dict

#returns dictionary with key == isoform id and value == gene id
def read_class_protein_coding():
    class_file = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                protein_coding = new_line[27]
                if protein_coding == "coding":
                    final = "protein_coding"
                elif protein_coding == "non_coding":
                    final = "nonprotein_coding"
                class_dict.update({isoform_id:final})
    return class_dict

#read bed file
def read_bed_files():
    bed_file_full = sys.argv[3]
    bed_file_weird = sys.argv[4]
    bed_full_dict = {}
    bed_weird_dict = {}
    bed_final_dict = {}
    with open(bed_file_full, 'r') as bed_full, open(bed_file_weird, 'r') as bed_weird:
        for line in bed_full:
            new_line = line.split()
            bed_full_dict.update({new_line[3]:[new_line[0], new_line[1], new_line[2], new_line[5]]})
        for l in bed_weird:
            new_line = l.split()
            bed_weird_dict.update({new_line[3]:[new_line[0], new_line[1], new_line[2], new_line[5]]})
    for key in bed_full_dict:
        if key in bed_weird_dict:
            bed_final_dict.update({key:bed_weird_dict[key]})
        elif key not in bed_weird_dict:
            bed_final_dict.update({key:bed_full_dict[key]})
    return bed_final_dict

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
    return exon_dict

def chrs_scaffolds():
    input_file = sys.argv[6]
    chr_list = []
    with open(input_file, 'r') as chrs:
        for line in chrs:
            new_line = line.strip("\n")
            chr_list.append(new_line)
    return chr_list


#begin to check concat file
#gene Features
#added isoform id (PB.XXXX.xx = PB.XXXX) to gene info
#returns dictionary with key == isoform and value == full gtf line
def check_gene_features():
    isoform_dict = read_class_genes()
    concat_file = read_concat_file()
    genes = concat_file["gene"]
    final_gene_dict = {}
    for gene in genes:
        new_line = gene.split("\t")
        chr = new_line[0]
        gene_info = new_line[8].split(";")
        gene_id = gene_info[0].split(" ")[1]
        protein_coding = gene_info[2].split(" ")[2]
        final_protein_coding = protein_coding.strip("'")
        for isoform in isoform_dict:
            isoform_gene = isoform_dict[isoform]
            if isoform_gene == gene_id:
                isoform_id = isoform
        full_isoform_gene_id = isoform_id.split(".")
        final_isoform_gene_id = full_isoform_gene_id[0] + "." + full_isoform_gene_id[1]
        final_gene_info_entry = " isoform_gene_id %s" % str(final_isoform_gene_id)
        final_gene_biotype = " gene_biotype %s;" % str(final_protein_coding)
        if len(gene_info) == 4:
            final_gene_feature = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;%s;%s;%s;\n" % (new_line[0], new_line[1], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], gene_info[0], gene_info[1], final_gene_biotype, final_gene_info_entry)
        elif len(gene_info) == 5:
            final_gene_feature = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;%s;%s;%s;%s;\n" % (new_line[0], new_line[1], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], gene_info[0], gene_info[1], final_gene_biotype, gene_info[3], final_gene_info_entry)
        if chr in final_gene_dict:
            final_gene_dict[chr].append(final_gene_feature)
        elif chr not in final_gene_dict:
            final_gene_dict.update({chr:[final_gene_feature]})
    return final_gene_dict

#transcript Features
#returns dictionary with key == isoform id and value == transcript feature gtf line
def check_transcript_features():
    concat_file = read_concat_file()
    transcripts = concat_file["transcript"]
    final_transcript_dict = {}
    for transcript in transcripts:
        new_line = transcript.split("\t")
        transcript_info = new_line[8].split(";")
        for val in transcript_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
        final_transcript_dict.update({isoform_id:transcript})
    return final_transcript_dict

#exon Features
#returns dictionary with key == isoform id and value == exon feature for gtf
def check_exon_features():
    concat_file = read_concat_file()
    protein_coding = read_class_protein_coding()
    exons = concat_file["exon"]
    final_exon_dict = {}
    exon_dict_finalized = {}
    for exon in exons:
        new_line = exon.split("\t")
        exon_info = new_line[8].split(";")
        for val in exon_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
        if isoform_id in final_exon_dict:
            final_exon_dict[isoform_id].append(exon)
        elif isoform_id not in final_exon_dict:
            final_exon_dict.update({isoform_id:[exon]})
    for key in final_exon_dict:
        set_key = list(set(final_exon_dict[key]))
        exon_numbers = []
        for x in set_key:
            info = x.split("\t")[8]
            info_2 = info.split(";")
            for v in info_2:
                if v.startswith(" exon_number"):
                    exon_numbers.append(v.split(" ")[2])
        if len(exon_numbers) != len(list(set(exon_numbers))):
            for exon in set_key:
                full_gene_info = exon.split("\t")[8]
                for i in full_gene_info.split(";"):
                    if i.startswith(" gene_biotype"):
                        split_i = i.split(" ")[2]
                if split_i == protein_coding[key]:
                    if key in exon_dict_finalized:
                        exon_dict_finalized[key].append(exon)
                    elif key not in exon_dict_finalized:
                        exon_dict_finalized.update({key:[exon]})
        elif len(exon_numbers) == len(list(set(exon_numbers))):
            exon_dict_finalized.update({key:set_key})
    return exon_dict_finalized

#fix minus strand first exon:
def fix_exons():
    exons = check_exon_features()
    bed_dict = read_bed_files()
    gtf_dict = read_isoseq_gtf()
    final_exons = {}
    count = 0
    for iso in exons:
        single_iso = exons[iso]
        for exon in single_iso:
            split_exon = exon.split("\t")
            if split_exon[6] == "-":
                iso_info = split_exon[8].split(" ")
                exon_number_index = iso_info.index("exon_number")
                if iso_info[exon_number_index+1] == "1;":
                    chr = split_exon[0]
                    source = split_exon[1]
                    feature = split_exon[2]
                    start = split_exon[3]
                    end = split_exon[4]
                    holder1 = split_exon[5]
                    strand = split_exon[6]
                    holder2 = split_exon[7]
                    gene_info = split_exon[8]
                    transcript_position = bed_dict[iso]
                    gtf_position = gtf_dict[iso]
                    final_start = gtf_position[len(gtf_position)-1][1]
                    final_end = transcript_position[2]
                    final = [chr, source, feature, final_start, final_end, holder1, strand, holder2, gene_info]
                    joined_final = "\t".join(final)
                    count += 1
                    if iso in final_exons:
                        final_exons[iso].append(joined_final)
                    elif iso not in final_exons:
                        final_exons.update({iso:[joined_final]})
                elif iso_info[exon_number_index+1] != "1;":
                    count +=1
                    if iso in final_exons:
                        final_exons[iso].append(exon)
                    elif iso not in final_exons:
                        final_exons.update({iso:[exon]})
            else:
                count +=1
                if iso in final_exons:
                    final_exons[iso].append(exon)
                elif iso not in final_exons:
                    final_exons.update({iso:[exon]})
    return final_exons

#start codon Feature
#returns dictionary with key == isoform id and value == start codon feature for gtf
def check_start_codon_features():
    concat_file = read_concat_file()
    start_codons = concat_file["start_codon"]
    final_start_codon_dict = {}
    for value in start_codons:
        new_value = value.split("\t")
        gene_info = new_value[8].split(";")
        for val in gene_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
                final_isoform_id = isoform_id.strip("'")
        if final_isoform_id in final_start_codon_dict:
            final_start_codon_dict[final_isoform_id].append(value)
        elif final_isoform_id not in final_start_codon_dict:
            final_start_codon_dict.update({final_isoform_id:[value]})
    return final_start_codon_dict



#stop codon Feature
#returns dictionary with key == isoform id and value == stop codon feature for gtf
def check_stop_codon_features():
    concat_file = read_concat_file()
    stop_codons = concat_file["stop_codon"]
    final_stop_codon_dict = {}
    for value in stop_codons:
        new_value = value.split("\t")
        gene_info = new_value[8].split(";")
        for val in gene_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
                final_isoform_id = isoform_id.strip("'")
        if final_isoform_id in final_stop_codon_dict:
            final_stop_codon_dict[final_isoform_id].append(value)
        elif final_isoform_id not in final_stop_codon_dict:
            final_stop_codon_dict.update({final_isoform_id:[value]})
    return final_stop_codon_dict


#CDS Feature
#returns dictionary with key == isoform id and value == cds feature for gtf
def check_CDS_features():
    concat_file = read_concat_file()
    cds_region = concat_file["CDS"]
    final_cds_dict = {}
    for value in cds_region:
        new_value = value.split("\t")
        gene_info = new_value[8].split(";")
        for val in gene_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
                final_isoform_id = isoform_id.strip("'")
        if final_isoform_id in final_cds_dict :
            final_cds_dict[final_isoform_id].append(value)
        elif final_isoform_id not in final_cds_dict :
            final_cds_dict.update({final_isoform_id:[value]})
    return final_cds_dict




#five prime utr Feature
#returns dictionary with key == isoform id and value == five prime utr feature for gtf
def check_five_prime_utr_features():
    concat_file = read_concat_file()
    five_prime_utr = concat_file["five_prime_utr"]
    final_five_prime_utr_dict = {}
    for value in five_prime_utr:
        new_value = value.split("\t")
        gene_info = new_value[8].split(";")
        for val in gene_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
                final_isoform_id = isoform_id.strip("'")
        if final_isoform_id in final_five_prime_utr_dict:
            final_five_prime_utr_dict[final_isoform_id].append(value)
        elif final_isoform_id not in final_five_prime_utr_dict:
            final_five_prime_utr_dict.update({final_isoform_id:[value]})
    return final_five_prime_utr_dict


#three prime utr Feature
#returns dictionary with key == isoform id and value == three prime utr feature for gtf
def check_three_prime_utr_features():
    concat_file = read_concat_file()
    three_prime_utr = concat_file["three_prime_utr"]
    final_three_prime_utr_dict = {}
    for value in three_prime_utr:
        new_value = value.split("\t")
        gene_info = new_value[8].split(";")
        for val in gene_info:
            if val.startswith(" isoform_id"):
                isoform_id = val.split(" ")[2]
                final_isoform_id = isoform_id.strip("'")
        if final_isoform_id in final_three_prime_utr_dict:
            final_three_prime_utr_dict[final_isoform_id].append(value)
        elif final_isoform_id not in final_three_prime_utr_dict:
            final_three_prime_utr_dict.update({final_isoform_id:[value]})
    return final_three_prime_utr_dict

#writing final file
def write_file():
    gene_features = check_gene_features()
    transcript_features = check_transcript_features()
    exon_features = fix_exons()
    start_codons = check_start_codon_features()
    stop_codons = check_stop_codon_features()
    five_prime_utrs = check_five_prime_utr_features()
    cds_regions = check_CDS_features()
    three_prime_utrs = check_three_prime_utr_features()
    chrs = chrs_scaffolds()
    output = sys.argv[7]
    transcript_list = []
    with open(output, 'a') as out:
        for chr in chrs:
            single_chr = gene_features[chr]
            for gene in single_chr:
                out.write(gene)
                split_gene = gene.split("\t")
                gene_info = split_gene[8].split(";")
                full_isoform_id = gene_info[len(gene_info)-2].split(" ")
                gene_isoform = full_isoform_id[2]
                for iso in transcript_features:
                    single_transcript = transcript_features[iso]
                    split_iso = iso.split(".")
                    single_gene_iso = split_iso[0] + "." + split_iso[1]
                    if gene_isoform == single_gene_iso:
                        if iso in transcript_list:
                            continue
                        elif iso not in transcript_list:
                            transcript_list.append(iso)
                            out.write(single_transcript)
                        single_exon_list = exon_features[iso]
                        for exon in single_exon_list:
                            out.write(exon)
                        if iso in start_codons:
                            single_start = start_codons[iso]
                            for start in single_start:
                                out.write(start)
                        if iso in stop_codons:
                            single_stop = stop_codons[iso]
                            for stop in single_stop:
                                out.write(stop)
                        if iso in cds_regions:
                            single_cds = cds_regions[iso]
                            for cds in single_cds:
                                out.write(cds)
                        if iso in five_prime_utrs:
                            single_five_prime = five_prime_utrs[iso]
                            for five_prime in single_five_prime:
                                out.write(five_prime)
                        if iso in three_prime_utrs:
                            single_three_prime = three_prime_utrs[iso]
                            for three_prime in single_three_prime:
                                out.write(three_prime)
        for key in transcript_features:
            if key not in transcript_list:
                out.write(transcript_features[key])
                single_exon_list = exon_features[key]
                for exon in single_exon_list:
                    out.write(exon)
                if key in start_codons:
                    single_start = start_codons[key]
                    for start in single_start:
                        out.write(start)
                if key in stop_codons:
                    single_stop = stop_codons[key]
                    for stop in single_stop:
                        out.write(stop)
                if key in cds_regions:
                    single_cds = cds_regions[key]
                    for cds in single_cds:
                        out.write(cds)
                if key in five_prime_utrs:
                    single_five_prime = five_prime_utrs[key]
                    for five_prime in single_five_prime:
                        out.write(five_prime)
                if key in three_prime_utrs:
                    single_three_prime = three_prime_utrs[key]
                    for three_prime in single_three_prime:
                        out.write(three_prime)


write_file()

#Analyzing FSM from isoseq and comparing tss and tts
#to run script: python3 FSM_TSS_TTS_differences.py <isoseq classification file> <ensembl classification file> <gtf file isoseq> <output tss differences> <tts differences>
#Author: Alice Naftaly, January 2020

import sys

#read in isoseq classification file
#returns dictionary with key = gene id and value == [isoform_id, chr_num, strand, length, num_exons, structural_category, gene_id, transcript_id, reference_length, ref_num_exons, diff_to_TSS, diff_to_TTS, subcategory, splice_junction_type, number_of_full_length_reads, coding_potential, orf_length, cds_length, cds_start, cds_end]
def read_isoseq_classification():
    input_file = sys.argv[1]
    isoseq_class = {}
    with open(input_file, 'r') as class_file:
        for line in class_file:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split("\t")
                isoform_id = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                length = new_line[3]
                num_exons = new_line[4]
                structural_category = new_line[5]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                reference_length = new_line[8]
                ref_num_exons = new_line[9]
                diff_to_TSS = new_line[10]
                diff_to_TTS = new_line[11]
                subcategory = new_line[12]
                splice_junction_type = new_line[14]
                number_of_full_length_reads = new_line[19]
                coding_potential = new_line[27]
                orf_length = new_line[28]
                cds_length = new_line[29]
                cds_start = new_line[30]
                cds_end = new_line[31]
                dict_value = [isoform_id, chr_num, strand, length, num_exons, structural_category, gene_id, transcript_id, reference_length, ref_num_exons, diff_to_TSS, diff_to_TTS, subcategory, splice_junction_type, number_of_full_length_reads, coding_potential, orf_length, cds_length, cds_start, cds_end]
                if gene_id in isoseq_class:
                    isoseq_class[gene_id].append(dict_value)
                elif gene_id not in isoseq_class:
                    isoseq_class.update({gene_id:[dict_value]})
    return isoseq_class

#reads ensembl classification file
#returns dictionary with key = gene id and value == [isoform_id, chr_num, strand, length, num_exons, structural_category, gene_id, transcript_id, subcategory, splice_junction_type, coding_potential, orf_length, cds_length, cds_start, cds_end]
def read_ensembl_classification():
    input_file = sys.argv[2]
    ensembl_class = {}
    with open(input_file, 'r') as class_file:
        for line in class_file:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                length = new_line[3]
                num_exons = new_line[4]
                structural_category = new_line[5]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                subcategory = new_line[12]
                splice_junction_type = new_line[14]
                coding_potential = new_line[27]
                orf_length = new_line[28]
                cds_length = new_line[29]
                cds_start = new_line[30]
                cds_end = new_line[31]
                dict_value = [isoform_id, chr_num, strand, length, num_exons, structural_category, gene_id, transcript_id, subcategory, splice_junction_type, coding_potential, orf_length, cds_length, cds_start, cds_end]
                if gene_id in ensembl_class:
                    ensembl_class[gene_id].append(dict_value)
                elif gene_id not in ensembl_class:
                    ensembl_class.update({gene_id:[dict_value]})
    return ensembl_class

#read gtf file from isoseq output
#returns dictionary with key == isoform id and value == [isoform id, chrom number, exon number, exon start, exon end]
def read_gtf_isoseq():
    gtf_file = sys.argv[3]
    gtf_dict = {}
    final_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            chrom = new_line[0]
            strand = new_line[6]
            exon_start = new_line[3]
            exon_end = new_line[4]
            isoform_id_full = new_line[9].strip(";")
            isoform_id = isoform_id_full.strip('\"')
            dict_value = [chrom, strand, exon_start, exon_end, isoform_id]
            if isoform_id in gtf_dict:
                gtf_dict[isoform_id].append(dict_value)
            elif isoform_id not in gtf_dict:
                gtf_dict.update({isoform_id:[dict_value]})
        for isoform in gtf_dict:
            single_isoform = gtf_dict[isoform]
            if len(single_isoform) == 1:
                single = single_isoform[0]
                #new_dict_value = [isoform, chrom, strand, exon_number, exon_start, exon_end]
                new_dict_value = [isoform, single[0], single[1], "1", single[2],single[3]]
                final_dict.update({isoform:new_dict_value})
            elif len(single_isoform) > 1:
                x = 1
                for value in single_isoform:
                    new_dict_value = [isoform, value[0], value[1], str(x), value[2],value[3]]
                    if isoform in final_dict:
                        final_dict[isoform].append(new_dict_value)
                        x += 1
                    elif isoform not in final_dict:
                        final_dict.update({isoform:[new_dict_value]})
                        x += 1
    return final_dict

#creating dictionary to convert between isoform id and gene id and isoform id and transcript id
def conversion_dictionary():
    isoseq_classification = read_isoseq_classification()
    conversion_dict = {}
    for gene in isoseq_classification:
        single_gene = isoseq_classification[gene]
        if len(single_gene) == 1:
            single = single_gene[0]
            gene_id = gene
            isoform_id = single[0]
            transcript_id = single[7]
            dict_value = [gene_id, transcript_id]
            conversion_dict.update({isoform_id:dict_value})
        elif len(single_gene) > 1:
            for value in single_gene:
                isoform_id = value[0]
                gene_id = gene
                transcript_id = value[7]
                dict_value = [gene_id, transcript_id]
                if isoform_id in conversion_dict:
                    conversion_dict[isoform_id].append(dict_value)
                elif isoform_id not in conversion_dict:
                    conversion_dict.update({isoform_id:[dict_value]})
    return conversion_dict


#calculating various aspects about isoseq isoforms

#examining number of isoforms with differential transcription start sites
#specifically looks at differences in start sites for full splice matches only
#returns final dictionary with key = gene id and value == [isoform id, transcript id, position different from TSS where 0 = same TSS and negative value means upstream of reference TSS and positive value means downstream of reference TSS]
def differences_start_sites():
    isoseq_class = read_isoseq_classification()
    differences_in_tss_dict = {}
    final_dict = {}
    for key in isoseq_class:
        if key.startswith("ENSGACG"):
            isoseq_value = isoseq_class[key]
            if len(isoseq_value) == 1:
                single_isoseq_value = isoseq_value[0]
                if single_isoseq_value[5] == "full-splice_match":
                    diff_to_tss = single_isoseq_value[10]
                    dict_value = [single_isoseq_value[0],single_isoseq_value[7],diff_to_tss]
                    if key in differences_in_tss_dict:
                        differences_in_tss_dict[key].append(dict_value)
                    elif key not in differences_in_tss_dict:
                        differences_in_tss_dict.update({key:[dict_value]})
            elif len(isoseq_value) > 1:
                for single in isoseq_value:
                    if single[5] == "full-splice_match":
                        diff_to_tss = single[10]
                        dict_value = [single[0],single[7],diff_to_tss]
                        if key in differences_in_tss_dict:
                            differences_in_tss_dict[key].append(dict_value)
                        elif key not in differences_in_tss_dict:
                            differences_in_tss_dict.update({key:[dict_value]})
    for gene in differences_in_tss_dict:
        single_key = differences_in_tss_dict[gene]
        if len(single_key) == 1:
            final_dict.update({gene:single_key[0]})
        elif len(single_key) > 1:
            for value in single_key:
                if gene in final_dict:
                    final_dict[gene].append(value)
                elif gene not in final_dict:
                    final_dict.update({gene:[value]})
    return final_dict




#plotting various values for differences in TSS compared to reference
#going to split differences in start sites as follows:
#0 differences = for any transcript that has no differences from reference sequence in start site
#1-5 nt differences = slight difference in start site
#6-10 nt differences
#10-25 nt differences
#25-50 nt differences
#50-100 nt differences
#100-500 nt differences
#500-1000 nt differences
#>1000 nt differences
def sort_differences_tss():
    diff_tss = differences_start_sites()
    differences = {}
    key_no_diff = "No.Differences"
    key_1_5_diff = "1.to.5.nt.Diff"
    key_5_10_diff = "5.to.10.nt.Diff"
    key_10_25_diff = "10.to.25.nt.Diff"
    key_25_50_diff = "25.to.50.nt.Diff"
    key_50_100_diff = "50.to.100.nt.Diff"
    key_100_500_diff = "100.to.500.nt.Diff"
    key_500_1000_diff = "500.to.1000.nt.Diff"
    key_greater_1000_diff = "Greater.than.1000.nt.Diff"
    for key in diff_tss:
        single_key = diff_tss[key]
        if len(single_key[0]) == 3:
            for single in single_key:
                #can be either upstream (-) or downstream (+); need to take absolute value of difference to tss
                diff = abs(int(single[2]))
                dict_value = [key, single[1],single[0],single[2]]
                if diff == 0:
                    if key_no_diff in differences:
                        differences[key_no_diff].append(dict_value)
                    elif key_no_diff not in differences:
                        differences.update({key_no_diff:[dict_value]})
                elif 1 <= diff <= 5:
                    if key_1_5_diff in differences:
                        differences[key_1_5_diff].append(dict_value)
                    elif key_1_5_diff not in differences:
                        differences.update({key_1_5_diff:[dict_value]})
                elif 5 < diff <= 10:
                    if key_5_10_diff in differences:
                        differences[key_5_10_diff].append(dict_value)
                    elif key_5_10_diff not in differences:
                        differences.update({key_5_10_diff:[dict_value]})
                elif 10 < diff <= 25:
                    if key_10_25_diff in differences:
                        differences[key_10_25_diff].append(dict_value)
                    elif key_10_25_diff not in differences:
                        differences.update({key_10_25_diff:[dict_value]})
                elif 25 < diff <= 50:
                    if key_25_50_diff in differences:
                        differences[key_25_50_diff].append(dict_value)
                    elif key_25_50_diff not in differences:
                        differences.update({key_25_50_diff:[dict_value]})
                elif 50 < diff <= 100:
                    if key_50_100_diff in differences:
                        differences[key_50_100_diff].append(dict_value)
                    elif key_50_100_diff not in differences:
                        differences.update({key_50_100_diff:[dict_value]})
                elif 100 < diff <= 500:
                    if key_100_500_diff in differences:
                        differences[key_100_500_diff].append(dict_value)
                    elif key_100_500_diff not in differences:
                        differences.update({key_100_500_diff:[dict_value]})
                elif 500 < diff <= 1000:
                    if key_500_1000_diff in differences:
                        differences[key_500_1000_diff].append(dict_value)
                    elif key_500_1000_diff not in differences:
                        differences.update({key_500_1000_diff:[dict_value]})
                elif diff > 1000:
                    if key_greater_1000_diff in differences:
                        differences[key_greater_1000_diff].append(dict_value)
                    elif key_greater_1000_diff not in differences:
                        differences.update({key_greater_1000_diff:[dict_value]})
        else:
            #can be either upstream (-) or downstream (+); need to take absolute value of difference to tss
            diff = abs(int(single_key[2]))
            dict_value = [key, single_key[1],single_key[0],single_key[2]]
            if diff == 0:
                if key_no_diff in differences:
                    differences[key_no_diff].append(dict_value)
                elif key_no_diff not in differences:
                    differences.update({key_no_diff:[dict_value]})
            elif 1 <= diff <= 5:
                if key_1_5_diff in differences:
                    differences[key_1_5_diff].append(dict_value)
                elif key_1_5_diff not in differences:
                    differences.update({key_1_5_diff:[dict_value]})
            elif 5 < diff <= 10:
                if key_5_10_diff in differences:
                    differences[key_5_10_diff].append(dict_value)
                elif key_5_10_diff not in differences:
                    differences.update({key_5_10_diff:[dict_value]})
            elif 10 < diff <= 25:
                if key_10_25_diff in differences:
                    differences[key_10_25_diff].append(dict_value)
                elif key_10_25_diff not in differences:
                    differences.update({key_10_25_diff:[dict_value]})
            elif 25 < diff <= 50:
                if key_25_50_diff in differences:
                    differences[key_25_50_diff].append(dict_value)
                elif key_25_50_diff not in differences:
                    differences.update({key_25_50_diff:[dict_value]})
            elif 50 < diff <= 100:
                if key_50_100_diff in differences:
                    differences[key_50_100_diff].append(dict_value)
                elif key_50_100_diff not in differences:
                    differences.update({key_50_100_diff:[dict_value]})
            elif 100 < diff <= 500:
                if key_100_500_diff in differences:
                    differences[key_100_500_diff].append(dict_value)
                elif key_100_500_diff not in differences:
                    differences.update({key_100_500_diff:[dict_value]})
            elif 500 < diff <= 1000:
                if key_500_1000_diff in differences:
                    differences[key_500_1000_diff].append(dict_value)
                elif key_500_1000_diff not in differences:
                    differences.update({key_500_1000_diff:[dict_value]})
            elif diff > 1000:
                if key_greater_1000_diff in differences:
                    differences[key_greater_1000_diff].append(dict_value)
                elif key_greater_1000_diff not in differences:
                    differences.update({key_greater_1000_diff:[dict_value]})
    return differences

#writing TSS differences to file
def write_tss_diff():
    diff_dict = sort_differences_tss()
    differences_types = ["No.Differences", "1.to.5.nt.Diff", "5.to.10.nt.Diff", "10.to.25.nt.Diff", "25.to.50.nt.Diff", "50.to.100.nt.Diff", "100.to.500.nt.Diff", "500.to.1000.nt.Diff", "Greater.than.1000.nt.Diff"]
    output = sys.argv[4]
    with open(output,'w') as out:
        header = "Difference.Type\tGene.ID\tTranscript.ID\tIsoform.ID\tDist.to.TSS\n"
        out.write(header)
        for diff_type in differences_types:
            single_type = diff_dict[diff_type]
            for value in single_type:
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(diff_type),str(value[0]),str(value[1]),str(value[2]),str(value[3]))
                out.write(final)

#examining number of isoforms with differential transcription end sites
#specifically looks at differences in end sites for full splice matches only
#returns final dictionary with key = gene id and value == [isoform id, transcript id, position different from TTS where 0 = same TTS and negative value means upstream of reference TTS and positive value means downstream of reference TTS]
def differences_end_sites():
    isoseq_class = read_isoseq_classification()
    differences_in_tts_dict = {}
    final_dict = {}
    for key in isoseq_class:
        if key.startswith("ENSGACG"):
            isoseq_value = isoseq_class[key]
            if len(isoseq_value) == 1:
                single_isoseq_value = isoseq_value[0]
                if single_isoseq_value[5] == "full-splice_match":
                    diff_to_tts = single_isoseq_value[11]
                    dict_value = [single_isoseq_value[0],single_isoseq_value[7],diff_to_tts]
                    if key in differences_in_tts_dict:
                        differences_in_tts_dict[key].append(dict_value)
                    elif key not in differences_in_tts_dict:
                        differences_in_tts_dict.update({key:[dict_value]})
            elif len(isoseq_value) > 1:
                for single in isoseq_value:
                    if single[5] == "full-splice_match":
                        diff_to_tts = single[11]
                        dict_value = [single[0],single[7],diff_to_tts]
                        if key in differences_in_tts_dict:
                            differences_in_tts_dict[key].append(dict_value)
                        elif key not in differences_in_tts_dict:
                            differences_in_tts_dict.update({key:[dict_value]})
    for gene in differences_in_tts_dict:
        single_key = differences_in_tts_dict[gene]
        if len(single_key) == 1:
            final_dict.update({gene:single_key[0]})
        elif len(single_key) > 1:
            for value in single_key:
                if gene in final_dict:
                    final_dict[gene].append(value)
                elif gene not in final_dict:
                    final_dict.update({gene:[value]})
    return final_dict


#plotting various values for differences in TTS compared to reference
#going to split differences in end sites as follows:
#0 differences = for any transcript that has no differences from reference sequence in end site
#1-5 nt differences = slight difference in end site
#6-10 nt differences
#10-25 nt differences
#25-50 nt differences
#50-100 nt differences
#100-500 nt differences
#500-1000 nt differences
#>1000 nt differences
def sort_differences_tts():
    diff_tts = differences_end_sites()
    differences = {}
    key_no_diff = "No.Differences"
    key_1_5_diff = "1.to.5.nt.Diff"
    key_5_10_diff = "5.to.10.nt.Diff"
    key_10_25_diff = "10.to.25.nt.Diff"
    key_25_50_diff = "25.to.50.nt.Diff"
    key_50_100_diff = "50.to.100.nt.Diff"
    key_100_500_diff = "100.to.500.nt.Diff"
    key_500_1000_diff = "500.to.1000.nt.Diff"
    key_greater_1000_diff = "Greater.than.1000.nt.Diff"
    for key in diff_tts:
        single_key = diff_tts[key]
        if len(single_key[0]) == 3:
            for single in single_key:
                #can be either upstream (-) or downstream (+); need to take absolute value of difference to tss
                diff = abs(int(single[2]))
                dict_value = [key, single[1],single[0],single[2]]
                if diff == 0:
                    if key_no_diff in differences:
                        differences[key_no_diff].append(dict_value)
                    elif key_no_diff not in differences:
                        differences.update({key_no_diff:[dict_value]})
                elif 1 <= diff <= 5:
                    if key_1_5_diff in differences:
                        differences[key_1_5_diff].append(dict_value)
                    elif key_1_5_diff not in differences:
                        differences.update({key_1_5_diff:[dict_value]})
                elif 5 < diff <= 10:
                    if key_5_10_diff in differences:
                        differences[key_5_10_diff].append(dict_value)
                    elif key_5_10_diff not in differences:
                        differences.update({key_5_10_diff:[dict_value]})
                elif 10 < diff <= 25:
                    if key_10_25_diff in differences:
                        differences[key_10_25_diff].append(dict_value)
                    elif key_10_25_diff not in differences:
                        differences.update({key_10_25_diff:[dict_value]})
                elif 25 < diff <= 50:
                    if key_25_50_diff in differences:
                        differences[key_25_50_diff].append(dict_value)
                    elif key_25_50_diff not in differences:
                        differences.update({key_25_50_diff:[dict_value]})
                elif 50 < diff <= 100:
                    if key_50_100_diff in differences:
                        differences[key_50_100_diff].append(dict_value)
                    elif key_50_100_diff not in differences:
                        differences.update({key_50_100_diff:[dict_value]})
                elif 100 < diff <= 500:
                    if key_100_500_diff in differences:
                        differences[key_100_500_diff].append(dict_value)
                    elif key_100_500_diff not in differences:
                        differences.update({key_100_500_diff:[dict_value]})
                elif 500 < diff <= 1000:
                    if key_500_1000_diff in differences:
                        differences[key_500_1000_diff].append(dict_value)
                    elif key_500_1000_diff not in differences:
                        differences.update({key_500_1000_diff:[dict_value]})
                elif diff > 1000:
                    if key_greater_1000_diff in differences:
                        differences[key_greater_1000_diff].append(dict_value)
                    elif key_greater_1000_diff not in differences:
                        differences.update({key_greater_1000_diff:[dict_value]})
        else:
            diff = abs(int(single_key[2]))
            dict_value = [key, single_key[1],single_key[0],single_key[2]]
            if diff == 0:
                if key_no_diff in differences:
                    differences[key_no_diff].append(dict_value)
                elif key_no_diff not in differences:
                    differences.update({key_no_diff:[dict_value]})
            elif 1 <= diff <= 5:
                if key_1_5_diff in differences:
                    differences[key_1_5_diff].append(dict_value)
                elif key_1_5_diff not in differences:
                    differences.update({key_1_5_diff:[dict_value]})
            elif 5 < diff <= 10:
                if key_5_10_diff in differences:
                    differences[key_5_10_diff].append(dict_value)
                elif key_5_10_diff not in differences:
                    differences.update({key_5_10_diff:[dict_value]})
            elif 10 < diff <= 25:
                if key_10_25_diff in differences:
                    differences[key_10_25_diff].append(dict_value)
                elif key_10_25_diff not in differences:
                    differences.update({key_10_25_diff:[dict_value]})
            elif 25 < diff <= 50:
                if key_25_50_diff in differences:
                    differences[key_25_50_diff].append(dict_value)
                elif key_25_50_diff not in differences:
                    differences.update({key_25_50_diff:[dict_value]})
            elif 50 < diff <= 100:
                if key_50_100_diff in differences:
                    differences[key_50_100_diff].append(dict_value)
                elif key_50_100_diff not in differences:
                    differences.update({key_50_100_diff:[dict_value]})
            elif 100 < diff <= 500:
                if key_100_500_diff in differences:
                    differences[key_100_500_diff].append(dict_value)
                elif key_100_500_diff not in differences:
                    differences.update({key_100_500_diff:[dict_value]})
            elif 500 < diff <= 1000:
                if key_500_1000_diff in differences:
                    differences[key_500_1000_diff].append(dict_value)
                elif key_500_1000_diff not in differences:
                    differences.update({key_500_1000_diff:[dict_value]})
            elif diff > 1000:
                if key_greater_1000_diff in differences:
                    differences[key_greater_1000_diff].append(dict_value)
                elif key_greater_1000_diff not in differences:
                    differences.update({key_greater_1000_diff:[dict_value]})
    return differences


#writing TSS differences to file
def write_tts_diff():
    diff_dict = sort_differences_tts()
    differences_types = ["No.Differences", "1.to.5.nt.Diff", "5.to.10.nt.Diff", "10.to.25.nt.Diff", "25.to.50.nt.Diff", "50.to.100.nt.Diff", "100.to.500.nt.Diff", "500.to.1000.nt.Diff", "Greater.than.1000.nt.Diff"]
    output = sys.argv[5]
    with open(output,'w') as out:
        header = "Difference.Type\tGene.ID\tTranscript.ID\tIsoform.ID\tDist.to.TSS\n"
        out.write(header)
        for diff_type in differences_types:
            single_type = diff_dict[diff_type]
            for value in single_type:
                final = "%s\t%s\t%s\t%s\t%s\n" % (str(diff_type),str(value[0]),str(value[1]),str(value[2]),str(value[3]))
                out.write(final)


#write all functions:
def write_all():
    tss_diffs = write_tss_diff()
    tts_diffs = write_tts_diff()

write_all()

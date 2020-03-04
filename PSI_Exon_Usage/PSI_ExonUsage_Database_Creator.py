#Creating sequence junction database to examine exon usage following methods of Razeto-Barry, Diaz, & Vasquez 2012 and Barbosa-Morais et al 2012
#This script works by separating each gene/transcript into exon trios where the first and third exon are the constituitive exons and the 2nd exon is the cassette exon (may be spliced out)
#   this means all single or 2 exon isoforms will be excluded, the first and last exon will not be examined for potential loss or differences
#1. Make dictionary of exon trios for both isoseq isoforms and all annotated genes from ensembl; will create separate functions so I can differentiate between only isoseq, only ensembl, and the two combined
#to run script: python3 PSI_ExonUsage_Database_Creator.py <isoseq classification file> <isoseq gtf file> <ensembl gtf file>



import sys

#read isoseq classification file
#returns dictionary with key = isoform id and value == [gene_id, transcript_id, isoform_id, chr_num, strand, num_exons]
def read_isoseq_classification():
    class_file = sys.argv[1]
    iso_class_dict = {}
    with open(class_file, 'r') as class_input:
        for line in class_input:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                chr_num = new_line[1]
                strand = new_line[2]
                num_exons = int(new_line[4])
                gene_id = new_line[6]
                transcript_id = new_line[7]
                if num_exons > 2:
                    dict_value = [gene_id, transcript_id, isoform_id, chr_num, strand, num_exons]
                    iso_class_dict.update({isoform_id:dict_value})
    return iso_class_dict


#read isoseq gtf file
#returns dictionary with key = isoform_id, value == [chr_num, strand, exon_start, exon_end]
def read_isoseq_gtf():
    gtf_file = sys.argv[2]
    gtf_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split("\t")
            chr_num = new_line[0]
            exon_start = new_line[3]
            exon_end = new_line[4]
            strand = new_line[6]
            transcript_info = new_line[8].split(" ")
            isoform_id_full = transcript_info[1].strip("\"")
            isoform_id_partial = isoform_id_full.strip(";")
            final_isoform_id = isoform_id_partial.strip("\"")
            dict_value = [chr_num, strand, exon_start, exon_end]
            if final_isoform_id in gtf_dict:
                gtf_dict[final_isoform_id].append(dict_value)
            elif final_isoform_id not in gtf_dict:
                gtf_dict.update({final_isoform_id:[dict_value]})
    return gtf_dict


#combining isoseq classification and gtf data
#returns dictionary with key = gene_id and value == [gene_id, transcript_id, isoform, chr_num, strand, exon_list]*for each isoform that matches to this gene
def combine_isoseq():
    class_dict = read_isoseq_classification()
    gtf_dict = read_isoseq_gtf()
    combined_dict = {}
    for isoform in class_dict:
        if isoform in gtf_dict:
            single_class = class_dict[isoform]
            single_gtf = gtf_dict[isoform]
            gene_id = single_class[0]
            transcript_id = single_class[1]
            chr_num = single_class[3]
            strand = single_class[4]
            num_exons = single_class[5]
            exon_list = []
            for value in single_gtf:
                exon_start = int(value[2])
                exon_end = int(value[3])
                exon_positions = [exon_start, exon_end]
                exon_list.append(exon_positions)
            dict_value = [gene_id, transcript_id, isoform, chr_num, strand, exon_list]
            if gene_id in combined_dict:
                combined_dict[gene_id].append(dict_value)
            elif gene_id not in combined_dict:
                combined_dict.update({gene_id:[dict_value]})
    return combined_dict


#read ensembl gtf to create a dictionary similar to combine_isoseq
#returns dictionary key = gene id and value == [gene id, transcript id, chr num, strand, [exon start, exon end]]... repeated for each exon in a gene/transcript
def read_ensembl_gtf():
    ensembl_gtf_file = sys.argv[3]
    ensembl_id_dict = {}
    ensembl_exon_dict = {}
    ensembl_gtf_dict = {}
    with open(ensembl_gtf_file, 'r') as ensembl_gtf:
        for line in ensembl_gtf:
            new_line = line.split("\t")
            chr_num = new_line[0]
            feature = new_line[2]
            strand = new_line[6]
            exon_list = []
            if feature == "transcript":
                id_info = new_line[8].split(" ")
                gene_id_full = id_info[1].strip(";")
                gene_id = gene_id_full.strip("\"")
                transcript_id_full = id_info[5].strip(";")
                transcript_id = transcript_id_full.strip("\"")
                dict_value = [gene_id, transcript_id, chr_num, strand]
                if gene_id in ensembl_id_dict:
                    ensembl_id_dict[gene_id].append(dict_value)
                elif gene_id not in ensembl_id_dict:
                    ensembl_id_dict.update({gene_id:[dict_value]})
            elif feature == "exon":
                exon_start = int(new_line[3])
                exon_end = int(new_line[4])
                exon_positions = [transcript_id, exon_start, exon_end]
                if gene_id in ensembl_exon_dict:
                    ensembl_exon_dict[gene_id].append(exon_positions)
                elif gene_id not in ensembl_exon_dict:
                    ensembl_exon_dict.update({gene_id:[exon_positions]})
        for gene in ensembl_id_dict:
            if gene in ensembl_exon_dict:
                single_gene_id = ensembl_id_dict[gene]
                single_gene_exon = ensembl_exon_dict[gene]
                if len(single_gene_id) == 1:
                    single = single_gene_id[0]
                    if len(single_gene_exon) == 1:
                        single_exon = single_gene_exon[0]
                        if single[1] == single_exon[0]:
                            final_value = [gene, single[1], single[2], single[3], single_exon[1:3]]
                            ensembl_gtf_dict.update({gene:final_value})
                    elif len(single_gene_exon) > 1:
                        for value in single_gene_exon:
                            if single[1] == value[0]:
                                final_value = [gene, single[1], single[2], single[3], value[1:3]]
                                if gene in ensembl_gtf_dict:
                                    ensembl_gtf_dict[gene].append(final_value)
                                elif gene not in ensembl_gtf_dict:
                                    ensembl_gtf_dict.update({gene:[final_value]})
                elif len(single_gene_id) > 1:
                    for val in single_gene_id:
                        if len(single_gene_exon) == 1:
                            single_exon = single_gene_exon[0]
                            if val[1] == single_exon[0]:
                                final_value = [gene, val[1], val[2], val[3], single_exon[1:3]]
                                ensembl_gtf_dict.update({gene:final_value})
                        elif len(single_gene_exon) > 1:
                            for v in single_gene_exon:
                                if val[1] == v[0]:
                                    final_value = [gene, val[1], val[2], val[3], v[1:3]]
                                    if gene in ensembl_gtf_dict:
                                        ensembl_gtf_dict[gene].append(final_value)
                                    elif gene not in ensembl_gtf_dict:
                                        ensembl_gtf_dict.update({gene:[final_value]})
    return ensembl_gtf_dict


#combine exons into one dictionary value
#returns dictionary with key = gene id and value == [gene id, transcript id, chr num, strand, [exons]] where [exons] == [[exon 1 start, exon 1 end], [exon 2 start, exon 2 end], [exon 3 start, exon 3 end]...]
def reduce_ensembl_dict():
    ensembl_dict = read_ensembl_gtf()
    reduced_dict = {}
    final_dict = {}
    for gene in ensembl_dict:
        single_ensembl_exon = ensembl_dict[gene]
        if len(single_ensembl_exon) == 5 and isinstance(single_ensembl_exon[0], list) == False:
            reduced_dict.update({gene:single_ensembl_exon})
        else:
            transcript_dict = {}
            for exon in single_ensembl_exon:
                gene_id = exon[0]
                transcript_id = exon[1]
                chr_num = exon[2]
                strand = exon[3]
                single_exon = exon[4]
                if transcript_id in transcript_dict:
                    transcript_dict[transcript_id].append(single_exon)
                elif transcript_id not in transcript_dict:
                    transcript_dict.update({transcript_id:[single_exon]})
            for t in transcript_dict:
                final_value = [gene_id, t, chr_num, strand, transcript_dict[t]]
                if gene in reduced_dict:
                    reduced_dict[gene].append(final_value)
                elif gene not in reduced_dict:
                    reduced_dict.update({gene:[final_value]})
    for key in reduced_dict:
        single_key = reduced_dict[key]
        if len(single_key) == 5 and isinstance(single_key[0], list) == False:
            exon_list = single_key[4]
            if len(single_key[4]) == 2:
                continue
            else:
                final_dict.update({key:single_key})
        else:
            for value in single_key:
                exon_list = value[4]
                if len(exon_list) <= 2:
                    continue
                elif len(exon_list) > 2:
                    if key in final_dict:
                        final_dict[key].append(value)
                    elif key not in final_dict:
                        final_dict.update({key:[value]})
    return final_dict


#combine ensembl and isoseq exon dictionaries
#will need to remove ensembl genes/transcript that are matches to Isoseq genes/transcript (basically if gene and transcript are present in ensembl and isoseq remove the ensembl entry)
def combine_ensembl_isoseq():
    isoseq_list = []
    ensembl_list = []
    ensembl_dict = reduce_ensembl_dict()
    isoseq_dict = combine_isoseq()
    final_combined_dict = {}
    #loops through the ensembl dictionary
    for key in ensembl_dict:
        single_ensembl = ensembl_dict[key]
        #if there is only 1 ensembl entry
        if len(single_ensembl) == 1:
            single = single_ensembl[0]
            ens_gene_id = single[0]
            ens_transcript_id = single[1]
            ens_chr_num = single[2]
            ens_strand = single[3]
            ens_exons = single[4]
            #loop through isoseq dict using ensembl key that has only 1 entry
            if key in isoseq_dict:
                single_isoseq = isoseq_dict[key]
                #if isoseq only has 1 entry for the ensembl key with only 1 entry
                if len(single_isoseq) == 1:
                    final_isoseq_single = single_isoseq[0]
                    iso_gene_id = final_isoseq_single[0]
                    iso_transcript_id = final_isoseq_single[1]
                    iso_isoform_id = final_isoseq_single[2]
                    iso_chr_num = final_isoseq_single[3]
                    iso_strand = final_isoseq_single[4]
                    iso_exons = final_isoseq_single[5]
                    #if ensembl gene id and isoseq gene id match (this is to check matches between ensembl and isoseq to determine which to keep)
                    #if ensembl exons are the same as the isoseq exons (except first exon start (pos 0) and last exon end (pos len(exon_list)), only the isoseq exons will be kept
                    #if the exons are different by more than 2 exons that are not the first exon start and last exon end or if one of the exons in the middle have a different start or end, both sets of exons will be kept
                    if (ens_gene_id == iso_gene_id) and (ens_transcript_id == iso_transcript_id):
                        #examines positive strand first, do not have to change order on exons (as exons are listed numerically in ascending order)
                        if ens_strand == iso_strand and ens_strand == "+":
                            all_ens_exons = [item for sublist in ens_exons for item in sublist]
                            all_iso_exons = [item for sublist in iso_exons for item in sublist]
                            sorted_all_ens_exons = sorted(all_ens_exons, key=int, reverse = False)
                            sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = False)
                            #looks at exon lists that have the same number of exons
                            if len(all_ens_exons) == len(all_iso_exons):
                                pos_shared = 0
                                index_not_shared = []
                                for index, val in enumerate(all_ens_exons):
                                    if val == all_iso_exons[index]:
                                        pos_shared += 1
                                    else:
                                        index_not_shared.append(index)
                                #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                if len(index_not_shared) <= 2:
                                    if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                        ensembl_list.append(single)
                                        if key in final_combined_dict:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                    elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                        ensembl_list.append(single)
                                        if key in final_combined_dict:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                    else:
                                        if key in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict[key].append(single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict.update({key:[single]})
                                                final_combined_dict[key].append(final_isoseq_single)
                                #if differences are greater than 2 exons
                                else:
                                    if key in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict[key].append(single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                    elif key not in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict.update({key:[single]})
                                            final_combined_dict[key].append(final_isoseq_single)
                            #looks at exon pairs that differ in exon number
                            else:
                                if key in final_combined_dict:
                                    if single in ensembl_list:
                                        final_combined_dict[key].append(final_isoseq_single)
                                    elif single not in ensembl_list:
                                        ensembl_list.append(single)
                                        final_combined_dict[key].append(single)
                                        final_combined_dict[key].append(final_isoseq_single)
                                elif key not in final_combined_dict:
                                    if single in ensembl_list:
                                        final_combined_dict.update({key:[final_isoseq_single]})
                                    elif single not in ensembl_list:
                                        ensembl_list.append(single)
                                        final_combined_dict.update({key:[single]})
                                        final_combined_dict[key].append(final_isoseq_single)
                        #looks at exon pairs of - strand
                        elif ens_strand == iso_strand and ens_strand == "-":
                            sorted_all_ens_exons_1 = sorted(ens_exons, key = lambda x: x[0], reverse=True)
                            sorted_all_iso_exons_1 = sorted(iso_exons, key = lambda x: x[0], reverse=True)
                            all_ens_exons = [item for sublist in sorted_all_ens_exons_1 for item in sublist]
                            all_iso_exons = [item for sublist in sorted_all_iso_exons_1 for item in sublist]
                            sorted_all_ens_exons = sorted(all_ens_exons, key=int)
                            sorted_all_iso_exons = sorted(all_iso_exons, key=int)
                            if len(sorted_all_ens_exons) == len(all_iso_exons):
                                pos_shared = 0
                                index_not_shared = []
                                for index, v in enumerate(sorted_all_ens_exons):
                                    if v == sorted_all_iso_exons[index]:
                                        pos_shared += 1
                                    else:
                                        index_not_shared.append(index)
                                #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                if len(index_not_shared) <= 2:
                                    if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                        if single in ensembl_list:
                                            for check in final_combined_dict:
                                                single_check = final_combined_dict[check]
                                                if len(single_check) > 1:
                                                    for check_value in single_check:
                                                        if check_value == single:
                                                            final_combined_dict[key].remove(single)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                        if key in final_combined_dict:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                    elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                        if single in ensembl_list:
                                            for check in final_combined_dict:
                                                single_check = final_combined_dict[check]
                                                if len(single_check) > 1:
                                                    for check_value in single_check:
                                                        if check_value == single:
                                                            final_combined_dict[key].remove(single)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                        if key in final_combined_dict:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                    else:
                                        if key in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict[key].append(single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict.update({key:[single]})
                                                final_combined_dict[key].append(final_isoseq_single)
                                #if differences are greater than 2 exons
                                else:
                                    if key in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict[key].append(single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                    elif key not in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict.update({key:[single]})
                                            final_combined_dict[key].append(final_isoseq_single)
                            #looks at exon pairs that differ in exon number
                            else:
                                if key in final_combined_dict:
                                    if single in ensembl_list:
                                        final_combined_dict[key].append(final_isoseq_single)
                                    elif single not in ensembl_list:
                                        ensembl_list.append(single)
                                        final_combined_dict[key].append(single)
                                        final_combined_dict[key].append(final_isoseq_single)
                                elif key not in final_combined_dict:
                                    if single in ensembl_list:
                                        final_combined_dict.update({key:[final_isoseq_single]})
                                    elif single not in ensembl_list:
                                        ensembl_list.append(single)
                                        final_combined_dict.update({key:[single]})
                                        final_combined_dict[key].append(final_isoseq_single)
                    #examines transcripts that do not match in transcript id
                    elif ens_transcript_id != iso_transcript_id and iso_transcript_id != "novel":
                        if key in final_combined_dict:
                            if single in ensembl_list:
                                final_combined_dict[key].append(final_isoseq_single)
                            elif single not in ensembl_list:
                                ensembl_list.append(single)
                                final_combined_dict[key].append(single)
                                final_combined_dict[key].append(final_isoseq_single)
                        elif key not in final_combined_dict:
                            if single in ensembl_list:
                                final_combined_dict.update({key:[final_isoseq_single]})
                            elif single not in ensembl_list:
                                ensembl_list.append(single)
                                final_combined_dict.update({key:[single]})
                                final_combined_dict[key].append(final_isoseq_single)
                    elif ens_transcript_id != iso_transcript_id and iso_transcript_id == "novel":
                        if key in final_combined_dict:
                            if single in ensembl_list:
                                final_combined_dict[key].append(final_isoseq_single)
                            elif single not in ensembl_list:
                                ensembl_list.append(single)
                                final_combined_dict[key].append(single)
                                final_combined_dict[key].append(final_isoseq_single)
                        elif key not in final_combined_dict:
                            if single in ensembl_list:
                                final_combined_dict.update({key:[final_isoseq_single]})
                            elif single not in ensembl_list:
                                ensembl_list.append(single)
                                final_combined_dict.update({key:[single]})
                                final_combined_dict[key].append(final_isoseq_single)
                #if isoseq key has more than 1 entry
                else:
                    for isoform in single_isoseq:
                        iso_gene_id = isoform[0]
                        iso_transcript_id = isoform[1]
                        iso_isoform_id = isoform[2]
                        iso_chr_num = isoform[3]
                        iso_strand = isoform[4]
                        iso_exons = isoform[5]
                        #if ensembl gene id and isoseq gene id match (this is to check matches between ensembl and isoseq to determine which to keep)
                        #if ensembl exons are the same as the isoseq exons (except first exon start (pos 0) and last exon end (pos len(exon_list)), only the isoseq exons will be kept
                        #if the exons are different by more than 2 exons that are not the first exon start and last exon end or if one of the exons in the middle have a different start or end, both sets of exons will be kept
                        if ens_gene_id == iso_gene_id and ens_transcript_id == iso_transcript_id:
                            #examines positive strand first, do not have to change order on exons (as exons are listed numerically in ascending order)
                            if ens_strand == iso_strand and ens_strand == "+":
                                all_ens_exons = [item for sublist in ens_exons for item in sublist]
                                all_iso_exons = [item for sublist in iso_exons for item in sublist]
                                sorted_all_ens_exons = sorted(all_ens_exons, key=int, reverse = False)
                                sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = False)
                                #looks at exon lists that have the same number of exons
                                if len(sorted_all_ens_exons) == len(sorted_all_iso_exons):
                                    pos_shared = 0
                                    index_not_shared = []
                                    for index, val in enumerate(sorted_all_ens_exons):
                                        if val == all_iso_exons[index]:
                                            pos_shared += 1
                                        else:
                                            index_not_shared.append(index)
                                    #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                    if len(index_not_shared) <= 2:
                                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                            if single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == single:
                                                                final_combined_dict[key].remove(single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[isoform]})
                                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                            if single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == single:
                                                                final_combined_dict[key].remove(single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[isoform]})
                                        else:
                                            if key in final_combined_dict:
                                                if single in ensembl_list:
                                                    final_combined_dict[key].append(isoform)
                                                elif single not in ensembl_list:
                                                    ensembl_list.append(single)
                                                    final_combined_dict[key].append(single)
                                                    final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                if single in ensembl_list:
                                                    final_combined_dict.update({key:[isoform]})
                                                elif single not in ensembl_list:
                                                    ensembl_list.append(single)
                                                    final_combined_dict.update({key:[single]})
                                                    final_combined_dict[key].append(isoform)
                                    #if differences are greater than 2 exons
                                    else:
                                        if key in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict[key].append(isoform)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict[key].append(single)
                                                final_combined_dict[key].append(isoform)
                                        elif key not in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict.update({key:[isoform]})
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict.update({key:[single]})
                                                final_combined_dict[key].append(isoform)
                                #looks at exon pairs that differ in exon number
                                else:
                                    if key in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict[key].append(isoform)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict[key].append(single)
                                            final_combined_dict[key].append(isoform)
                                    elif key not in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict.update({key:[isoform]})
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict.update({key:[single]})
                                            final_combined_dict[key].append(isoform)
                            #looks at exon pairs of - strand
                            elif ens_strand == iso_strand and ens_strand == "-":
                                sorted_all_ens_exons_1 = sorted(ens_exons, key = lambda x: x[0], reverse=True)
                                sorted_all_iso_exons_1 = sorted(iso_exons, key = lambda x: x[0], reverse=True)
                                all_ens_exons = [item for sublist in sorted_all_ens_exons_1 for item in sublist]
                                all_iso_exons = [item for sublist in sorted_all_iso_exons_1 for item in sublist]
                                sorted_all_ens_exons = sorted(all_ens_exons, key=int)
                                sorted_all_iso_exons = sorted(all_iso_exons, key=int)
                                if len(sorted_all_ens_exons) == len(all_iso_exons):
                                    pos_shared = 0
                                    index_not_shared = []
                                    for index, v in enumerate(sorted_all_ens_exons):
                                        if v == sorted_all_iso_exons[index]:
                                            pos_shared += 1
                                        else:
                                            index_not_shared.append(index)
                                    #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                    if len(index_not_shared) <= 2:
                                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                            if single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == single:
                                                                final_combined_dict[key].remove(single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[isoform]})
                                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                            if single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == single:
                                                                final_combined_dict[key].remove(single)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[isoform]})
                                        else:
                                            if key in final_combined_dict:
                                                if single in ensembl_list:
                                                    final_combined_dict[key].append(isoform)
                                                elif single not in ensembl_list:
                                                    ensembl_list.append(single)
                                                    final_combined_dict[key].append(single)
                                                    final_combined_dict[key].append(isoform)
                                            elif key not in final_combined_dict:
                                                if single in ensembl_list:
                                                    final_combined_dict.update({key:[isoform]})
                                                elif single not in ensembl_list:
                                                    ensembl_list.append(single)
                                                    final_combined_dict.update({key:[single]})
                                                    final_combined_dict[key].append(isoform)
                                    #if differences are greater than 2 exons
                                    else:
                                        if key in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict[key].append(isoform)
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict[key].append(single)
                                                final_combined_dict[key].append(isoform)
                                        elif key not in final_combined_dict:
                                            if single in ensembl_list:
                                                final_combined_dict.update({key:[isoform]})
                                            elif single not in ensembl_list:
                                                ensembl_list.append(single)
                                                final_combined_dict.update({key:[single]})
                                                final_combined_dict[key].append(isoform)
                                #looks at exon pairs that differ in exon number
                                else:
                                    if key in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict[key].append(isoform)
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict[key].append(single)
                                            final_combined_dict[key].append(isoform)
                                    elif key not in final_combined_dict:
                                        if single in ensembl_list:
                                            final_combined_dict.update({key:[isoform]})
                                        elif single not in ensembl_list:
                                            ensembl_list.append(single)
                                            final_combined_dict.update({key:[single]})
                                            final_combined_dict[key].append(isoform)
                        #examines transcripts that do not match in transcript id
                        elif ens_transcript_id != iso_transcript_id and iso_transcript_id != "novel":
                            if key in final_combined_dict:
                                if single in ensembl_list:
                                    final_combined_dict[key].append(isoform)
                                elif single not in ensembl_list:
                                    ensembl_list.append(single)
                                    final_combined_dict[key].append(single)
                                    final_combined_dict[key].append(isoform)
                            elif key not in final_combined_dict:
                                if single in ensembl_list:
                                    final_combined_dict.update({key:[isoform]})
                                elif single not in ensembl_list:
                                    ensembl_list.append(single)
                                    final_combined_dict.update({key:[single]})
                                    final_combined_dict[key].append(isoform)
                        elif ens_transcript_id != iso_transcript_id and iso_transcript_id == "novel":
                            if key in final_combined_dict:
                                if single in ensembl_list:
                                    final_combined_dict[key].append(isoform)
                                elif single not in ensembl_list:
                                    ensembl_list.append(single)
                                    final_combined_dict[key].append(single)
                                    final_combined_dict[key].append(isoform)
                            elif key not in final_combined_dict:
                                if single in ensembl_list:
                                    final_combined_dict.update({key:[isoform]})
                                elif single not in ensembl_list:
                                    ensembl_list.append(single)
                                    final_combined_dict.update({key:[single]})
                                    final_combined_dict[key].append(isoform)
            #if the ensembl gene is not in isoseq
            #else:
            #    final_combined_dict.update({key:[single_ensembl]})
        #if ensembl gene has multiple entries
        else:
            for ens_single in single_ensembl:
                ens_gene_id = ens_single[0]
                ens_transcript_id = ens_single[1]
                ens_chr_num = ens_single[2]
                ens_strand = ens_single[3]
                ens_exons = ens_single[4]
                #if the ensembl gene is in the isoseq dictionary
                if key in isoseq_dict:
                    single_isoseq = isoseq_dict[key]
                    #if the isoseq entry only has one entry
                    if len(single_isoseq) == 1:
                        final_isoseq_single = single_isoseq[0]
                        iso_gene_id = final_isoseq_single[0]
                        iso_transcript_id = final_isoseq_single[1]
                        iso_isoform_id = final_isoseq_single[2]
                        iso_chr_num = final_isoseq_single[3]
                        iso_strand = final_isoseq_single[4]
                        iso_exons = final_isoseq_single[5]
                        #if ensembl gene id and isoseq gene id match (this is to check matches between ensembl and isoseq to determine which to keep)
                        #if ensembl exons are the same as the isoseq exons (except first exon start (pos 0) and last exon end (pos len(exon_list)), only the isoseq exons will be kept
                        #if the exons are different by more than 2 exons that are not the first exon start and last exon end or if one of the exons in the middle have a different start or end, both sets of exons will be kept
                        if ens_gene_id == iso_gene_id and ens_transcript_id == iso_transcript_id:
                            #examines positive strand first, do not have to change order on exons (as exons are listed numerically in ascending order)
                            if ens_strand == iso_strand and ens_strand == "+":
                                all_ens_exons = [item for sublist in ens_exons for item in sublist]
                                all_iso_exons = [item for sublist in iso_exons for item in sublist]
                                sorted_all_ens_exons = sorted(all_ens_exons, key=int, reverse = False)
                                sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = False)
                                #looks at exon lists that have the same number of exons
                                if len(sorted_all_ens_exons) == len(sorted_all_iso_exons):
                                    pos_shared = 0
                                    index_not_shared = []
                                    for index, val in enumerate(sorted_all_ens_exons):
                                        if val == all_iso_exons[index]:
                                            pos_shared += 1
                                        else:
                                            index_not_shared.append(index)
                                    #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                    if len(index_not_shared) <= 2:
                                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                            isoseq_list.append(final_isoseq_single)
                                            if ens_single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == ens_single:
                                                                final_combined_dict[key].remove(ens_single)
                                            elif ens_single not in ensembl_list:
                                                ensembl_list.append(ens_single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                            isoseq_list.append(final_isoseq_single)
                                            if ens_single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == ens_single:
                                                                final_combined_dict[key].remove(ens_single)
                                            elif ens_single not in ensembl_list:
                                                ensembl_list.append(ens_single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                        else:
                                            if key in final_combined_dict:
                                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict[key].append(final_isoseq_single)
                                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict[key].append(ens_single)
                                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict[key].append(ens_single)
                                                    final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict.update({key:[final_isoseq_single]})
                                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                    final_combined_dict[key].append(final_isoseq_single)
                                    #if differences are greater than 2 exons
                                    else:
                                        if key in final_combined_dict:
                                            if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict[key].append(ens_single)
                                            elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict[key].append(ens_single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                            elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict.update({key:[ens_single]})
                                            elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict.update({key:[ens_single]})
                                                final_combined_dict[key].append(final_isoseq_single)
                                #looks at exon pairs that differ in exon number
                                else:
                                    if key in final_combined_dict:
                                        if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                            continue
                                        elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            final_combined_dict[key].append(ens_single)
                                        elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict[key].append(ens_single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                    elif key not in final_combined_dict:
                                        if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                            continue
                                        elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                        elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            final_combined_dict.update({key:[ens_single]})
                                        elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict.update({key:[ens_single]})
                                            final_combined_dict[key].append(final_isoseq_single)
                            #looks at exon pairs of - strand
                            elif ens_strand == iso_strand and ens_strand == "-":
                                sorted_all_ens_exons_1 = sorted(ens_exons, key = lambda x: x[0], reverse=True)
                                sorted_all_iso_exons_1 = sorted(iso_exons, key = lambda x: x[0], reverse=True)
                                all_ens_exons = [item for sublist in sorted_all_ens_exons_1 for item in sublist]
                                all_iso_exons = [item for sublist in sorted_all_iso_exons_1 for item in sublist]
                                sorted_all_ens_exons = sorted(all_ens_exons, key=int)
                                sorted_all_iso_exons = sorted(all_iso_exons, key=int)
                                if len(sorted_all_ens_exons) == len(all_iso_exons):
                                    pos_shared = 0
                                    index_not_shared = []
                                    for index, v in enumerate(sorted_all_ens_exons):
                                        if v == sorted_all_iso_exons[index]:
                                            pos_shared += 1
                                        else:
                                            index_not_shared.append(index)
                                    #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                    if len(index_not_shared) <= 2:
                                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                            isoseq_list.append(final_isoseq_single)
                                            if ens_single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == ens_single:
                                                                final_combined_dict[key].remove(ens_single)
                                            elif ens_single not in ensembl_list:
                                                ensembl_list.append(ens_single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                            isoseq_list.append(final_isoseq_single)
                                            if ens_single in ensembl_list:
                                                for check in final_combined_dict:
                                                    single_check = final_combined_dict[check]
                                                    if len(single_check) > 1:
                                                        for check_value in single_check:
                                                            if check_value == ens_single:
                                                                final_combined_dict[key].remove(ens_single)
                                            elif ens_single not in ensembl_list:
                                                ensembl_list.append(ens_single)
                                            if key in final_combined_dict:
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                        else:
                                            if key in final_combined_dict:
                                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict[key].append(final_isoseq_single)
                                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict[key].append(ens_single)
                                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict[key].append(ens_single)
                                                    final_combined_dict[key].append(final_isoseq_single)
                                            elif key not in final_combined_dict:
                                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict.update({key:[final_isoseq_single]})
                                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(final_isoseq_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                    final_combined_dict[key].append(final_isoseq_single)
                                    #if differences are greater than 2 exons
                                    else:
                                        if key in final_combined_dict:
                                            if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                            elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict[key].append(ens_single)
                                            elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict[key].append(ens_single)
                                                final_combined_dict[key].append(final_isoseq_single)
                                        elif key not in final_combined_dict:
                                            if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict.update({key:[final_isoseq_single]})
                                            elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict.update({key:[ens_single]})
                                            elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(final_isoseq_single)
                                                final_combined_dict.update({key:[ens_single]})
                                                final_combined_dict[key].append(final_isoseq_single)
                                #looks at exon pairs that differ in exon number
                                else:
                                    if key in final_combined_dict:
                                        if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                            continue
                                        elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                        elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            final_combined_dict[key].append(ens_single)
                                        elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict[key].append(ens_single)
                                            final_combined_dict[key].append(final_isoseq_single)
                                    elif key not in final_combined_dict:
                                        if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                            continue
                                        elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict.update({key:[final_isoseq_single]})
                                        elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            final_combined_dict.update({key:[ens_single]})
                                        elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                            ensembl_list.append(ens_single)
                                            isoseq_list.append(final_isoseq_single)
                                            final_combined_dict.update({key:[ens_single]})
                                            final_combined_dict[key].append(final_isoseq_single)
                        #examines transcripts that do not match in transcript id
                        elif ens_transcript_id != iso_transcript_id and iso_transcript_id != "novel":
                            if key in final_combined_dict:
                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                    continue
                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict[key].append(final_isoseq_single)
                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    final_combined_dict[key].append(ens_single)
                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict[key].append(ens_single)
                                    final_combined_dict[key].append(final_isoseq_single)
                            elif key not in final_combined_dict:
                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                    continue
                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict.update({key:[final_isoseq_single]})
                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    final_combined_dict.update({key:[ens_single]})
                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict.update({key:[ens_single]})
                                    final_combined_dict[key].append(final_isoseq_single)
                        elif ens_transcript_id != iso_transcript_id and iso_transcript_id == "novel":
                            if key in final_combined_dict:
                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                    continue
                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict[key].append(final_isoseq_single)
                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    final_combined_dict[key].append(ens_single)
                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict[key].append(ens_single)
                                    final_combined_dict[key].append(final_isoseq_single)
                            elif key not in final_combined_dict:
                                if ens_single in ensembl_list and final_isoseq_single in isoseq_list:
                                    continue
                                elif ens_single in ensembl_list and final_isoseq_single not in isoseq_list:
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict.update({key:[final_isoseq_single]})
                                elif ens_single not in ensembl_list and final_isoseq_single in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    final_combined_dict.update({key:[ens_single]})
                                elif ens_single not in ensembl_list and final_isoseq_single not in isoseq_list:
                                    ensembl_list.append(ens_single)
                                    isoseq_list.append(final_isoseq_single)
                                    final_combined_dict.update({key:[ens_single]})
                                    final_combined_dict[key].append(final_isoseq_single)
                    #if isoseq entry has more than 1 value
                    else:
                        for iso in single_isoseq:
                            iso_gene_id = iso[0]
                            iso_transcript_id = iso[1]
                            iso_isoform_id = iso[2]
                            iso_chr_num = iso[3]
                            iso_strand = iso[4]
                            iso_exons = iso[5]
                            #if ensembl gene id and isoseq gene id match (this is to check matches between ensembl and isoseq to determine which to keep)
                            #if ensembl exons are the same as the isoseq exons (except first exon start (pos 0) and last exon end (pos len(exon_list)), only the isoseq exons will be kept
                            #if the exons are different by more than 2 exons that are not the first exon start and last exon end or if one of the exons in the middle have a different start or end, both sets of exons will be kept
                            if ens_gene_id == iso_gene_id and ens_transcript_id == iso_transcript_id:
                                #examines positive strand first, do not have to change order on exons (as exons are listed numerically in ascending order)
                                if ens_strand == iso_strand and ens_strand == "+":
                                    all_ens_exons = [item for sublist in ens_exons for item in sublist]
                                    all_iso_exons = [item for sublist in iso_exons for item in sublist]
                                    sorted_all_ens_exons = sorted(all_ens_exons, key=int, reverse = False)
                                    sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = False)
                                    #looks at exon lists that have the same number of exons
                                    if len(sorted_all_ens_exons) == len(sorted_all_iso_exons):
                                        pos_shared = 0
                                        index_not_shared = []
                                        for index, val in enumerate(sorted_all_ens_exons):
                                            if val == all_iso_exons[index]:
                                                pos_shared += 1
                                            else:
                                                index_not_shared.append(index)
                                        #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                        if len(index_not_shared) <= 2:
                                            if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                                if ens_single in ensembl_list:
                                                    for check in final_combined_dict:
                                                        single_check = final_combined_dict[check]
                                                        if len(single_check) > 1:
                                                            for check_value in single_check:
                                                                if check_value == ens_single:
                                                                    final_combined_dict[key].remove(ens_single)
                                                elif ens_single not in ensembl_list:
                                                    ensembl_list.append(ens_single)
                                                if key in final_combined_dict and iso in isoseq_list:
                                                    continue
                                                elif key in final_combined_dict and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict and iso in isoseq_list:
                                                    continue
                                                elif key not in final_combined_dict and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict.update({key:[iso]})
                                            elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                                if ens_single in ensembl_list:
                                                    for check in final_combined_dict:
                                                        single_check = final_combined_dict[check]
                                                        if len(single_check) > 1:
                                                            for check_value in single_check:
                                                                if check_value == ens_single:
                                                                    final_combined_dict[key].remove(ens_single)
                                                elif ens_single not in ensembl_list:
                                                    ensembl_list.append(ens_single)
                                                if key in final_combined_dict and iso in isoseq_list:
                                                    continue
                                                elif key in final_combined_dict and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict and iso in isoseq_list:
                                                    continue
                                                elif key not in final_combined_dict and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict.update({key:[iso]})
                                            else:
                                                if key in final_combined_dict:
                                                    if ens_single in ensembl_list and iso in isoseq_list:
                                                        continue
                                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                                        isoseq_list.append(iso)
                                                        final_combined_dict[key].append(iso)
                                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        final_combined_dict[key].append(ens_single)
                                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        isoseq_list.append(iso)
                                                        final_combined_dict[key].append(ens_single)
                                                        final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict:
                                                    if ens_single in ensembl_list and iso in isoseq_list:
                                                        continue
                                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                                        isoseq_list.append(iso)
                                                        final_combined_dict.update({key:[iso]})
                                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        final_combined_dict.update({key:[ens_single]})
                                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        isoseq_list.append(iso)
                                                        final_combined_dict.update({key:[ens_single]})
                                                        final_combined_dict[key].append(iso)
                                        #if differences are greater than 2 exons
                                        else:
                                            if key in final_combined_dict:
                                                if ens_single in ensembl_list and iso in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(iso)
                                                elif ens_single not in ensembl_list and iso in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict[key].append(ens_single)
                                                elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(ens_single)
                                                    final_combined_dict[key].append(iso)
                                            elif key not in final_combined_dict:
                                                if ens_single in ensembl_list and iso in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict.update({key:[iso]})
                                                elif ens_single not in ensembl_list and iso in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(iso)
                                                    final_combined_dict.update({key:[ens_single]})
                                                    final_combined_dict[key].append(iso)
                                    #looks at exon pairs that differ in exon number
                                    else:
                                        if key in final_combined_dict:
                                            if ens_single in ensembl_list and iso in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and iso not in isoseq_list:
                                                isoseq_list.append(iso)
                                                final_combined_dict[key].append(iso)
                                            elif ens_single not in ensembl_list and iso in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict[key].append(ens_single)
                                            elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(iso)
                                                final_combined_dict[key].append(ens_single)
                                                final_combined_dict[key].append(iso)
                                        elif key not in final_combined_dict:
                                            if ens_single in ensembl_list and iso in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and iso not in isoseq_list:
                                                isoseq_list.append(iso)
                                                final_combined_dict.update({key:[iso]})
                                            elif ens_single not in ensembl_list and iso in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict.update({key:[ens_single]})
                                            elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(iso)
                                                final_combined_dict.update({key:[ens_single]})
                                                final_combined_dict[key].append(iso)
                                #looks at exon pairs of - strand
                                elif ens_strand == iso_strand and ens_strand == "-":
                                    sorted_all_ens_exons_1 = sorted(ens_exons, key = lambda x: x[0], reverse=True)
                                    sorted_all_iso_exons_1 = sorted(iso_exons, key = lambda x: x[0], reverse=True)
                                    all_ens_exons = [item for sublist in sorted_all_ens_exons_1 for item in sublist]
                                    all_iso_exons = [item for sublist in sorted_all_iso_exons_1 for item in sublist]
                                    sorted_all_ens_exons = sorted(all_ens_exons, key=int)
                                    sorted_all_iso_exons = sorted(all_iso_exons, key=int)
                                    if len(sorted_all_ens_exons) == len(all_iso_exons):
                                        pos_shared = 0
                                        index_not_shared = []
                                        for index, v in enumerate(sorted_all_ens_exons):
                                            if v == sorted_all_iso_exons[index]:
                                                pos_shared += 1
                                            else:
                                                index_not_shared.append(index)
                                        #examines pairs (ensembl and isoseq exons) that differ in less than 2 exons
                                        if len(index_not_shared) <= 2:
                                            if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(all_ens_exons)-1):
                                                isoseq_list.append(iso)
                                                if ens_single in ensembl_list:
                                                    for check in final_combined_dict:
                                                        single_check = final_combined_dict[check]
                                                        if len(single_check) > 1:
                                                            for check_value in single_check:
                                                                if check_value == ens_single:
                                                                    final_combined_dict[key].remove(ens_single)
                                                elif ens_single not in ensembl_list:
                                                    ensembl_list.append(ens_single)
                                                if key in final_combined_dict:
                                                    final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict:
                                                    final_combined_dict.update({key:[iso]})
                                            elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(all_ens_exons)-1:
                                                isoseq_list.append(iso)
                                                if ens_single in ensembl_list:
                                                    for check in final_combined_dict:
                                                        single_check = final_combined_dict[check]
                                                        if len(single_check) > 1:
                                                            for check_value in single_check:
                                                                if check_value == ens_single:
                                                                    final_combined_dict[key].remove(ens_single)
                                                elif ens_single not in ensembl_list:
                                                    ensembl_list.append(ens_single)
                                                if key in final_combined_dict:
                                                    final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict:
                                                    final_combined_dict.update({key:[iso]})
                                            else:
                                                if key in final_combined_dict:
                                                    if ens_single in ensembl_list and iso in isoseq_list:
                                                        continue
                                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                                        isoseq_list.append(iso)
                                                        final_combined_dict[key].append(iso)
                                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        final_combined_dict[key].append(ens_single)
                                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        isoseq_list.append(iso)
                                                        final_combined_dict[key].append(ens_single)
                                                        final_combined_dict[key].append(iso)
                                                elif key not in final_combined_dict:
                                                    if ens_single in ensembl_list and iso in isoseq_list:
                                                        continue
                                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                                        isoseq_list.append(iso)
                                                        final_combined_dict.update({key:[iso]})
                                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        final_combined_dict.update({key:[ens_single]})
                                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                        ensembl_list.append(ens_single)
                                                        isoseq_list.append(iso)
                                                        final_combined_dict.update({key:[ens_single]})
                                                        final_combined_dict[key].append(iso)
                                        #if differences are greater than 2 exons
                                        else:
                                            if key in final_combined_dict:
                                                if ens_single in ensembl_list and iso in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and iso not in isoseq_list:
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(iso)
                                                elif ens_single not in ensembl_list and iso in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict[key].append(ens_single)
                                                elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(iso)
                                                    final_combined_dict[key].append(ens_single)
                                                    final_combined_dict[key].append(iso)
                                            elif key not in final_combined_dict:
                                                if ens_single in ensembl_list and iso in isoseq_list:
                                                    continue
                                                elif ens_single in ensembl_list and iso not in isoseq_list:
                                                    isoseq_list.append()
                                                    final_combined_dict.update({key:[iso]})
                                                elif ens_single not in ensembl_list and iso in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    final_combined_dict.update({key:[ens_single]})
                                                elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                    ensembl_list.append(ens_single)
                                                    isoseq_list.append(iso)
                                                    final_combined_dict.update({key:[ens_single]})
                                                    final_combined_dict[key].append(iso)
                                    #looks at exon pairs that differ in exon number
                                    else:
                                        if key in final_combined_dict:
                                            if ens_single in ensembl_list and iso in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and iso not in isoseq_list:
                                                isoseq_list.append(iso)
                                                final_combined_dict[key].append(iso)
                                            elif ens_single not in ensembl_list and iso in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict[key].append(ens_single)
                                            elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(iso)
                                                final_combined_dict[key].append(ens_single)
                                                final_combined_dict[key].append(iso)
                                        elif key not in final_combined_dict:
                                            if ens_single in ensembl_list and iso in isoseq_list:
                                                continue
                                            elif ens_single in ensembl_list and iso not in isoseq_list:
                                                isoseq_list.append(iso)
                                                final_combined_dict.update({key:[iso]})
                                            elif ens_single not in ensembl_list and iso in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                final_combined_dict.update({key:[ens_single]})
                                            elif ens_single not in ensembl_list and iso not in isoseq_list:
                                                ensembl_list.append(ens_single)
                                                isoseq_list.append(iso)
                                                final_combined_dict.update({key:[ens_single]})
                                                final_combined_dict[key].append(iso)
                            #examines transcripts that do not match in transcript id
                            elif ens_transcript_id != iso_transcript_id and iso_transcript_id != "novel":
                                if key in final_combined_dict:
                                    if ens_single in ensembl_list and iso in isoseq_list:
                                        continue
                                    #elif ens_single in ensembl_list and iso not in isoseq_list:
                                    #    isoseq_list.append(iso)
                                    #    final_combined_dict[key].append(iso)
                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                        ensembl_list.append(ens_single)
                                        final_combined_dict[key].append(ens_single)
                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                        ensembl_list.append(ens_single)
                                        #isoseq_list.append(iso)
                                        final_combined_dict[key].append(ens_single)
                                        #final_combined_dict[key].append(iso)
                                elif key not in final_combined_dict:
                                    if ens_single in ensembl_list and iso in isoseq_list:
                                        continue
                                    #elif ens_single in ensembl_list and iso not in isoseq_list:
                                    #    isoseq_list.append(iso)
                                    #    final_combined_dict.update({key:[iso]})
                                    elif ens_single not in ensembl_list and iso in isoseq_list:
                                        ensembl_list.append(ens_single)
                                        final_combined_dict.update({key:[ens_single]})
                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                        ensembl_list.append(ens_single)
                                        #isoseq_list.append(iso)
                                        final_combined_dict.update({key:[ens_single]})
                                        #final_combined_dict[key].append(iso)
                            elif ens_transcript_id != iso_transcript_id and iso_transcript_id == "novel":
                                if key in final_combined_dict:
                                    if ens_single in ensembl_list and iso in isoseq_list:
                                        continue
                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                        isoseq_list.append(iso)
                                        final_combined_dict[key].append(iso)
                                    #elif ens_single not in ensembl_list and iso in isoseq_list:
                                    #    ensembl_list.append(ens_single)
                                    #    final_combined_dict[key].append(ens_single)
                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                        #ensembl_list.append(ens_single)
                                        isoseq_list.append(iso)
                                        #final_combined_dict[key].append(ens_single)
                                        final_combined_dict[key].append(iso)
                                elif key not in final_combined_dict:
                                    if ens_single in ensembl_list and iso in isoseq_list:
                                        continue
                                    elif ens_single in ensembl_list and iso not in isoseq_list:
                                        isoseq_list.append(iso)
                                        final_combined_dict.update({key:[iso]})
                                    #elif ens_single not in ensembl_list and iso in isoseq_list:
                                    #    ensembl_list.append(ens_single)
                                    #    final_combined_dict.update({key:[ens_single]})
                                    elif ens_single not in ensembl_list and iso not in isoseq_list:
                                        #ensembl_list.append(ens_single)
                                        isoseq_list.append(iso)
                                        #final_combined_dict.update({key:[ens_single]})
                                        final_combined_dict.update({key:[iso]})
                #if the ensembl gene is not in isoseq
                #else:
                #    final_combined_dict.update({key:[single_ensembl]})
    for key2 in isoseq_dict:
        if key2 not in ensembl_dict:
            single_key2 = isoseq_dict[key2]
            final_combined_dict.update({key2:single_key2})
    return final_combined_dict



#need to run through dictionary and reduce any isoforms that are the same except for first exon start and last exon end will take the most upstream start site and the farthest downstream end site
#also will make sure here aren't too many [] around dictionary values; want single values to have [[single value]] so that the length == 1.
def reduce_combined_dict():
    combined_dict = combine_ensembl_isoseq()
    final_dictionary = {}
    gene_list = []
    duplicates_list = []
    non_duplicates_list = []
    final_duplicates = []
    final_non_duplicates = []
    for gene in combined_dict:
        single_gene = combined_dict[gene]
        #if there is only 1 isoform for a gene, this adds the gene to the final dictionary without extra []; though some isoforms might have 2 sets of [] which will need to be removed in the next function
        if len(single_gene) == 1:
            final_dictionary.update({gene:single_gene[0]})
        #if there is more than 1 isoform for a gene, loop through each isoform
        else:
            for value in single_gene:
                if len(value) == 1:
                    final_dictionary.update({gene:value})
                elif len(value) == 5 or len(value) == 6:
                    gene_list.append(value)
                else:
                    for val in value:
                        gene_list.append(val)
    x = 0
    while x < len(gene_list):
        test_isoform = gene_list[x]
        if len(test_isoform) == 5:
            test_gene_id = test_isoform[0]
            test_transcript_id = test_isoform[1]
            test_chr_num = test_isoform[2]
            test_strand = test_isoform[3]
            test_exons = test_isoform[4]
        elif len(test_isoform) == 6:
            test_gene_id = test_isoform[0]
            test_transcript_id = test_isoform[1]
            test_isoform_id = test_isoform[2]
            test_chr_num = test_isoform[3]
            test_strand = test_isoform[4]
            test_exons = test_isoform[5]
        for index, isoform in enumerate(gene_list):
            if index == x:
                continue
            else:
                isoform_exons = isoform[len(isoform)-1]
                if len(test_exons) == len(isoform_exons):
                    if test_strand == "+":
                        all_test_exons = [item for sublist in test_exons for item in sublist]
                        all_iso_exons = [item for sublist in isoform_exons for item in sublist]
                        sorted_all_test_exons = sorted(all_test_exons, key=int, reverse = False)
                        sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = False)
                        pos_shared = 0
                        index_not_shared = []
                        for ind, v in enumerate(sorted_all_test_exons):
                            if v == sorted_all_iso_exons[ind]:
                                pos_shared += 1
                            else:
                                index_not_shared.append(ind)
                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(sorted_all_test_exons)-1):
                            if test_isoform in duplicates_list:
                                continue
                            elif test_isoform not in duplicates_list:
                                duplicates_list.append(test_isoform)
                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(sorted_all_test_exons)-1:
                            if test_isoform in duplicates_list:
                                continue
                            elif test_isoform not in duplicates_list:
                                duplicates_list.append(test_isoform)
                    elif test_strand == "-":
                        all_test_exons = [item for sublist in test_exons for item in sublist]
                        all_iso_exons = [item for sublist in isoform_exons for item in sublist]
                        sorted_all_test_exons = sorted(all_test_exons, key=int, reverse = True)
                        sorted_all_iso_exons = sorted(all_iso_exons, key=int, reverse = True)
                        pos_shared = 0
                        index_not_shared = []
                        for ind, v in enumerate(sorted_all_test_exons):
                            if v == sorted_all_iso_exons[ind]:
                                pos_shared += 1
                            else:
                                index_not_shared.append(ind)
                        if len(index_not_shared) == 1 and (index_not_shared[0] == 0 or index_not_shared[0] == len(sorted_all_test_exons)-1):
                            if test_isoform in duplicates_list:
                                continue
                            elif test_isoform not in duplicates_list:
                                duplicates_list.append(test_isoform)
                        elif len(index_not_shared) == 2 and index_not_shared[0] == 0 and index_not_shared[1] == len(sorted_all_test_exons)-1:
                            if test_isoform in duplicates_list:
                                continue
                            elif test_isoform not in duplicates_list:
                                duplicates_list.append(test_isoform)
        x += 1
    dups_to_remove = []
    for key in combined_dict:
        single_gene = combined_dict[key]
        if len(single_gene) > 1:
            for val in single_gene:
                if len(val) == 5 or len(val) == 6:
                    exons_list_num = len(val[len(val)-1])
                    if val in duplicates_list and exons_list_num in dups_to_remove:
                        continue
                    elif val in duplicates_list and exons_list_num not in dups_to_remove:
                        dups_to_remove.append(exons_list_num)
                        if key in final_dictionary:
                            final_dictionary[key].append(val)
                        elif key not in final_dictionary:
                            final_dictionary.update({key:[val]})
                    elif val not in duplicates_list:
                        if key in final_dictionary:
                            final_dictionary[key].append(val)
                        elif key not in final_dictionary:
                            final_dictionary.update({key:[val]})
                else:
                    for v in val:
                        exons_list_num = len(v[len(v)-1])
                        if v in duplicates_list and exons_list_num in dups_to_remove:
                            continue
                        elif v in duplicates_list and exons_list_num not in dups_to_remove:
                            dups_to_remove.append(exons_list_num)
                            if key in final_dictionary:
                                final_dictionary[key].append(v)
                            elif key not in final_dictionary:
                                final_dictionary.update({key:[v]})
                        elif v not in duplicates_list:
                            if key in final_dictionary:
                                final_dictionary[key].append(v)
                            elif key not in final_dictionary:
                                final_dictionary.update({key:[v]})
        dups_to_remove = []
    return final_dictionary

#converting reduced combined dictionary so all + strand genes increase in exon positions and all - strand genes decrease in exon positions
def convert_direction_combined_dict():
    combined_dict = reduce_combined_dict()
    converted_dict = {}
    for key in combined_dict:
        single_key = combined_dict[key]
        #if there is one isoform per gene and the isoform is encompassed by [[]]
        if len(single_key) == 1:
            single = single_key[0]
            if len(single) == 5:
                strand = single[3]
                exons = single[4]
                single_key_info = single_key[0][0:len(single)-1]
            elif len(single_key) == 6:
                strand = single_key[4]
                exons = single_key[5]
                single_key_info = single_key[0][0:len(single)-1]
            #if strand == +, the exons should increase in position
            if strand == "+":
                new_exons = sorted(exons, key=lambda x: x[0], reverse = False)
                final = single_key_info + [new_exons]
                converted_dict.update({key:final})
            #if strand == -, the exons should decrease in position
            elif strand == "-":
                final_exons = []
                #this first sorts all exon pairs (start, end) in descending order
                new_exons = sorted(exons, key=lambda x: x[0], reverse = True)
                #because for - strand genes, the second position in a pair is actually the exon start, need to flip all exon pairs
                for exon in new_exons:
                    exon.reverse()
                    final_exons.append(exon)
                final = single_key_info + [final_exons]
                converted_dict.update({key:final})
        #if there is only one isoform for a gene that is not encompassed by [[]]
        elif (len(single_key) == 6 or len(single_key) == 5) and isinstance(single_key[0], list) == False:
            if len(single_key) == 5:
                strand = single_key[3]
                exons = single_key[4]
                single_key_info = single_key[0:len(single_key)-1]
            elif len(single_key) == 6:
                strand = single_key[4]
                exons = single_key[5]
                single_key_info = single_key[0:len(single_key)-1]
            #if strand == +, the exons should increase in position
            if strand == "+":
                new_exons = sorted(exons, key=lambda x: x[0], reverse = False)
                final = single_key_info + [new_exons]
                converted_dict.update({key:final})
            #if strand == -, the exons should decrease in position
            elif strand == "-":
                final_exons = []
                #this first sorts all exon pairs (start, end) in descending order
                new_exons = sorted(exons, key=lambda x: x[0], reverse = True)
                #because for - strand genes, the second position in a pair is actually the exon start, need to flip all exon pairs
                for exon in new_exons:
                    exon.reverse()
                    final_exons.append(exon)
                final = single_key_info + [final_exons]
                converted_dict.update({key:final})
        #if there are multiple isoforms for a gene
        else:
            #reads through each isoform one at a time
            for single in single_key:
                if len(single) == 5:
                    strand = single[3]
                    exons = single[4]
                    single_info = single[0:len(single)-1]
                elif len(single) == 6:
                    strand = single[4]
                    exons = single[5]
                    single_info = single[0:len(single)-1]
                #if strand == +, the exons should increase in position
                if strand == "+":
                    new_exons = sorted(exons, key=lambda x: x[0], reverse = False)
                    final = single_info + [new_exons]
                    if key in converted_dict:
                        converted_dict[key].append(final)
                    elif key not in converted_dict:
                        converted_dict.update({key:[final]})
                #if strand == -, the exons should decrease in position
                elif strand == "-":
                    final_exons = []
                    #this first sorts all exon pairs (start, end) in descending order
                    new_exons = sorted(exons, key=lambda x: x[0], reverse = True)
                    #because for - strand genes, the second position in a pair is actually the exon start, need to flip all exon pairs
                    for exon in new_exons:
                        exon.reverse()
                        final_exons.append(exon)
                    final = single_info + [final_exons]
                    if key in converted_dict:
                        converted_dict[key].append(final)
                    elif key not in converted_dict:
                        converted_dict.update({key:[final]})
    return converted_dict

#need to create exon trios for all transcripts with 3 or more exons
#create trio data base for all exons
#returns a dictionary with key = gene_id and value == [chr_num, strand, C1 start, C1 end, A start, A end, C2 start, C2 end]*number of exons-2
def create_trios():
    transcript_dict = convert_direction_combined_dict()
    trios_list = []
    trios_dict = {}
    for transcript in transcript_dict:
        single_transcript = transcript_dict[transcript]
        if len(single_transcript) == 1:
            single = single_transcript[0]
            if len(single) == 5:
                transcript_id = single[1]
                chr_num = single[2]
                strand = single[3]
            elif len(single) == 6:
                transcript_id = single[2]
                chr_num = single[3]
                strand = single[4]
            gene_id = single[0]
            exons = single[len(single)-1]
            x = 0
            #iterates through exons to create trios
            while x < len(exons)-2:
                exon_trio = [exons[x], exons[x+1], exons[x+2]]
                dict_value = [transcript_id, chr_num, strand, exon_trio]
                if gene_id in trios_dict:
                    trios_dict[gene_id].append(dict_value)
                elif gene_id not in trios_dict:
                    trios_dict.update({gene_id:[dict_value]})
                x += 1
        elif (len(single_transcript) == 5 or len(single_transcript) == 6) and isinstance(single_transcript[0], list) == False:
            if len(single_transcript) == 5:
                transcript_id = single_transcript[1]
                chr_num = single_transcript[2]
                strand = single_transcript[3]
            elif len(single_transcript) == 6:
                transcript_id = single_transcript[2]
                chr_num = single_transcript[3]
                strand = single_transcript[4]
            gene_id = single_transcript[0]
            exons = single_transcript[len(single_transcript)-1]
            x = 0
            #iterates through exons to create trios
            while x < len(exons)-2:
                exon_trio = [exons[x], exons[x+1], exons[x+2]]
                dict_value = [transcript_id, chr_num, strand, exon_trio]
                if gene_id in trios_dict:
                    trios_dict[gene_id].append(dict_value)
                elif gene_id not in trios_dict:
                    trios_dict.update({gene_id:[dict_value]})
                x += 1
        else:
            for single in single_transcript:
                if len(single) == 5:
                    transcript_id = single[1]
                    chr_num = single[2]
                    strand = single[3]
                elif len(single) == 6:
                    transcript_id = single[2]
                    chr_num = single[3]
                    strand = single[4]
                gene_id = single[0]
                exons = single[len(single)-1]
                x = 0
                #iterates through exons to create trios
                while x < len(exons)-2:
                    exon_trio = [exons[x], exons[x+1], exons[x+2]]
                    dict_value = [transcript_id, chr_num, strand, exon_trio]
                    if gene_id in trios_dict:
                        trios_dict[gene_id].append(dict_value)
                    elif gene_id not in trios_dict:
                        trios_dict.update({gene_id:[dict_value]})
                    x += 1
    return trios_dict


#need to filter trios dict to remove duplicates of exons (i.e. where two or more transcripts share same exon trio)
#will also need to address the differences in start and end positions for the first and last exon for transcripts from the same gene
#need to filter trios dict to get all internal exon trios and remove duplicates or condense pairs where first or last exon is different

#this function first filters out complete duplicates of exon trios from combined reduced dictionary
#only removes exact duplicates (where all three exon start and ends for a trio match another exon trio)
def filter_trios_duplicates():
    trios_dict = create_trios()
    filtered_trios_dict = {}
    exon_trio_list = []
    filtered_dict = {}
    for gene in trios_dict:
        single_gene = trios_dict[gene]
        #reads through exon trios for each gene (all combined from every isoform)
        for single in single_gene:
            chr_num = single[1]
            strand = single[2]
            exon_trio = single[3]
            exon_trio_tuple = []
            for exon_pair in exon_trio:
                exon_tuple = tuple(exon_pair)
                exon_trio_tuple.append(exon_tuple)
            exon_trio_list.append(exon_trio_tuple)
        exon_trio_tuple_list = [tuple(t) for t in exon_trio_list]
        #removes duplicate exon trios
        set_tuple = list(set(exon_trio_tuple_list))
        final_dictionary_value = [chr_num, strand, set_tuple]
        filtered_dict.update({gene:final_dictionary_value})
        exon_trio_list = []
    return filtered_dict


#now addressing start and end differences for first and last exon in a transcript
#first need a separate dictionary with the first exon (start, end) and the last exon (start, end) for every transcript
#dictionaries need to contain trios not just single exon
def pull_first_last_exon():
    combined_dict = convert_direction_combined_dict()
    first_exon_dict = {}
    last_exon_dict = {}
    for gene in combined_dict:
        single_gene = combined_dict[gene]
        #there are no instances of single isoform/trio that is bracketed by two sets of brackets = making length == 1, so don't have to ask for single_gene of length 1
        #this catches all of this genes that only have one isoform
        #ensembl genes (not matched to isoseq) will have a length of 5 and isoseq isoforms will be of length 6 (this is because the isoform ID is added)
        if (len(single_gene) == 5 or len(single_gene) == 6) and isinstance(single_gene[0], list) == False:
            if len(single_gene) == 5:
                strand = single_gene[3]
            elif len(single_gene) == 6:
                strand = single_gene[4]
            exons = single_gene[len(single_gene)-1]
            exon_trio_0 = exons[0:3]
            exon_trio_last = exons[len(exons)-3:len(exons)]
            exon_pair_0 = exons[0]
            exon_pair_last = exons[len(exons)-1]
            #exon positions should increase
            if strand == "+":
                #if exon pair at position 0 is the first exon, it should be less than the last exon start position
                if int(exon_pair_0[0]) < int(exon_pair_last[0]):
                    first_exon_trio = exon_trio_0
                    last_exon_trio = exon_trio_last
                #if exon pair at position 0 is greater than the exon pair at the last position, then the exon pair at position 0 is the last exon
                elif int(exon_pair_0[0]) > int(exon_pair_last[0]):
                    first_exon_trio = exon_trio_last
                    last_exon_trio = exon_trio_0
                if gene in first_exon_dict:
                    first_exon_dict[gene].append(first_exon_trio)
                elif gene not in first_exon_dict:
                    first_exon_dict.update({gene:[first_exon_trio]})
                if gene in last_exon_dict:
                    last_exon_dict[gene].append(last_exon_trio)
                elif gene not in last_exon_dict:
                    last_exon_dict.update({gene:[last_exon_trio]})
            #exon positions should decrease
            elif strand == "-":
                #if exon pair at position 0 is the first exon, it should be greater than the last exon start position
                if int(exon_pair_0[0]) > int(exon_pair_last[0]):
                    first_exon_trio = exon_trio_0
                    last_exon_trio = exon_trio_last
                #if exon pair at position 0 is less than the exon pair at the last position, then the exon pair at position 0 is the last exon
                elif int(exon_pair_0[0]) < int(exon_pair_last[0]):
                    first_exon_trio = exon_trio_last
                    last_exon_trio = exon_trio_0
                if gene in first_exon_dict:
                    first_exon_dict[gene].append(first_exon_trio)
                elif gene not in first_exon_dict:
                    first_exon_dict.update({gene:[first_exon_trio]})
                if gene in last_exon_dict:
                    last_exon_dict[gene].append(last_exon_trio)
                elif gene not in last_exon_dict:
                    last_exon_dict.update({gene:[last_exon_trio]})
        else:
            for s in single_gene:
                if len(s) == 5:
                    strand = s[3]
                elif len(s) == 6:
                    strand = s[4]
                exons = s[len(s)-1]
                exon_trio_0 = exons[0:3]
                exon_trio_last = exons[len(exons)-3:len(exons)]
                exon_pair_0 = exons[0]
                exon_pair_last = exons[len(exons)-1]
                #exon positions should increase
                if strand == "+":
                    #if exon pair at position 0 is the first exon, it should be less than the last exon start position
                    if int(exon_pair_0[0]) < int(exon_pair_last[0]):
                        first_exon_trio = exon_trio_0
                        last_exon_trio = exon_trio_last
                    #if exon pair at position 0 is greater than the exon pair at the last position, then the exon pair at position 0 is the last exon
                    elif int(exon_pair_0[0]) > int(exon_pair_last[0]):
                        first_exon_trio = exon_trio_last
                        last_exon_trio = exon_trio_0
                    if gene in first_exon_dict:
                        first_exon_dict[gene].append(first_exon_trio)
                    elif gene not in first_exon_dict:
                        first_exon_dict.update({gene:[first_exon_trio]})
                    if gene in last_exon_dict:
                        last_exon_dict[gene].append(last_exon_trio)
                    elif gene not in last_exon_dict:
                        last_exon_dict.update({gene:[last_exon_trio]})
                #exon positions should decrease
                elif strand == "-":
                    #if exon pair at position 0 is the first exon, it should be greater than the last exon start position
                    if int(exon_pair_0[0]) > int(exon_pair_last[0]):
                        first_exon_trio = exon_trio_0
                        last_exon_trio = exon_trio_last
                    #if exon pair at position 0 is less than the exon pair at the last position, then the exon pair at position 0 is the last exon
                    elif int(exon_pair_0[0]) < int(exon_pair_last[0]):
                        first_exon_trio = exon_trio_last
                        last_exon_trio = exon_trio_0
                    if gene in first_exon_dict:
                        first_exon_dict[gene].append(first_exon_trio)
                    elif gene not in first_exon_dict:
                        first_exon_dict.update({gene:[first_exon_trio]})
                    if gene in last_exon_dict:
                        last_exon_dict[gene].append(last_exon_trio)
                    elif gene not in last_exon_dict:
                        last_exon_dict.update({gene:[last_exon_trio]})
    return first_exon_dict, last_exon_dict


#need to filter first and last exons before filtering trios
#this filtering is collapsing varying start sites and end sites to retain the start/stop position that is the farthest upstream/downstream
def filter_first_exon():
    first_exons, last_exons = pull_first_last_exon()
    filtered_first_exons = {}
    for gene in first_exons:
        single_gene = first_exons[gene]
        #if there is only 1 isoform for a gene, there will be only 1 first exon
        #then the exon start/end can immediately be added to the final dictionary
        if len(single_gene) == 1:
            single = single_gene[0]
            filtered_first_exons.update({gene:single})
        #if there is more than 1 isoform for a gene, there will be many exons that count as the start exon
        #need to collapse down any starting exons that end at the same position (alternative start sites), but only if the rest of the trio matches
        elif len(single_gene) > 1:
            #need to pull out all exon positions except first to figure out if these are the same
            non_first_exon_start = []
            for single in single_gene:
                non_start_pos = single[1:len(single)]
                first_exon_end_pos = single[0][1]
                non_start_pos.insert(0, [first_exon_end_pos])
                combined = [tuple(t) for t in non_start_pos]
                non_first_exon_start.append(combined)
            tuple_non_first_exon_start = [tuple(u) for u in non_first_exon_start]
            set_tuple_non_first_exon_start = set(tuple_non_first_exon_start)
            #if the length of the single gene value and the set of the tuple non first exon start list, this means that the exon trios are different for all of the first trios, so all should be kept
            if len(single_gene) == len(set_tuple_non_first_exon_start):
                for exon_trio in single_gene:
                    if gene in filtered_first_exons:
                        filtered_first_exons[gene].append(exon_trio)
                    elif gene not in filtered_first_exons:
                        filtered_first_exons.update({gene:[exon_trio]})
            #if the length of the set_tuple_non_first_exon_start == 0, this means all of the first exons cover the same trios except for the starting position of the first exon.
            elif len(set_tuple_non_first_exon_start) == 1:
                exon_starts = []
                exon_ends = []
                #iterates through the first exons from all the isoforms
                for s in single_gene:
                    exon_starts.append(s[0][0])
                    exon_ends.append(s[0][1])
                #want to pull the start position with that is the farthest upstream of gene for this
                #examines a single exon pair to determine direction
                #if gene is on + strand, single[0][0] < single[0][1]
                if int(s[0][0]) < int(s[0][1]):
                    final_start = min(exon_starts)
                #if gene is on - strand, single[0][0] > single[0][1]
                elif int(single[0][0]) > int(single[0][1]):
                    final_start = max(exon_starts)
                final_first_exon_full = [final_start, exon_ends[0]]
                final_exon_trio = [final_first_exon_full, s[1], s[2]]
                filtered_first_exons.update({gene:final_exon_trio})
            #if the length of the set_tuple doesn't equal 1 or the length of the single gene, this means some of the trios should be condensed and other should be kept
            else:
                #pulls all indeces that represent a duplicate value from the list of values with [[first exon end pos], [second exon start, second exon end pos], [third exon start, third exon end pos]]
                duplicates_index = [i for i, x in enumerate(non_first_exon_start) if non_first_exon_start.count(x) > 1]
                exon_starts = []
                exon_ends = []
                #iterates through the first exons from all the isoforms
                for v in single_gene:
                    exon_starts.append(v[0][0])
                    exon_ends.append(v[0][1])
                possible_exon_starts = []
                for index, value in enumerate(exon_starts):
                    #if starting value is in duplicates index, add the ending value to a new list
                    if index in duplicates_index:
                        possible_exon_starts.append(exon_starts[index])
                    #if the starting value is not in the duplicates list, can just add start, end combo to filtered exon dictionary
                    else:
                        final_first_exon_full = [value, exon_ends[index]]
                        final_exon_trio = [final_first_exon_full, single_gene[index][1], single_gene[index][2]]
                        if gene in filtered_first_exons:
                            filtered_first_exons[gene].append(final_exon_trio)
                        elif gene not in filtered_first_exons:
                            filtered_first_exons.update({gene:[final_exon_trio]})
                #provides exon end values that were duplicated
                duplicates_value = list(set([exon_ends[a] for a in duplicates_index]))
                #loops through duplicates values
                for val in duplicates_value:
                    #provides the index of duplicate values in the exon ends; this is mainly for when there is more than one value that was duplicated
                    single_dup_index_list = [b for b in range(len(exon_ends)) if exon_ends[b] == val]
                    single_dup_starts_list = []
                    single_dup_ends_list = []
                    for ind in single_dup_index_list:
                        single_dup_starts_list.append(exon_starts[ind])
                        single_dup_ends_list.append(exon_ends[ind])
                    #if start is < than end, this means the gene is on the + strand
                    if single_dup_starts_list[0] < single_dup_ends_list[0]:
                        final_start = min(single_dup_starts_list)
                        final_first_exon_full = [final_start, single_dup_ends_list[0]]
                    #if start is > than end, this means the gene is on the - strand
                    elif single_dup_starts_list[0] > single_dup_ends_list[0]:
                        final_start = max(single_dup_starts_list)
                        final_first_exon_full = [final_start, single_dup_ends_list[0]]
                    final_exon_trio = [final_first_exon_full, single_gene[ind][1], single_gene[ind][2]]
                    if gene in filtered_first_exons:
                        filtered_first_exons[gene].append(final_exon_trio)
                    elif gene not in filtered_first_exons:
                        filtered_first_exons.update({gene:[final_exon_trio]})
    return filtered_first_exons


def filter_last_exon():
    first_exons, last_exons = pull_first_last_exon()
    filtered_last_exons = {}
    for gene in last_exons:
        single_gene = last_exons[gene]
        print(single_gene)
        #if there is only 1 isoform for a gene, there will be only 1 last exon
        #then the exon start/end can immediately be added to the final dictionary
        '''if len(single_gene) == 1:
            single = single_gene[0]
            filtered_last_exons.update({gene:single})
        #if there is more than 1 isoform for a gene, there will be many exons that count as the last exon
        #need to collapse down any ending exons that start at the same position (alternative end sites), but only if the rest of the trio matches
        elif len(single_gene) > 1:
            #need to pull out all exon positions except last to figure out if these are the same
            non_last_exon_end = []
            for single in single_gene:
                print(single)
                non_last_pos = single[0:len(single)-1]
                last_exon_end_pos = single[0][1]
                non_start_pos.insert(0, [first_exon_end_pos])
                combined = [tuple(t) for t in non_start_pos]
                non_first_exon_start.append(combined)
            tuple_non_first_exon_start = [tuple(u) for u in non_first_exon_start]
            set_tuple_non_first_exon_start = set(tuple_non_first_exon_start)
            #if the length of the single gene value and the set of the tuple non first exon start list, this means that the exon trios are different for all of the first trios, so all should be kept
            if len(single_gene) == len(set_tuple_non_first_exon_start):
                for exon_trio in single_gene:
                    if gene in filtered_first_exons:
                        filtered_first_exons[gene].append(exon_trio)
                    elif gene not in filtered_first_exons:
                        filtered_first_exons.update({gene:[exon_trio]})
            #if the length of the set_tuple_non_first_exon_start == 0, this means all of the first exons cover the same trios except for the starting position of the first exon.
            elif len(set_tuple_non_first_exon_start) == 1:
                exon_starts = []
                exon_ends = []
                #iterates through the first exons from all the isoforms
                for s in single_gene:
                    exon_starts.append(s[0][0])
                    exon_ends.append(s[0][1])
                #want to pull the start position with that is the farthest upstream of gene for this
                #examines a single exon pair to determine direction
                #if gene is on + strand, single[0][0] < single[0][1]
                if int(s[0][0]) < int(s[0][1]):
                    final_start = min(exon_starts)
                #if gene is on - strand, single[0][0] > single[0][1]
                elif int(single[0][0]) > int(single[0][1]):
                    final_start = max(exon_starts)
                final_first_exon_full = [final_start, exon_ends[0]]
                final_exon_trio = [final_first_exon_full, s[1], s[2]]
                filtered_first_exons.update({gene:final_exon_trio})
            #if the length of the set_tuple doesn't equal 1 or the length of the single gene, this means some of the trios should be condensed and other should be kept
            else:
                #pulls all indeces that represent a duplicate value from the list of values with [[first exon end pos], [second exon start, second exon end pos], [third exon start, third exon end pos]]
                duplicates_index = [i for i, x in enumerate(non_first_exon_start) if non_first_exon_start.count(x) > 1]
                exon_starts = []
                exon_ends = []
                #iterates through the first exons from all the isoforms
                for v in single_gene:
                    exon_starts.append(v[0][0])
                    exon_ends.append(v[0][1])
                possible_exon_starts = []
                for index, value in enumerate(exon_starts):
                    #if starting value is in duplicates index, add the ending value to a new list
                    if index in duplicates_index:
                        possible_exon_starts.append(exon_starts[index])
                    #if the starting value is not in the duplicates list, can just add start, end combo to filtered exon dictionary
                    else:
                        final_first_exon_full = [value, exon_ends[index]]
                        final_exon_trio = [final_first_exon_full, single_gene[index][1], single_gene[index][2]]
                        if gene in filtered_first_exons:
                            filtered_first_exons[gene].append(final_exon_trio)
                        elif gene not in filtered_first_exons:
                            filtered_first_exons.update({gene:[final_exon_trio]})
                #provides exon end values that were duplicated
                duplicates_value = list(set([exon_ends[a] for a in duplicates_index]))
                #loops through duplicates values
                for val in duplicates_value:
                    #provides the index of duplicate values in the exon ends; this is mainly for when there is more than one value that was duplicated
                    single_dup_index_list = [b for b in range(len(exon_ends)) if exon_ends[b] == val]
                    single_dup_starts_list = []
                    single_dup_ends_list = []
                    for ind in single_dup_index_list:
                        single_dup_starts_list.append(exon_starts[ind])
                        single_dup_ends_list.append(exon_ends[ind])
                    #if start is < than end, this means the gene is on the + strand
                    if single_dup_starts_list[0] < single_dup_ends_list[0]:
                        final_start = min(single_dup_starts_list)
                        final_first_exon_full = [final_start, single_dup_ends_list[0]]
                    #if start is > than end, this means the gene is on the - strand
                    elif single_dup_starts_list[0] > single_dup_ends_list[0]:
                        final_start = max(single_dup_starts_list)
                        final_first_exon_full = [final_start, single_dup_ends_list[0]]
                    final_exon_trio = [final_first_exon_full, single_gene[ind][1], single_gene[ind][2]]
                    if gene in filtered_first_exons:
                        filtered_first_exons[gene].append(final_exon_trio)
                    elif gene not in filtered_first_exons:
                        filtered_first_exons.update({gene:[final_exon_trio]})
    return filtered_first_exons'''

filter_last_exon()

#need to only filter out first and last exon where the rest of the trio is the same
def filter_trios_first_exon():
    filtered_dict = filter_trios_duplicates()
    filtered_first_exons = filter_first_exon()
    final_first_exons = {}
    final_remaining_exons_dict = {}
    for gene in filtered_dict:
        single_gene = filtered_dict[gene]
        single_first_exons = filtered_first_exons[gene]
        chr_num = single_gene[0]
        strand = single_gene[1]
        single_gene_exon_trios = single_gene[2]
        first_exon_test_list = []
        remaining_exon_trio_list = []
        for single_trio in single_gene_exon_trios:
            single_trio_list = [list(ele) for ele in single_trio]
            first_exon = single_trio_list[0]
            first_exon_start = first_exon[0]
            first_exon_end = first_exon[1]
            if len(single_first_exons) == 2 and isinstance(single_first_exons[0], list) == False:
                if first_exon_end == single_first_exons[1]:
                    first_exon_test_list.append(single_trio_list)
            elif isinstance(single_first_exons[0], list) == True:
                for value in single_first_exons:
                    if first_exon_end == value[1]:
                        first_exon_test_list.append(single_trio_list)
        #reduce first exon
        if isinstance(single_first_exons[0], list) == False:
            for exon_trio in first_exon_test_list:
                first_exon_pair = exon_trio[0]
                #if first_exon_pair exactly matches the first exon from single_first_exons; add to dictionary
                #since the single_first_exons contains the most upstream alternative start site, do not have to sort through the rest of the values
                if first_exon_pair == single_first_exons:
                    final_exon_trio = exon_trio
                    final_first_exons.update({gene:final_exon_trio})
        elif isinstance(single_first_exons[0], list) == True:
            for exon_pair in single_first_exons:
                for exon_trio in first_exon_test_list:
                    first_exon_pair = exon_trio[0]
                    #if first_exon_pair exactly matches the first exon from single_first_exons; add to dictionary
                    if first_exon_pair == exon_pair:
                        final_exon_trio = exon_trio
                        if gene in final_first_exons:
                            final_first_exons[gene].append(final_exon_trio)
                        elif gene not in final_first_exons:
                            final_first_exons.update({gene:[final_exon_trio]})
    return final_first_exons


def filter_trios_last_exon():
    filtered_dict = filter_trios_duplicates()
    filtered_first_exons = filter_trios_first_exon()
    filtered_last_exons = filter_last_exon()
    for gene in filtered_dict:
        single_gene = filtered_dict[gene]
        final_first_exons = filtered_first_exons[gene]
        single_last_exons = filtered_last_exons[gene]
        chr_num = single_gene[0]
        strand = single_gene[1]
        single_gene_exon_trios = single_gene[2]
        last_exon_test_list = []
        remaining_exon_trio_list = []
        for single_trio in single_gene_exon_trios:
            single_trio_list = [list(ele) for ele in single_trio]
            last_exon = single_trio_list[2]
            last_exon_start = last_exon[0]
            last_exon_end = last_exon[1]
            if len(single_last_exons) == 2 and isinstance(single_last_exons[0], list) == False:
                if last_exon_start == single_last_exons[0]:
                    last_exon_test_list.append(single_trio_list)
                else:
                    remaining_exon_trio_list.append(single_trio_list)
            elif isinstance(single_last_exons[0], list) == True:
                for value in single_last_exons:
                    if last_exon_start == value[0]:
                        last_exon_test_list.append(single_trio_list)
                    else:
                        remaining_exon_trio_list.append(single_trio_list)
        #removing exon trios from remaining exon trio list if in final_first_exon list and last_exon_test_list
        if len(final_first_exons[0]) == 2:
            if final_first_exons in remaining_exon_trio_list:
                remaining_exon_trio_list.remove(final_first_exons)
        elif len(final_first_exons[0]) != 2:
            for trio in final_first_exons:
                if trio in remaining_exon_trio_list:
                    remaining_exon_trio_list.remove(trio)
        for trio_2 in last_exon_test_list:
            if trio_2 in remaining_exon_trio_list:
                remaining_exon_trio_list.remove(trio_2)
        tuple_remaining_exon_trios = []
        for val in remaining_exon_trio_list:
            tuple_val = [tuple(t) for t in val]
            tuple_remaining_exon_trios.append(tuple_val)
        final_tuple_remaining_exons = [tuple(x) for x in tuple_remaining_exon_trios]
        set_remaining_exon_trios = list(set(final_tuple_remaining_exons))
        #convert remaining exon trios tuples to list of lists
        final_remaining_exons_list = []
        for tup in set_remaining_exon_trios:
            single_list = [list(g) for g in tup]
            final_remaining_exons_list.append(single_list)
        #reduce last exon
        final_last_exons = []
        if isinstance(single_last_exons[0], list) == False:
            for exon_trio in last_exon_test_list:
                last_exon_pair = exon_trio[2]
                #if last_exon_pair exactly matches the last exon from single_last_exons; add to dictionary
                #since the single_first_exons contains the most upstream alternative start site, do not have to sort through the rest of the values
                if last_exon_pair == single_last_exons:
                    final_last_exons= exon_trio
        elif isinstance(single_last_exons[0], list) == True:
            for exon_pair in single_last_exons:
                for exon_trio in last_exon_test_list:
                    last_exon_pair = exon_trio[2]
                    #if last_exon_pair exactly matches the last exon from single_first_exons; add to dictionary
                    if last_exon_pair == exon_pair:
                        final_exon_trio = exon_trio
                        final_last_exons.append(final_exon_trio)
        print(gene)
        print(chr_num)
        print(strand)
        print(len(final_first_exons))
        print(len(final_last_exons))
        print(len(final_remaining_exons_list))


#filter_trios_last_exon()

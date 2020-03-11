##removing duplicates identified from BLAST results in combined sexes isoseq data when run as blastp
#if isoforms have the same number of amino acids between isoforms, these will be considered duplicates
#to run script: python3 Remove_duplicates_protein_BLAST.py <protein blast results for combined sexes> <combined sexes faa file> <collapsed isoforms output file> <removed isoforms output file> 
#Author: Alice Naftaly, March 2020

import sys

#read in blast output
#returns dictionary with key = isoform id and value = alignment
def read_blast():
    blast_output = sys.argv[1]
    blast_dict = {}
    with open(blast_output, 'r') as blast_results:
        for line in blast_results:
            if line.startswith("#"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                if isoform_id in blast_dict:
                    blast_dict[isoform_id].append(new_line)
                elif isoform_id not in blast_dict:
                    blast_dict.update({isoform_id:[new_line]})
    print("Read BLAST Output File")
    return blast_dict

#read in faa File
#want isoform id and number of amino acids in isoform
def read_faa():
    faa_file = sys.argv[2]
    faa_dict = {}
    with open(faa_file, 'r') as amino_acid_seqs:
        for line in amino_acid_seqs:
            if line.startswith(">"):
                new_line = line.split("\t")
                isoform_id = new_line[0].strip(">")
                line_info = new_line[1].split("|")
                num_amino_acids = line_info[2].strip("_aa")
                cds_start = line_info[4]
                cds_end = line_info[5].strip("\n")
                final = [num_amino_acids, cds_start, cds_end]
                faa_dict.update({isoform_id:final})
    print("Read faa file")
    return faa_dict


#identifying isoforms that code for the same amino acid
def find_duplicates():
    blast_out = read_blast()
    amino_acid_dict = read_faa()
    all_isoforms = []
    matches_dict = {}
    x = 0
    for isoform in blast_out:
        single_blast = blast_out[isoform]
        single_isoform_aa = amino_acid_dict[isoform]
        single_isoform_aa_length = int(single_isoform_aa[0])
        isoform_number = isoform.strip("PB.")
        split_isoform = isoform_number.split(".")
        final_isoform = int(split_isoform[0])
        matches = []
        for match in single_blast:
            alignment_match_isoform = match[2]
            alignment_length = int(match[4])
            if alignment_length == single_isoform_aa_length:
                match_isoform_aa_length = int(amino_acid_dict[alignment_match_isoform][0])
                if match_isoform_aa_length == single_isoform_aa_length:
                    matches.append(alignment_match_isoform)
        matches.insert(0, isoform)
        set_matches = list(set(matches))
        if isoform in all_isoforms:
            continue
        elif isoform not in all_isoforms:
            if len(set_matches) == 1:
                all_isoforms.append(set_matches[0])
                if final_isoform in matches_dict:
                    matches_dict[final_isoform].append(set_matches[0])
                elif final_isoform not in matches_dict:
                    matches_dict.update({final_isoform:[set_matches[0]]})
            elif len(set_matches) > 1:
                for iso_match in set_matches:
                    all_isoforms.append(iso_match)
                if final_isoform in matches_dict:
                    matches_dict[final_isoform].append(set_matches)
                elif final_isoform not in matches_dict:
                    matches_dict.update({final_isoform:[set_matches]})
    final_dictionary = {}
    for key in matches_dict:
        single_matches = matches_dict[key]
        if len(single_matches) == 1 and isinstance(single_matches[0], list) == False:
            if key in final_dictionary:
                final_dictionary[key].append(single_matches[0])
            elif key not in final_dictionary:
                final_dictionary.update({key:[single_matches[0]]})
        elif len(single_matches) == 1 and isinstance(single_matches[0], list) == True:
            if key in final_dictionary:
                final_dictionary[key].append(single_matches)
            elif key not in final_dictionary:
                final_dictionary.update({key:single_matches})
        else:
            for value in single_matches:
                if isinstance(value, str) == True:
                    new_value = [value]
                    if key in final_dictionary:
                        final_dictionary[key].append(new_value)
                    elif key not in final_dictionary:
                        final_dictionary.update({key:[new_value]})
                else:
                    if key in final_dictionary:
                        final_dictionary[key].append(value)
                    elif key not in final_dictionary:
                        final_dictionary.update({key:[value]})
    print("Created matches dictionary")
    return final_dictionary


#want to have the same genes be written together
#also will add isoforms that match more than 1 gene to a separate list to be removed
def order_dictionary():
    matches = find_duplicates()
    genes_to_remove = []
    genes_to_keep = []
    gene_numbers = sorted(matches.keys())
    for gene in gene_numbers:
        single_gene = matches[gene]
        if len(single_gene) == 1 and isinstance(single_gene[0], list) == False:
            genes_to_keep.append(gene)
        elif len(single_gene) == 1 and isinstance(single_gene[0], list) == True:
            single = single_gene[0]
            counts = 0
            for value in single:
                split_value = value.split(".")
                ind_gene = int(split_value[1])
                if ind_gene == gene:
                    counts += 1
            if counts == len(single):
                genes_to_keep.append(gene)
            elif counts != len(single):
                genes_to_remove.append(single)
        else:
            flattened_list = [item for sublist in single_gene for item in sublist]
            counts = 0
            for val in flattened_list:
                split_val = val.split(".")
                ind_gene = int(split_val[1])
                if ind_gene == gene:
                    counts += 1
            if counts == len(flattened_list):
                genes_to_keep.append(gene)
            elif counts != len(flattened_list):
                genes_to_remove.append(flattened_list)
    #flatten genes to remove list
    flattened_list_genes_to_remove = [a for sub_item in genes_to_remove for a in sub_item]
    final_genes_to_remove = list(set(flattened_list_genes_to_remove))
    return genes_to_keep, final_genes_to_remove


#write genes to keep to a file in order
def write_genes_to_keep():
    output = sys.argv[3]
    genes_to_keep, genes_to_remove = order_dictionary()
    matches = find_duplicates()
    with open(output, 'a') as out:
        for gene in genes_to_keep:
            if gene in matches:
                single_gene = matches[gene]
                if len(single_gene) == 1 and isinstance(single_gene[0], list) == False:
                    if single_gene[0] in genes_to_remove:
                        continue
                    elif single_gene[0] not in genes_to_remove:
                        out.write(single_gene[0] + "\n")
                elif len(single_gene) == 1 and isinstance(single_gene[0], list) == True:
                    keep_genes = []
                    for value in single_gene[0]:
                        if value in genes_to_remove:
                            continue
                        elif value not in genes_to_remove:
                            keep_genes.append(value)
                    final = "\t".join(keep_genes)
                    out.write(final + "\n")
                elif len(single_gene) > 1:
                    for val in single_gene:
                        if len(val) == 1:
                            if val[0] in genes_to_remove:
                                continue
                            elif val[0] not in genes_to_remove:
                                out.write(val[0] + "\n")
                        elif len(val) > 1:
                            kept_genes = []
                            for v in val:
                                if v in genes_to_remove:
                                    continue
                                elif v not in genes_to_remove:
                                    kept_genes.append(v)
                            final = "\t".join(kept_genes)
                            out.write(final + "\n")

#write isoforms to remove from the comparison
def write_genes_removed():
    output = sys.argv[4]
    genes_to_keep, genes_to_remove = order_dictionary()
    with open(output, 'a') as out:
        for value in genes_to_remove:
            out.write(value + "\n")

#calls both functions to write
def call():
    kept_genes = write_genes_to_keep()
    removed_genes = write_genes_removed()

call()

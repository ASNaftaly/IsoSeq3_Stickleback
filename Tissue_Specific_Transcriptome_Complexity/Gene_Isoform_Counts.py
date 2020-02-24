#Calculating various stats for single tissue analyses after filtering.
#Want to count the total number of genes (will need to convert gene id from single tissue gene to combined sexes gene), how many isoforms, the number of isoforms per gene, how many annotated vs novel genes
#to run script: python3 Gene_Isoform_Counts.py <read single tissue analysis classification file> <combined sexes classification file> <output, isoform counts with Gene.ID\tIsoform.Counts> <output, exon counts with Isoform.ID\tExon.Counts>
#Author: Alice Naftaly, February 2020

import sys

#read single tissue classification file
#returns dictionary with key = new isoform id and value == whole line
def read_single_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                new_isoform = new_line[0]
                class_dict.update({new_isoform:new_line})
    return class_dict

#read combined sexes classification file
#return dictionary with key = isoform id and value = gene id
def read_combined_class_genes():
    class_file = sys.argv[2]
    combined_class_dict = {}
    with open(class_file, 'r') as combined_class:
        for line in combined_class:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[6]
                combined_class_dict.update({isoform:gene_id})
    return combined_class_dict

#return dictionary with key = isoform id and value = transcript id
def read_combined_class_transcripts():
    class_file = sys.argv[2]
    combined_class_dict = {}
    with open(class_file, 'r') as combined_class:
        for line in combined_class:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                transcript_id = new_line[7]
                combined_class_dict.update({isoform:transcript_id})
    return combined_class_dict

#Get gene count
#sorts through genes and makes sure the ids match the combined sexes analysis gene id
#Then provides count of total genes (novel and annotated) and prints out summary statements
#returns dictionary with key == gene id and value == whole line
def gene_count():
    class_dict = read_single_class()
    combined_class = read_combined_class_genes()
    gene_dict = {}
    excluded_isoforms = 0
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        single_combined_gene = combined_class[isoform]
        gene_id = single_isoform[7]
        if single_combined_gene == gene_id:
            if gene_id in gene_dict:
                gene_dict[gene_id].append(single_isoform)
            elif gene_id not in gene_dict:
                gene_dict.update({gene_id:[single_isoform]})
        elif gene_id.startswith("novel") and single_combined_gene.startswith("novel"):
            if single_combined_gene in gene_dict:
                gene_dict[single_combined_gene].append(single_isoform)
            elif single_combined_gene not in gene_dict:
                gene_dict.update({single_combined_gene:[single_isoform]})
        else:
            excluded_isoforms += 1
    return gene_dict

#only writes summary statement
def gene_count_summary():
    class_dict = read_single_class()
    combined_class = read_combined_class_genes()
    gene_dict = {}
    excluded_isoforms = 0
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        single_combined_gene = combined_class[isoform]
        gene_id = single_isoform[7]
        if single_combined_gene == gene_id:
            if gene_id in gene_dict:
                gene_dict[gene_id].append(single_isoform)
            elif gene_id not in gene_dict:
                gene_dict.update({gene_id:[single_isoform]})
        elif gene_id.startswith("novel") and single_combined_gene.startswith("novel"):
            if single_combined_gene in gene_dict:
                gene_dict[single_combined_gene].append(single_isoform)
            elif single_combined_gene not in gene_dict:
                gene_dict.update({single_combined_gene:[single_isoform]})
        else:
            excluded_isoforms += 1
    print("Number of genes identified")
    print(len(gene_dict))
    print("Number of genes/isoforms excluded because of incorrect match")
    print(excluded_isoforms)
    annotated_genes = 0
    novel_genes = 0
    other_genes = 0
    for key in gene_dict:
        if key.startswith("ENSGACG"):
            annotated_genes += 1
        elif key.startswith("novel"):
            novel_genes += 1
        else:
            other_genes += 1
    print("Number of annotated genes")
    print(annotated_genes)
    print("Number of novel genes")
    print(novel_genes)
    print("Genes that don't fit either category above")
    print(other_genes)

#count number of isoforms per gene
#counts the number of isoforms per gene and returns summary statements
#also returns key = gene and value == number of isoforms per gene
#will write this dictionary to a new output file
def count_isoforms():
    gene_dict = gene_count()
    isoform_count_dict = {}
    isoform_dict = {}
    for gene in gene_dict:
        single_gene = gene_dict[gene]
        if len(single_gene) == 1:
            isoform_count_dict.update({gene:"1"})
        elif len(single_gene) > 1:
            num_isoforms = len(single_gene)
            isoform_count_dict.update({gene:str(num_isoforms)})
    one_isoform_total = 0
    two_isoforms_total = 0
    three_isoforms_total = 0
    four_isoforms_total = 0
    five_isoforms_total = 0
    six_or_more_isoforms_total = 0
    for key in isoform_count_dict:
        single_key = isoform_count_dict[key]
        if int(single_key) == 1:
            one_isoform_total += 1
        elif int(single_key) == 2:
            two_isoforms_total += 1
        elif int(single_key) == 3:
            three_isoforms_total += 1
        elif int(single_key) == 4:
            four_isoforms_total += 1
        elif int(single_key) == 5:
            five_isoforms_total += 1
        elif int(single_key) >= 6:
            six_or_more_isoforms_total += 1
    print("Number of genes with 1 isoform")
    print(one_isoform_total)
    print("Number of genes with 2 isoforms")
    print(two_isoforms_total)
    print("Number of genes with 3 isoforms")
    print(three_isoforms_total)
    print("Number of genes with 4 isoforms")
    print(four_isoforms_total)
    print("Number of genes with 5 isoforms")
    print(five_isoforms_total)
    print("Number of genes with 6 or more isoforms")
    print(six_or_more_isoforms_total)
    return isoform_count_dict


#number of exons per isoforms
#prints summary statements
#returns exon count dictionary with key = isoform ID and value = number of exons
def exon_count():
    genes = gene_count()
    exon_dict = {}
    for gene in genes:
        single_gene = genes[gene]
        if len(single_gene) == 1:
            single = single_gene[0]
            isoform = single[0]
            num_exons = single[5]
            exon_dict.update({isoform:num_exons})
        elif len(single_gene) > 1:
            for single in single_gene:
                isoform = single[0]
                num_exons = single[5]
                exon_dict.update({isoform:num_exons})
    one_exon_total = 0
    two_to_five_exons_total = 0
    six_to_ten_exons_total = 0
    eleven_to_fifteen_exons_total = 0
    sixteen_to_twenty_exons_total = 0
    more_than_twenty_exons_total = 0
    for key in exon_dict:
        single_key = exon_dict[key]
        if int(single_key) == 1:
            one_exon_total += 1
        elif 2 <= int(single_key) <= 5:
            two_to_five_exons_total += 1
        elif 6 <= int(single_key) <= 10:
            six_to_ten_exons_total += 1
        elif 11 <= int(single_key) <= 15:
            eleven_to_fifteen_exons_total += 1
        elif 16 <= int(single_key) <= 20:
            sixteen_to_twenty_exons_total += 1
        elif int(single_key) > 20:
            more_than_twenty_exons_total += 1
    print("Number of isoforms with 1 exon")
    print(one_exon_total)
    print("Number of isoforms with 2-5 exons")
    print(two_to_five_exons_total)
    print("Number of isoforms with 6-10 exons")
    print(six_to_ten_exons_total)
    print("Number of isoforms with 11-15 exons")
    print(eleven_to_fifteen_exons_total)
    print("Number of isoforms with 16-20 exons")
    print(sixteen_to_twenty_exons_total)
    print("Number of isoforms with more than 20 exons")
    print(more_than_twenty_exons_total)
    return exon_dict

#write isoform counts to output files
def write_isoform_counts():
    isoform_counts = count_isoforms()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Gene.ID\tIsoform.Count\n"
        out.write(header)
        for isoform in isoform_counts:
            single_isoform = isoform_counts[isoform]
            final = "%s\t%s\n" % (str(isoform),str(single_isoform))
            out.write(final)

#write exon counts to output files
def write_exon_counts():
    exon_counts = exon_count()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Isoform.ID\tExon.Count\n"
        out.write(header)
        for isoform in exon_counts:
            single_isoform = exon_counts[isoform]
            final = "%s\t%s\n" % (str(isoform),str(single_isoform))
            out.write(final)

#call all functions:
def call():
    gene_counts = gene_count_summary()
    isoform_counts = write_isoform_counts()
    exon_counts = write_exon_counts()

call()

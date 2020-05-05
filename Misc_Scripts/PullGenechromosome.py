#pulling gene ids and chromosomes from classification files
#will need to read in classification files for Isoseq (combined sexes) and gtf from ensembl
#will then get gene_id and chr number and sort through these so there is no redundancy and write this to an output file
#to run script: PullGenechromosome.py <ensembl gtf file> <isoseq classification file> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read ensembl gtf file
#only need chromosome number and gene id
def read_ensembl_gtf():
    ensembl_file = sys.argv[1]
    ensembl_dict = {}
    with open(ensembl_file, 'r') as ensembl_gtf:
        for line in ensembl_gtf:
            new_line = line.split("\t")
            type = new_line[2]
            if type == "gene":
                chr_num = new_line[0]
                gene_info = new_line[8].split(";")
                gene_id_full = gene_info[0].split(" ")
                gene_id = gene_id_full[1].strip("\"")
                ensembl_dict.update({gene_id:chr_num})
    return ensembl_dict



#read isoseq classification file
def read_isoseq_class():
    isoseq_file = sys.argv[2]
    isoseq_dict = {}
    with open(isoseq_file, 'r') as isoseq_class:
        for line in isoseq_class:
            new_line = line.split("\t")
            chr_num = new_line[1]
            gene_id = new_line[6]
            isoseq_dict.update({gene_id:chr_num})
    return isoseq_dict


#combined ensembl and isoseq dictionaries
#reduces ensembl and isoseq dictionaries so only 1 copy of a gene is present in the final dictionary
def combine():
    ensembl_dict = read_ensembl_gtf()
    isoseq_dict = read_isoseq_class()
    final_dict = {}
    for gene in ensembl_dict:
        if gene in isoseq_dict:
            single_gene = ensembl_dict[gene]
            final_dict.update({gene:single_gene})
        elif gene not in isoseq_dict:
            single_gene = ensembl_dict[gene]
            final_dict.update({gene:single_gene})
    for gene2 in isoseq_dict:
        if gene2 not in ensembl_dict:
            single_gene = isoseq_dict[gene2]
            final_dict.update({gene2:single_gene})
    return final_dict


#write output
def write():
    final_dict = combine()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Gene.ID\tChr.Num\n"
        out.write(header)
        for gene in final_dict:
            chr_num = final_dict[gene]
            final = "%s\t%s\n" % (str(gene), str(chr_num))
            out.write(final)

write()

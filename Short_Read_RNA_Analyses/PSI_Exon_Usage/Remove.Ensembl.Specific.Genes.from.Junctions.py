#pulling isoseq specific junctions by removing ensembl specific genes from all junctions file
#will read in ensembl specific gene ids and all junctions
#to run script: python3 Remove.Ensembl.Specific.Genes.from.Junctions.py <ensembl specific genes file> <all junctions file> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read in ensembl specific genes file
#returns list of ensembl gene ids
def read_ensembl_genes():
    ensembl_file = sys.argv[1]
    ensembl_gene_list = []
    with open(ensembl_file, 'r') as ensembl:
        for line in ensembl:
            new_line = line.strip("\n")
            ensembl_gene_list.append(new_line)
    return ensembl_gene_list


#read in all junctions file
#returns dictionary with key == gene id and value == line for each junction
def read_junctions():
    junctions_file = sys.argv[2]
    junctions_dict = {}
    with open(junctions_file, 'r') as junctions:
        for line in junctions:
            new_line = line.split()
            junction_id = new_line[3].split("_")
            gene_id = junction_id[0]
            if gene_id in junctions_dict:
                junctions_dict[gene_id].append(line)
            elif gene_id not in junctions_dict:
                junctions_dict.update({gene_id:[line]})
    return junctions_dict


#remove ensembl specific genes
#returns dictionary with key == gene and value == lien for each junction
def remove_ensembl_genes():
    ensembl_genes = read_ensembl_genes()
    junctions_dict = read_junctions()
    for gene in ensembl_genes:
        if gene in junctions_dict:
            del junctions_dict[gene]
    return junctions_dict


#write new junction dict (with just isoseq specific genes)
def write_new_junctions_file():
    junctions_dict = remove_ensembl_genes()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for junction in junctions_dict:
            single_junction = junctions_dict[junction]
            for value in single_junction:
                out.write(value)

write_new_junctions_file()

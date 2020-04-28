#combination of Pull_GO_terms_isoseq.py and Pull_GO_terms_novel_genes.py
#will read in classification file with all isoforms for a specific analysis, ensembl genes with GO terms, and Novel isoforms with GO terms from interproscan
#to run script: python3 Pull_GO_terms_all.py <classification file> <ensembl genes go terms file> <novel genes go terms file> <output file>
#Author: Alice Naftaly, April 2020

import sys


#read classification file
#need gene id and transcript id
#returns dictionary with key == "Ensembl_Annotated" or "Novel" and value == list of gene ids
def read_class_file():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                gene_id = new_line[6]
                transcript_id = new_line[7]
                if gene_id.startswith("ENSGACG") and transcript_id.startswith("ENSGACT"):
                    if "Ensembl_Annotated" in class_dict:
                        class_dict["Ensembl_Annotated"].append(gene_id)
                    elif "Ensembl_Annotated" not in class_dict:
                        class_dict.update({"Ensembl_Annotated":[gene_id]})
                elif gene_id.startswith("ENSGACG") and transcript_id == "novel":
                    if "Novel" in class_dict:
                        class_dict["Novel"].append(gene_id)
                    elif "Novel" not in class_dict:
                        class_dict.update({"Novel":[gene_id]})
                elif gene_id.startswith("novelGene") and transcript_id == "novel":
                    if "Novel" in class_dict:
                        class_dict["Novel"].append(gene_id)
                    elif "Novel" not in class_dict:
                        class_dict.update({"Novel":[gene_id]})
    return class_dict


#read in ensembl annotated go terms
#returns dictionary with key == gene id and value == list of go terms for that gene id
def read_ensembl_go_terms():
    ensembl_file = sys.argv[2]
    ensembl_dict = {}
    with open(ensembl_file, 'r') as ensembl:
        for line in ensembl:
            new_line = line.split()
            gene_id = new_line[0]
            go_terms = new_line[1].split(",")
            ensembl_dict.update({gene_id:go_terms})
    return ensembl_dict


#read in novel annotated go terms
#returns dictionary with key == gene id and value == list of go terms for that gene id
def read_novel_go_terms():
    novel_file = sys.argv[3]
    novel_dict = {}
    with open(novel_file, 'r') as novel:
        for line in novel:
            new_line = line.split()
            gene_id = new_line[0]
            go_terms = new_line[1].split(",")
            novel_dict.update({gene_id:go_terms})
    return novel_dict

#pull go terms for all isoforms in classification file
#returns dictionary with key == gene id and value == list of GO terms for gene
def combine():
    class_dict = read_class_file()
    ensembl_dict = read_ensembl_go_terms()
    novel_dict = read_novel_go_terms()
    annotated_isoforms = class_dict["Ensembl_Annotated"]
    novel_isoforms = class_dict["Novel"]
    combined_dict = {}
    final_dict = {}
    for ensembl_gene in ensembl_dict:
        single_gene = ensembl_dict[ensembl_gene]
        if ensembl_gene in annotated_isoforms:
            if ensembl_gene in combined_dict:
                combined_dict[ensembl_gene].append(single_gene)
            elif ensembl_gene not in combined_dict:
                combined_dict.update({ensembl_gene:[single_gene]})
    for novel_gene in novel_dict:
        single_gene = novel_dict[novel_gene]
        if novel_gene in novel_isoforms:
            if novel_gene in combined_dict:
                combined_dict[novel_gene].append(single_gene)
            elif novel_gene not in combined_dict:
                combined_dict.update({novel_gene:[single_gene]})
    for key in combined_dict:
        single_key = combined_dict[key][0]
        set_key = list(set(single_key))
        final_dict.update({key:set_key})
    return final_dict



#write output
def write():
    final_dict = combine()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for gene in final_dict:
            single_gene = final_dict[gene]
            go_terms = ",".join(single_gene)
            final = gene + "\t" + go_terms + "\n"
            out.write(final)


write()

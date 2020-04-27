#combined annotated ensembl gene go terms and novel gene GO terms
#this script will be used by Randomize_GO_terms_isoseq.py to create random distribution
#to run script: python3 Combined_Ensembl_Novel_GO_terms.py <ensembl table as csv> <novel genes go terms file for combined sexes analysis> <combined sexes classification file> <output file; format: Gene.ID \t List of GO terms
#Author: Alice Naftaly, March 2020

import sys

#first pull genes and go terms from ensembl table
#only need gene id and go terms
#returns dictionary with key = gene id and value == list of go terms
def read_ensembl_table():
    ensembl_table_file = sys.argv[1]
    go_term_dict = {}
    final_go_terms = {}
    with open(ensembl_table_file, 'r') as ensembl_table:
        for line in ensembl_table:
            new_line = line.split(",")
            gene_id = new_line[0]
            go_term = new_line[len(new_line)-1].strip("\n")
            if gene_id in go_term_dict:
                go_term_dict[gene_id].append(go_term)
            elif gene_id not in go_term_dict:
                go_term_dict.update({gene_id:[go_term]})
    #remove genes with no go terms or empty go term values
    for key in go_term_dict:
        single_key = go_term_dict[key]
        list_of_go_terms = []
        for value in single_key:
            if value.startswith("GO:"):
                list_of_go_terms.append(value)
        if len(list_of_go_terms) > 0:
            final_go_terms.update({key:list_of_go_terms})
    return final_go_terms

#need to add the novel genes detected to the ensembl annotated genes
#read GO terms for novel genes
#returns dictionary with key == isoform id and value == [list of go terms]
def read_novel_GO_terms():
    go_terms_file = sys.argv[2]
    go_term_dict = {}
    with open(go_terms_file, 'r') as go_terms:
        for line in go_terms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                single_go_term = new_line[2]
                if isoform_id in go_term_dict:
                    go_term_dict[isoform_id].append(single_go_term)
                elif isoform_id not in go_term_dict:
                    go_term_dict.update({isoform_id:[single_go_term]})
    return go_term_dict

#read combined sexes classification file to get appropriate gene ids
#returns gene dict where key = isoform id and value == gene id
def pull_combined_sexes_genes():
    class_file = sys.argv[3]
    gene_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                gene_id = new_line[6]
                if gene_id.startswith("novel"):
                    gene_dict.update({isoform_id:gene_id})
    return gene_dict


#combine genes and GO terms for novel genes
#returns a dictionary with key = gene id and value == list of go terms
#note: the number of genes this pulls is relatively low since these are the novel genes; there are GO terms for novel transcripts that match to known ensembl genes, but there are far fewer GO terms for completely novel genes.
def combine_novel_genes_GO_terms():
    novel_go_terms = read_novel_GO_terms()
    gene_ids = pull_combined_sexes_genes()
    novel_gene_go_terms = {}
    for isoform in gene_ids:
        single_gene = gene_ids[isoform]
        if isoform in novel_go_terms:
            single_isoform = novel_go_terms[isoform]
            novel_gene_go_terms.update({single_gene:single_isoform})
    return novel_gene_go_terms

#need to combine go terms from ensembl annotations and novel genes for combined sexes
#returns one dictionary with key == gene and value == list of go terms for that gene
def combine_ensembl_novel_genes():
    ensembl_genes = read_ensembl_table()
    novel_genes = combine_novel_genes_GO_terms()
    combined_final_dict = {}
    combined_final_dict.update(ensembl_genes)
    combined_final_dict.update(novel_genes)
    return combined_final_dict

#write combined dictionary to output file
#final format should be gene id \t list of go terms as list separated by commas
def write():
    final_dict = combine_ensembl_novel_genes()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for gene in final_dict:
            single_gene = final_dict[gene]
            go_term_list = ",".join(single_gene)
            final = "%s\t%s\n" % (str(gene), go_term_list)
            out.write(final)
write()

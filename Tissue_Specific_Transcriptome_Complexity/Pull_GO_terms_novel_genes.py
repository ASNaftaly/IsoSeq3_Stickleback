#Pull GO terms for ensembl novel genes
#need input files that list isoform ids 1 per line, file from interproscan with isoform\tcategory\tGO.term
#will pull all go terms from each isoform in file, if no go terms are available, will skip that gene
#final file should have the final format:
#Gene \t list of GO terms separated by ,
#to run script: Pull_GO_terms_novel_genes.py <converted classification file> <combined sexes classification file> <all tissues novel genes go terms file from interproscan and BLAST2GO> <output file; format = Gene.ID\tGo.Terms>
#Author: Alice Naftaly, March 2020

import sys

#reads in isoform ids and gene ids from file
#will only keep novel gene; reads from classification file with convertd isoform ids
#returns isoform ids in a list
def read_isoform_ids_converted():
    class_file = sys.argv[1]
    all_gene_ids = []
    isoform_list = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                #converted isoform id
                isoform_id = new_line[0]
                #gene_id
                gene_id = new_line[7]
                if gene_id.startswith("novel"):
                    isoform_list.append(isoform_id)
    return isoform_list

#read combined sexes classification file to get appropriate gene ids
#returns gene dict where key = isoform id and value == gene id
def pull_combined_sexes_genes():
    class_file = sys.argv[2]
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

#combine gene ids and isoform ids
#returns dictionary with key == isoform id and value == gene id
def combine():
    isoform_ids = read_isoform_ids_converted()
    combined_genes = pull_combined_sexes_genes()
    combined_dict = {}
    for isoform in isoform_ids:
        if isoform in combined_genes:
            single_gene = combined_genes[isoform]
            combined_dict.update({isoform:single_gene})
    return combined_dict

#read GO terms for novel genes
#returns dictionary with key == isoform id and value == [list of go terms]
def read_GO_terms():
    go_terms_file = sys.argv[3]
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


#match genes to go terms
#returns dictionary with key == gene id and value
def match_genes_go_terms():
    go_terms = read_GO_terms()
    isoform_dict = combine()
    matched_isoforms = {}
    for isoform in isoform_dict:
        if isoform in go_terms:
            single_isoform = go_terms[isoform]
            gene_id = isoform_dict[isoform]
            if gene_id in matched_isoforms:
                matched_isoforms[gene_id].append(single_isoform)
            elif gene_id not in matched_isoforms:
                matched_isoforms.update({gene_id:[single_isoform]})
    return matched_isoforms


#write matched genes and GO terms to a file
#final output = Gene.Id\tGO.Terms (in a list separated by commas)
def write():
    matched_genes = match_genes_go_terms()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for key in matched_genes:
            single_key = matched_genes[key][0]
            go_terms = ",".join(single_key)
            final = "%s\t%s\n" % (str(key), go_terms)
            out.write(final)

write()

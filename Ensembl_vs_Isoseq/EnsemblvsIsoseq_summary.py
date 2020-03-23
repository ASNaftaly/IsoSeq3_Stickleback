#ensembl vs combined sexes analysis
#want to count the number of isoforms per gene by splitting into the following categories: 1 isoform, 2-3 isoforms, 4-5 isoforms, more than 6 isoforms
#also want this to print out the number of coding transcripts vs noncoding transcripts
#also print out structural categories to compare coding vs noncoding
#to run script: EnsemblvsIsoseq_summary.py <ensembl classification file> <combined sexes classification file with duplicates removed> >> <log of output; print statements>


import sys

#read in ensembl classificaiton file
#returns dictionary with key == gene id and value == [transcript id, splice type, coding]
def read_ensembl_class():
    ensembl_file = sys.argv[1]
    ensembl_dict = {}
    with open(ensembl_file, 'r') as ensembl:
        for line in ensembl:
            if line.startswith("EN"):
                new_line = line.split("\t")
                gene_id = new_line[6]
                transcript_id = new_line[7]
                splice_type = new_line[5]
                coding = new_line[27]
                final = [transcript_id, splice_type, coding]
                if gene_id in ensembl_dict:
                    ensembl_dict[gene_id].append(final)
                elif gene_id not in ensembl_dict:
                    ensembl_dict.update({gene_id:[final]})
    return ensembl_dict


#read in dups removed combined sexes classification file
def read_isoseq_classification():
    isoseq_file = sys.argv[2]
    isoseq_dict = {}
    with open(isoseq_file, 'r') as isoseq:
        for line in isoseq:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                splice_type = new_line[5]
                coding = new_line[27]
                final = [isoform_id, transcript_id, splice_type, coding]
                if gene_id in isoseq_dict:
                    isoseq_dict[gene_id].append(final)
                elif gene_id not in isoseq_dict:
                    isoseq_dict.update({gene_id:[final]})
    return isoseq_dict

#print ensembl stats
def ensembl_counts():
    ensembl_dict = read_ensembl_class()
    isoform_list = []
    for gene in ensembl_dict:
        single_gene = ensembl_dict[gene]
        for value in single_gene:
            isoform_list.append(value[0])
    num_genes = len(ensembl_dict)
    num_isoforms = len(set(isoform_list))
    print("Number of Genes in Ensembl Annotations, build 97")
    print(num_genes)
    print("Number of Isoforms in Ensembl Annotations, build 97")
    print(num_isoforms)

#print number of non duplicate genes and isoforms:
def isoseq_counts():
    isoseq_dict = read_isoseq_classification()
    isoform_list = []
    for gene in isoseq_dict:
        single_gene = isoseq_dict[gene]
        for value in single_gene:
            isoform_list.append(value[0])
    num_genes = len(isoseq_dict)
    num_isoforms = len(set(isoform_list))
    print("Number of Genes in Isoseq Combined Sexes Analysis after Duplicates removed")
    print(num_genes)
    print("Number of Isoforms in Isoseq Combined Sexes after Duplicates removed")
    print(num_isoforms)


#how many isoforms are coding
def overall_coding_counts():
    isoseq_dict = read_isoseq_classification()
    coding_list = []
    non_coding_list = []
    for gene in isoseq_dict:
        single_gene = isoseq_dict[gene]
        for value in single_gene:
            code = value[3]
            if code == "coding":
                coding_list.append(value[0])
            elif code == "non_coding":
                non_coding_list.append(value[0])
    num_coding = len(coding_list)
    num_non_coding = len(non_coding_list)
    print("Number of isoforms that are coding")
    print(num_coding)
    print("Number of isoforms that are not coding")
    print(num_non_coding)


#count number of isfoorms per gene for both ensembl and isoseq files
def num_isoforms_per_gene():
    ensembl_dict = read_ensembl_class()
    isoseq_dict = read_isoseq_classification()
    ensembl_1_isoform = []
    isoseq_1_isoform = []
    ensembl_2_3_isoforms = []
    isoseq_2_3_isoforms = []
    ensembl_4_5_isoforms = []
    isoseq_4_5_isoforms = []
    ensembl_more_than_6_isoforms = []
    isoseq_more_than_6_isoforms = []
    for gene_1 in ensembl_dict:
        single_gene_1 = ensembl_dict[gene_1]
        if len(single_gene_1) == 1:
            ensembl_1_isoform.append(gene_1)
        elif 2 <= len(single_gene_1) <= 3:
            ensembl_2_3_isoforms.append(gene_1)
        elif 4 <= len(single_gene_1) <= 5:
            ensembl_4_5_isoforms.append(gene_1)
        elif len(single_gene_1) >= 6:
            ensembl_more_than_6_isoforms.append(gene_1)
    for gene_2 in isoseq_dict:
        single_gene_2 = isoseq_dict[gene_2]
        if len(single_gene_2) == 1:
            isoseq_1_isoform.append(gene_2)
        elif 2 <= len(single_gene_2) <= 3:
            isoseq_2_3_isoforms.append(gene_2)
        elif 4 <= len(single_gene_2) <= 5:
            isoseq_4_5_isoforms.append(gene_2)
        elif len(single_gene_2) >= 6:
            isoseq_more_than_6_isoforms.append(gene_2)
    print("Number of transcripts per Ensembl gene \n")
    print("Number of Genes with 1 isoform:")
    print(len(ensembl_1_isoform))
    print("Number of Genes with 2-3 isoforms")
    print(len(ensembl_2_3_isoforms))
    print("Number of Genes with 4-5 isoforms")
    print(len(ensembl_4_5_isoforms))
    print("Number of Genes with 6 or more isoforms")
    print(len(ensembl_more_than_6_isoforms))
    print("Number of transcripts per Isoseq gene \n")
    print("Number of Genes with 1 isoform:")
    print(len(isoseq_1_isoform))
    print("Number of Genes with 2-3 isoforms")
    print(len(isoseq_2_3_isoforms))
    print("Number of Genes with 4-5 isoforms")
    print(len(isoseq_4_5_isoforms))
    print("Number of Genes with 6 or more isoforms")
    print(len(isoseq_more_than_6_isoforms))

#splice type counts along with coding vs noncoding for isoseq isoforms
def splice_type_counts():
    isoseq_dict = read_isoseq_classification()
    splice_type_dict = {}
    for gene in isoseq_dict:
        single_gene = isoseq_dict[gene]
        for value in single_gene:
            splice_type = value[2]
            coding_potential = value[3]
            if splice_type in splice_type_dict:
                splice_type_dict[splice_type].append(coding_potential)
            elif splice_type not in splice_type_dict:
                splice_type_dict.update({splice_type:[coding_potential]})
    print("Number of isoforms in each splice type and coding/noncoding split\n")
    for key in splice_type_dict:
        coding_counts = 0
        non_coding_counts = 0
        print(key)
        single_key = splice_type_dict[key]
        print("Number of isoforms in splice type")
        print(len(single_key))
        for value in single_key:
            if value == "coding":
                coding_counts += 1
            elif value == "non_coding":
                non_coding_counts += 1
        print("Number of coding/noncoding counts")
        print("Number Coding")
        print(coding_counts)
        print("Number NonCoding")
        print(non_coding_counts)


#call all functions:
def call():
    ensembl_stats = ensembl_counts()
    general_counts = isoseq_counts()
    coding_counts = overall_coding_counts()
    isoform_counts = num_isoforms_per_gene()
    splice_counts = splice_type_counts()

call()

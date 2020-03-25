#pull novel genes and isoforms for Novel isoforms portion of isoseq paper
#will read in combined sexes classification file without duplicates
#to run script: python3 Novel_Isoform_Breakdown.py <combined sexes classification file>
#

import sys

#read in classification file
#returns dictionary with key = isoform id and value == [gene id, transcript id, splice type, exon type, canonical, coding potential]
def read_combined_classification():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform_id = new_line[0]
                splice_type = new_line[5]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                exon_type = new_line[12]
                canonical = new_line[14]
                coding_potential = new_line[27]
                final = [gene_id, transcript_id, splice_type, exon_type, canonical, coding_potential]
                class_dict.update({isoform_id:final})
    return class_dict


#count novel isoforms
#this counts both novel genes and novel transcripts
def total_count_novel_isoforms():
    class_dict = read_combined_classification()
    annotated_counts = []
    novel_counts = []
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        #if gene is novel
        if single_isoform[0].startswith("novel"):
            novel_counts.append(isoform)
        elif single_isoform[0].startswith("EN") and single_isoform[1].startswith("novel"):
            novel_counts.append(isoform)
        elif single_isoform[0].startswith("EN") and single_isoform[1].startswith("EN"):
            annotated_counts.append(isoform)
    print("Number of novel isoforms")
    print(len(novel_counts))
    print("Number of annotated isoforms")
    print(len(annotated_counts))
    print("\n\n\n")



#Novel genes break down
def novel_genes_breakdown():
    class_dict = read_combined_classification()
    novel_genes = {}
    novel_genes_counts = 0
    #creates a new dictionary with key == novel gene ids and value == [isoform, transcript id, splice type, exon type, canonical, coding potential]
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        if single_isoform[0].startswith("novel"):
            novel_genes_counts += 1
            final = [isoform, single_isoform[1], single_isoform[2], single_isoform[3], single_isoform[4], single_isoform[5]]
            if single_isoform[0] in novel_genes:
                novel_genes[single_isoform[0]].append(final)
            elif single_isoform[0] not in novel_genes:
                novel_genes.update({single_isoform[0]:[final]})
    print("Total number of novel genes")
    print(novel_genes_counts)
    print("\n")
    coding_dict = {}
    noncoding_dict = {}
    #making 2 new dictionaries for coding and noncoding potential where each value == [isoform, gene, splice type, exon type, canonical]
    for gene in novel_genes:
        single_gene = novel_genes[gene]
        for value in single_gene:
            isoform = value[0]
            splice_type = value[2]
            exon_type = value[3]
            canonical = value[4]
            coding_potential = value[5]
            final = [isoform, gene, splice_type, exon_type, canonical]
            if coding_potential == "coding":
                if isoform in coding_dict:
                    coding_dict[isoform].append(final)
                elif isoform not in coding_dict:
                    coding_dict.update({isoform:[final]})
            elif coding_potential == "non_coding":
                if isoform in noncoding_dict:
                    noncoding_dict[isoform].append(final)
                elif isoform not in noncoding_dict:
                    noncoding_dict.update({isoform:[final]})
    print("Number of coding novel genes")
    print(len(coding_dict))
    print("Number of noncoding novel genes")
    print(len(noncoding_dict))
    print("\n")
    coding_splice_dict = {}
    noncoding_splice_dict = {}
    coding_exon_dict = {}
    noncoding_exon_dict = {}
    coding_canonical_counts = 0
    coding_noncanonical_counts = 0
    noncoding_canonical_counts = 0
    noncoding_noncanonical_counts = 0
    #breakdown of coding and noncoding genes based on splice type and exon type and canonical
    for coding_gene in coding_dict:
        single_coding_gene = coding_dict[coding_gene]
        for val in single_coding_gene:
            splice_type = val[2]
            exon_type = val[3]
            can_vs_noncan = val[4]
            if splice_type in coding_splice_dict:
                coding_splice_dict[splice_type].append("1")
            elif splice_type not in coding_splice_dict:
                coding_splice_dict.update({splice_type:["1"]})
            if exon_type in coding_exon_dict:
                coding_exon_dict[exon_type].append("1")
            elif exon_type not in coding_exon_dict:
                coding_exon_dict.update({exon_type:["1"]})
            if can_vs_noncan == "canonical":
                coding_canonical_counts += 1
            elif can_vs_noncan == "non_canonical":
                coding_noncanonical_counts += 1
    for non_coding_gene in noncoding_dict:
        single_noncoding_gene = noncoding_dict[non_coding_gene]
        for v in single_noncoding_gene:
            splice_type = v[2]
            exon_type = v[3]
            can_vs_noncan = v[4]
            if splice_type in noncoding_splice_dict:
                noncoding_splice_dict[splice_type].append("1")
            elif splice_type not in noncoding_splice_dict:
                noncoding_splice_dict.update({splice_type:["1"]})
            if exon_type in noncoding_exon_dict:
                noncoding_exon_dict[exon_type].append("1")
            elif exon_type not in noncoding_exon_dict:
                noncoding_exon_dict.update({exon_type:["1"]})
            if can_vs_noncan == "canonical":
                noncoding_canonical_counts += 1
            elif can_vs_noncan == "non_canonical":
                noncoding_noncanonical_counts += 1
    #print out coding genes breakdown
    print("Number of splice types for novel coding genes")
    for key_1 in coding_splice_dict:
        print(key_1)
        print(len(coding_splice_dict[key_1]))
    print("Number of exon types for novel coding genes")
    for key_2 in coding_exon_dict:
        print(key_2)
        print(len(coding_exon_dict[key_2]))
    print("Number of canonical classification for novel coding genes")
    print("Canonical")
    print(coding_canonical_counts)
    print("Non_canonical")
    print(coding_noncanonical_counts)
    print("\n")
    #print out noncoding genes breakdown
    print("Number of splice types for novel noncoding genes")
    for key_3 in noncoding_splice_dict:
        print(key_3)
        print(len(noncoding_splice_dict[key_3]))
    print("Number of exon types for novel noncoding genes")
    for key_4 in noncoding_exon_dict:
        print(key_4)
        print(len(noncoding_exon_dict[key_4]))
    print("Number of canonical classification for novel noncoding genes")
    print("Canonical")
    print(noncoding_canonical_counts)
    print("Non_canonical")
    print(noncoding_noncanonical_counts)
    print("\n")


#Novel transcripts break down
def novel_transcripts_breakdown():
    class_dict = read_combined_classification()
    novel_transcripts = {}
    #creates a new dictionary with key == novel transcripts ids and value == [isoform, transcript id, splice type, exon type, canonical, coding potential]
    novel_transcripts_counts = 0
    for isoform in class_dict:
        single_isoform = class_dict[isoform]
        if single_isoform[0].startswith("EN") and single_isoform[1].startswith("novel"):
            final = [isoform, single_isoform[1], single_isoform[2], single_isoform[3], single_isoform[4], single_isoform[5]]
            novel_transcripts_counts += 1
            if single_isoform[1] in novel_transcripts:
                novel_transcripts[single_isoform[1]].append(final)
            elif single_isoform[1] not in novel_transcripts:
                novel_transcripts.update({single_isoform[1]:[final]})
    print("Total number of novel transcripts")
    print(novel_transcripts_counts)
    print("\n")
    coding_dict = {}
    coding_counts = 0
    noncoding_dict = {}
    non_coding_counts = 0
    #making 2 new dictionaries for coding and noncoding potential where each value == [isoform, gene, splice type, exon type, canonical]
    for gene in novel_transcripts:
        single_transcript = novel_transcripts[gene]
        for value in single_transcript:
            isoform = value[0]
            splice_type = value[2]
            exon_type = value[3]
            canonical = value[4]
            coding_potential = value[5]
            final = [isoform, splice_type, exon_type, canonical]
            if coding_potential == "coding":
                coding_counts += 1
                if isoform in coding_dict:
                    coding_dict[isoform].append(final)
                elif isoform not in coding_dict:
                    coding_dict.update({isoform:[final]})
            elif coding_potential == "non_coding":
                non_coding_counts += 1
                if isoform in noncoding_dict:
                    noncoding_dict[isoform].append(final)
                elif isoform not in noncoding_dict:
                    noncoding_dict.update({isoform:[final]})
    print("Number of coding novel transcripts")
    print(coding_counts)
    print("Number of noncoding novel transcripts")
    print(non_coding_counts)
    print("\n")
    coding_splice_dict = {}
    noncoding_splice_dict = {}
    coding_exon_dict = {}
    noncoding_exon_dict = {}
    coding_canonical_counts = 0
    coding_noncanonical_counts = 0
    noncoding_canonical_counts = 0
    noncoding_noncanonical_counts = 0
    #breakdown of coding and noncoding transcripts based on splice type and exon type and canonical
    for coding_gene in coding_dict:
        single_coding_gene = coding_dict[coding_gene]
        for val in single_coding_gene:
            splice_type = val[1]
            exon_type = val[2]
            can_vs_noncan = val[3]
            if splice_type in coding_splice_dict:
                coding_splice_dict[splice_type].append("1")
            elif splice_type not in coding_splice_dict:
                coding_splice_dict.update({splice_type:["1"]})
            if exon_type in coding_exon_dict:
                coding_exon_dict[exon_type].append("1")
            elif exon_type not in coding_exon_dict:
                coding_exon_dict.update({exon_type:["1"]})
            if can_vs_noncan == "canonical":
                coding_canonical_counts += 1
            elif can_vs_noncan == "non_canonical":
                coding_noncanonical_counts += 1
    for non_coding_gene in noncoding_dict:
        single_noncoding_gene = noncoding_dict[non_coding_gene]
        for v in single_noncoding_gene:
            splice_type = v[1]
            exon_type = v[2]
            can_vs_noncan = v[3]
            if splice_type in noncoding_splice_dict:
                noncoding_splice_dict[splice_type].append("1")
            elif splice_type not in noncoding_splice_dict:
                noncoding_splice_dict.update({splice_type:["1"]})
            if exon_type in noncoding_exon_dict:
                noncoding_exon_dict[exon_type].append("1")
            elif exon_type not in noncoding_exon_dict:
                noncoding_exon_dict.update({exon_type:["1"]})
            if can_vs_noncan == "canonical":
                noncoding_canonical_counts += 1
            elif can_vs_noncan == "non_canonical":
                noncoding_noncanonical_counts += 1
    #print out coding genes breakdown
    print("Number of splice types for novel coding transcripts")
    for key_1 in coding_splice_dict:
        print(key_1)
        print(len(coding_splice_dict[key_1]))
    print("Number of exon types for novel coding transcripts")
    for key_2 in coding_exon_dict:
        print(key_2)
        print(len(coding_exon_dict[key_2]))
    print("Number of canonical classification for novel coding transcripts")
    print("Canonical")
    print(coding_canonical_counts)
    print("Non_canonical")
    print(coding_noncanonical_counts)
    print("\n")
    #print out noncoding genes breakdown
    print("Number of splice types for novel noncoding transcripts")
    for key_3 in noncoding_splice_dict:
        print(key_3)
        print(len(noncoding_splice_dict[key_3]))
    print("Number of exon types for novel noncoding transcripts")
    for key_4 in noncoding_exon_dict:
        print(key_4)
        print(len(noncoding_exon_dict[key_4]))
    print("Number of canonical classification for novel noncoding transcripts")
    print("Canonical")
    print(noncoding_canonical_counts)
    print("Non_canonical")
    print(noncoding_noncanonical_counts)
    print("\n")



#call all functions
def call():
    total_novel_counts = total_count_novel_isoforms()
    total_novel_genes = novel_genes_breakdown()
    total_novel_transcripts = novel_transcripts_breakdown()

call()

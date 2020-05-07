#pulling sequences for exon trios to determine exon usage
#will analyze 1 chr or scaffold at a time
#- strand will need the reverse complement for the sequence (this will be a separate script)
#to run script: python3 PullSeqs_ExonUsage.py <chromosome or scaffold number> <Ensembl_adjusted fasta file> <gene IDs with chromosome numbers file> <exon trios file> <output file, trios fasta sequences> <output file, trios lengths for records>
#Author: Alice Naftaly, May 2020

import sys

#read in whole fasta file
#returns sequence of one chromosome or scaffold as list
def read_fasta():
    chr_num = sys.argv[1]
    fasta_file = sys.argv[2]
    fasta_dict = {}
    final_single_seq = []
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                header_line = line.strip("\n")
                final_header = header_line.strip(">")
            else:
                if final_header in fasta_dict:
                    fasta_dict[final_header].append(line.strip("\n"))
                elif final_header not in fasta_dict:
                    fasta_dict.update({final_header:[line.strip("\n")]})
    single_chr = fasta_dict[chr_num]
    for seq in single_chr:
        final_single_seq += seq
    return final_single_seq


#read gene ids and chromosome numbers file to figure out which genes need to be examined
#returns a list of genes that are on the chromosome
def read_geneIDs_chrNums():
    chr_num = sys.argv[1]
    input_file = sys.argv[3]
    gene_list = []
    with open(input_file, 'r') as info:
        for line in info:
            if line.startswith("Gene.ID"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                chr_number = new_line[1]
                if chr_num == chr_number:
                    gene_list.append(gene_id)
    return gene_list


#read exon trio file and pull exon trios for the genes to be examined on a single chromosome/scaffold
#returns dictionary with key == geneID_+/- strand and value = all exon trios for that gene
def read_exon_trios():
    exon_file = sys.argv[4]
    genes = read_geneIDs_chrNums()
    exon_trios_dict = {}
    with open(exon_file, 'r') as exons:
        for line in exons:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                strand = new_line[1]
                #exon trio = [exon 1 start, exon 1 end, exon 2 start, exon 2 end, exon 3 start, exon 3 end]
                exon_trio = [new_line[2], new_line[3],new_line[4],new_line[5],new_line[6], new_line[7]]
                if gene_id in genes:
                    new_id = gene_id + "_" + strand
                    if new_id in exon_trios_dict:
                        exon_trios_dict[new_id].append(exon_trio)
                    elif new_id not in exon_trios_dict:
                        exon_trios_dict.update({new_id:[exon_trio]})
    return exon_trios_dict


#pulling sequences for each trio
#returns dictionary with key = gene.id_chr.num_exon.trio.num_scenario and value == [exon junction sequence as a list, exon length (1 or 2), exon length (2 or 3)]
def pull_seqs():
    exon_trios = read_exon_trios()
    chr_num = sys.argv[1]
    sequences = read_fasta()
    exon_trio_seqs_dict = {}
    for gene in exon_trios:
        single_gene = exon_trios[gene]
        #this is to count the exon trios per gene (in output file)
        exon_trio_num = 1
        for exon_trio in single_gene:
            exon_1 = [int(exon_trio[0]), int(exon_trio[1])]
            exon_2 = [int(exon_trio[2]), int(exon_trio[3])]
            exon_3 = [int(exon_trio[4]), int(exon_trio[5])]
            exon_1_size = abs(exon_1[0] - exon_1[1])
            exon_2_size = abs(exon_2[0] - exon_2[1])
            exon_3_size = abs(exon_3[0] - exon_3[1])
            #each sequence needs a unique identifier: gene.id_chr.num_exon.trio.num_scenario
            #there will be 2 scenarios (no alternative splicing and with alternative splicing)
            #no alternative splicing will have 2 sets of sequences (exon 1-2 and exon 2-3)
            #the labels for these scenarios will be: no alternative splicing = no_as1 (exon 1-2) and no_as2 (exon 2-3) and alternative splicing will be as1
            #first scenario is that all 3 exons are used (no alternative splicing)
            #this means I need to pull 75bp from the end of the first exon, 75bp from the beginning of exon 2 as one exon pair
            #then need to pull 75bp from the end of the second exon and 75bp from the beginning of exon 3
            #if the exons are not 75bp long, will pull the entire exon
            #will write a separate file with the lengths of each exon trio
            #if both exons are greater than 75bp long
            if exon_1_size > 75 and exon_2_size > 75:
                #exon 1 end
                #determining direction
                #+ strand
                if exon_1[0] < exon_1[1]:
                    upstream_75bp_1 = exon_1[1] - 75
                    exon_1_start_pull = upstream_75bp_1 - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    upstream_75bp_1 = exon_1[1] + 75
                    exon_1_start_pull = upstream_75bp_1 + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 2 begining
                #+ strand
                if exon_2[0] < exon_2[1]:
                    downstream_75bp_2 = exon_2[0] + 75
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = downstream_75bp_2 - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    downstream_75bp_2 = exon_2[0] - 75
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = downstream_75bp_2 + 1
            #if exon 1 is less than 75bp long
            elif exon_1_size <= 75 and exon_2_size > 75:
                #exon 1 end
                #+ strand
                if exon_1[0] < exon_1[1]:
                    exon_1_start_pull = exon_1[0] - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    exon_1_start_pull = exon_1[0] + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 2 begining
                #+ strand
                if exon_2[0] < exon_2[1]:
                    downstream_75bp_2 = exon_2[0] + 75
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = downstream_75bp_2 - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    downstream_75bp_2 = exon_2[0] - 75
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = downstream_75bp_2 + 1
            #if exon 2 is less than 75bp long
            elif exon_1_size > 75 and exon_2_size <= 75:
                #exon 1 end
                #determining direction
                #+ strand
                if exon_1[0] < exon_1[1]:
                    upstream_75bp_1 = exon_1[1] - 75
                    exon_1_start_pull = upstream_75bp_1 - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    upstream_75bp_1 = exon_1[1] + 75
                    exon_1_start_pull = upstream_75bp_1 + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 2 begining
                #+ strand
                if exon_2[0] < exon_2[1]:
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = exon_2[1] + 1
            elif exon_1_size <= 75 and exon_2_size <= 75:
                #exon 1 end
                #+ strand
                if exon_1[0] < exon_1[1]:
                    exon_1_start_pull = exon_1[0] - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    exon_1_start_pull = exon_1[0] + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 2 begining
                #+ strand
                if exon_2[0] < exon_2[1]:
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = exon_2[1] + 1
            #pulling sequences for exon junction 1-2
            #+ strand
            if exon_1_start_pull < exon_1_end_pull:
                no_as1_exon_1_seq = sequences[exon_1_start_pull:exon_1_end_pull]
                no_as1_exon_2_seq = sequences[exon_2_start_pull:exon_2_end_pull]
            elif exon_1_start_pull > exon_1_end_pull:
                no_as1_exon_1_seq = sequences[exon_1_end_pull:exon_1_start_pull]
                no_as1_exon_2_seq = sequences[exon_2_end_pull:exon_2_start_pull]
            no_as1_sequence = no_as1_exon_1_seq + no_as1_exon_2_seq
            header_values = [gene, chr_num, str(exon_trio_num),"no_as1"]
            no_as1_header = ".".join(header_values)
            dict_value = [no_as1_sequence, exon_1_size, exon_2_size]
            exon_trio_seqs_dict.update({no_as1_header:dict_value})
            #exon junction 2-3
            if exon_2_size > 75 and exon_3_size > 75:
                #exon 2 end
                #determining direction
                #+ strand
                if exon_2[0] < exon_2[1]:
                    upstream_75bp_2 = exon_2[1] - 75
                    exon_2_start_pull = upstream_75bp_2 - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    upstream_75bp_2 = exon_2[1] + 75
                    exon_2_start_pull = upstream_75bp_2 + 1
                    exon_2_end_pull = exon_2[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    downstream_75bp_3 = exon_3[0] + 75
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = downstream_75bp_3 - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    downstream_75bp_3 = exon_3[0] - 75
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = downstream_75bp_3 + 1
            #if exon 2 is less than 75bp long
            elif exon_2_size <= 75 and exon_3_size > 75:
                #exon 2 end
                #+ strand
                if exon_2[0] < exon_2[1]:
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = exon_2[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    downstream_75bp_3 = exon_3[0] + 75
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = downstream_75bp_3 - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    downstream_75bp_3 = exon_3[0] - 75
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = downstream_75bp_3 + 1
            #if exon 3 is less than 75bp long
            elif exon_2_size > 75 and exon_3_size <= 75:
                #exon 2 end
                #determining direction
                #+ strand
                if exon_2[0] < exon_2[1]:
                    upstream_75bp_2 = exon_2[1] - 75
                    exon_2_start_pull = upstream_75bp_2 - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    upstream_75bp_2 = exon_2[1] + 75
                    exon_2_start_pull = upstream_75bp_2 + 1
                    exon_2_end_pull = exon_2[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = exon_3[1] - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = exon_3[1] + 1
            elif exon_2_size <= 75 and exon_3_size <= 75:
                #exon 3 end
                #+ strand
                if exon_2[0] < exon_2[1]:
                    exon_2_start_pull = exon_2[0] - 1
                    exon_2_end_pull = exon_2[1] - 1
                #- strand
                elif exon_2[0] > exon_2[1]:
                    exon_2_start_pull = exon_2[0] + 1
                    exon_2_end_pull = exon_2[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = exon_3[1] - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = exon_3[1] + 1
            #pulling sequences for exon junction 2-3
            #+ strand
            if exon_2_start_pull < exon_2_end_pull:
                no_as2_exon_2_seq = sequences[exon_2_start_pull:exon_2_end_pull]
                no_as2_exon_3_seq = sequences[exon_3_start_pull:exon_3_end_pull]
            elif exon_2_start_pull > exon_2_end_pull:
                no_as2_exon_2_seq = sequences[exon_2_end_pull:exon_1_start_pull]
                no_as2_exon_3_seq = sequences[exon_3_end_pull:exon_2_start_pull]
            no_as2_sequence = no_as2_exon_2_seq + no_as2_exon_3_seq
            header_values = [gene, chr_num, str(exon_trio_num),"no_as2"]
            no_as2_header = ".".join(header_values)
            dict_value = [no_as2_sequence, exon_2_size, exon_3_size]
            exon_trio_seqs_dict.update({no_as2_header:dict_value})
            #alternative splicing scenario (exons 1-3)
            if exon_1_size > 75 and exon_3_size > 75:
                #exon 1 end
                #determining direction
                #+ strand
                if exon_1[0] < exon_1[1]:
                    upstream_75bp_1 = exon_1[1] - 75
                    exon_1_start_pull = upstream_75bp_1 - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    upstream_75bp_1 = exon_1[1] + 75
                    exon_1_start_pull = upstream_75bp_1 + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    downstream_75bp_3 = exon_3[0] + 75
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = downstream_75bp_3 - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    downstream_75bp_3 = exon_3[0] - 75
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = downstream_75bp_3 + 1
            #if exon 1 is less than 75bp long
            elif exon_1_size <= 75 and exon_3_size > 75:
                #exon 1 end
                #+ strand
                if exon_1[0] < exon_1[1]:
                    exon_1_start_pull = exon_1[0] - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    exon_1_start_pull = exon_1[0] + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    downstream_75bp_3 = exon_3[0] + 75
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = downstream_75bp_3 - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    downstream_75bp_3 = exon_3[0] - 75
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = downstream_75bp_3 + 1
            #if exon 3 is less than 75bp long
            elif exon_1_size > 75 and exon_3_size <= 75:
                #exon 1 end
                #determining direction
                #+ strand
                if exon_1[0] < exon_1[1]:
                    upstream_75bp_1 = exon_1[1] - 75
                    exon_1_start_pull = upstream_75bp_1 - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    upstream_75bp_1 = exon_1[1] + 75
                    exon_1_start_pull = upstream_75bp_1 + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = exon_3[1] - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = exon_3[1] + 1
            #if both exons are less than 75bp long
            elif exon_1_size <= 75 and exon_3_size <= 75:
                #exon 1 end
                #+ strand
                if exon_1[0] < exon_1[1]:
                    exon_1_start_pull = exon_1[0] - 1
                    exon_1_end_pull = exon_1[1] - 1
                #- strand
                elif exon_1[0] > exon_1[1]:
                    exon_1_start_pull = exon_1[0] + 1
                    exon_1_end_pull = exon_1[1] + 1
                #exon 3 begining
                #+ strand
                if exon_3[0] < exon_3[1]:
                    exon_3_start_pull = exon_3[0] - 1
                    exon_3_end_pull = exon_3[1] - 1
                #- strand
                elif exon_3[0] > exon_3[1]:
                    exon_3_start_pull = exon_3[0] + 1
                    exon_3_end_pull = exon_3[1] + 1
            #pulling sequences for exon junction 1-3 - alternative splicing event
            #+ strand
            if exon_1_start_pull < exon_3_end_pull:
                as1_exon_1_seq = sequences[exon_1_start_pull:exon_1_end_pull]
                as1_exon_3_seq = sequences[exon_3_start_pull:exon_3_end_pull]
            elif exon_1_start_pull > exon_3_end_pull:
                as1_exon_1_seq = sequences[exon_1_end_pull:exon_1_start_pull]
                as1_exon_3_seq = sequences[exon_3_end_pull:exon_2_start_pull]
            as1_sequence = as1_exon_1_seq + as1_exon_3_seq
            header_values = [gene, chr_num, str(exon_trio_num),"as1"]
            as1_header = ".".join(header_values)
            dict_value = [as1_sequence, exon_1_size, exon_3_size]
            exon_trio_seqs_dict.update({as1_header:dict_value})
            exon_trio_num += 1
    return exon_trio_seqs_dict


#write exon trio fasta
#uses first value in each dictionary key
def write_exon_trio_fasta():
    sequences = pull_seqs()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for seq in sequences:
            single_value = sequences[seq]
            full_seq = single_value[0]
            final_seq = "".join(full_seq)
            header = ">" + seq + "\n"
            out.write(header)
            final = final_seq + "\n"
            out.write(final)


#write exon length for records
#uses 2nd and 3rd values from each dictionary key
def write_exon_trio_lengths():
    sequences = pull_seqs()
    all_lengths_dict = {}
    output = sys.argv[6]
    with open(output, 'a') as out:
        file_header = "Identifier\tExon.1.length\tExon.2.Length\tExon.3.Length\n"
        out.write(file_header)
        for seq in sequences:
            single_sequence = sequences[seq]
            header_split = seq.split(".")
            new_header_values = [header_split[0], header_split[1], header_split[2]]
            final_header = ".".join(new_header_values)
            first_exon_length = single_sequence[1]
            second_exon_length = single_sequence[2]
            dict_value = [first_exon_length, second_exon_length]
            if final_header in all_lengths_dict:
                all_lengths_dict[final_header].append(dict_value)
            elif final_header not in all_lengths_dict:
                all_lengths_dict.update({final_header:[dict_value]})
        for header in all_lengths_dict:
            single_value = all_lengths_dict[header]
            #single value format = [[exon 1 length, exon 2 length], [exon 2 length, exon 3 length], [exon 1 length, exon 3 length]]
            exon_1_length = single_value[0][0]
            exon_2_length = single_value[0][1]
            exon_3_length = single_value[1][1]
            final_value = "%s\t%s\t%s\t%s\n" % (str(header), str(exon_1_length), str(exon_2_length), str(exon_3_length))
            out.write(final_value)



#call all functions
def call():
    fasta_output = write_exon_trio_fasta()
    trios_lengths = write_exon_trio_lengths()

call()

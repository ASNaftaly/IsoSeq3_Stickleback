#Need to pull kmers from junction regions
#should be 51 kmers per junction where each kmer is 150bp long
#this script will be very simila rto PullSeqs_ExonUsage.py
#to run script: python3 PullSeqs_Junction.kmers.py <chromosome or scaffold number> <Ensembl_adjusted fasta file> <gene IDs with chromosome numbers file> <exon trios file> <output file, kmers fasta sequences>
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
#returns dictionary with key = gene.id_chr.num_exon.trio.num_scenario and value == [exon junction start, exon junction end]
def pull_junctions():
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
            #each sequence needs a unique identifier: gene.id_chr.num_exon.trio.num_scenario
            #there will be 2 scenarios (no alternative splicing and with alternative splicing)
            #no alternative splicing will have 2 sets of sequences (exon 1-2 and exon 2-3)
            #the labels for these scenarios will be: no alternative splicing = no_as1 (exon 1-2) and no_as2 (exon 2-3) and alternative splicing will be as1
            #first scenario is that all 3 exons are used (no alternative splicing)
            #this means I need to pull 100bp from the end of the first exon, 100bp from the beginning of exon 2 as one exon pair
            #then need to pull 100bp from the end of the second exon and 100bp from the beginning of exon 3
            #if the exons are not 100bp long, will pull the entire exon plus the upstream/downstream sequence to get to 100bp
            #will write a separate file with the lengths of each exon trio
            #pulling exon 1 sequences = only need the 100bp at the end of the exon
            if exon_1[0] < exon_1[1]:
                    upstream_100bp_1 = exon_1[1] - 100
                    exon_1_start_pull = upstream_100bp_1 - 1
                    exon_1_end_pull = exon_1[1] - 1
            #- strand
            elif exon_1[0] > exon_1[1]:
                upstream_100bp_1 = exon_1[1] + 100
                exon_1_start_pull = upstream_100bp_1 + 1
                exon_1_end_pull = exon_1[1] + 1
            #pulling exon 2 sequences = need both beginning and end of the exon
            #+ strand
            if exon_2[0] < exon_2[1]:
                downstream_100bp_2 = exon_2[0] + 100
                start_exon_2_start_pull = exon_2[0] - 1
                start_exon_2_end_pull = downstream_100bp_2 - 1
                upstream_100bp_2 = exon_2[1] - 100
                end_exon_2_start_pull = upstream_100bp_2 - 1
                end_exon_2_end_pull = exon_2[1] - 1
            #- strand
            elif exon_2[0] > exon_2[1]:
                downstream_100bp_2 = exon_2[0] - 100
                start_exon_2_start_pull = exon_2[0] + 1
                start_exon_2_end_pull = downstream_100bp_2 + 1
                upstream_100bp_2 = exon_2[1] + 100
                end_exon_2_start_pull = upstream_100bp_2 + 1
                end_exon_2_end_pull = exon_2[1] + 1
            #pulling exon 3 sequences, just need the beginning of the exon
            #+ strand
            if exon_3[0] < exon_3[1]:
                downstream_100bp_3 = exon_3[0] + 100
                exon_3_start_pull = exon_3[0] - 1
                exon_3_end_pull = downstream_100bp_3 - 1
            #- strand
            elif exon_3[0] > exon_3[1]:
                downstream_100bp_3 = exon_3[0] - 100
                exon_3_start_pull = exon_3[0] + 1
                exon_3_end_pull = downstream_100bp_3 + 1
            #pulling sequences for exon junction 1-2
            #+ strand
            if exon_1_start_pull < exon_1_end_pull:
                no_as1_exon_1_seq = sequences[exon_1_start_pull:exon_1_end_pull]
                no_as1_exon_2_seq = sequences[start_exon_2_start_pull:start_exon_2_end_pull]
            elif exon_1_start_pull > exon_1_end_pull:
                no_as1_exon_1_seq = sequences[exon_1_end_pull:exon_1_start_pull]
                no_as1_exon_2_seq = sequences[start_exon_2_end_pull:start_exon_2_start_pull]
            no_as1_sequence = no_as1_exon_1_seq + no_as1_exon_2_seq
            header_values = [gene, chr_num, str(exon_trio_num),"no_as1"]
            no_as1_header = ".".join(header_values)
            exon_trio_seqs_dict.update({no_as1_header:no_as1_sequence})
            #pulling sequences for exon junction 2-3
            #+ strand
            if end_exon_2_start_pull < end_exon_2_end_pull:
                no_as2_exon_2_seq = sequences[end_exon_2_start_pull:end_exon_2_end_pull]
                no_as2_exon_3_seq = sequences[exon_3_start_pull:exon_3_end_pull]
            elif end_exon_2_start_pull > end_exon_2_end_pull:
                no_as2_exon_2_seq = sequences[end_exon_2_end_pull:end_exon_2_start_pull]
                no_as2_exon_3_seq = sequences[exon_3_end_pull:exon_3_start_pull]
            no_as2_sequence = no_as2_exon_2_seq + no_as2_exon_3_seq
            header_values = [gene, chr_num, str(exon_trio_num),"no_as2"]
            no_as2_header = ".".join(header_values)
            exon_trio_seqs_dict.update({no_as2_header:no_as2_sequence})
            #pulling sequences for exon junction 1-3 - alternative splicing event
            #+ strand
            if exon_1_start_pull < exon_3_end_pull:
                as1_exon_1_seq = sequences[exon_1_start_pull:exon_1_end_pull]
                as1_exon_3_seq = sequences[exon_3_start_pull:exon_3_end_pull]
            elif exon_1_start_pull > exon_3_end_pull:
                as1_exon_1_seq = sequences[exon_1_end_pull:exon_1_start_pull]
                as1_exon_3_seq = sequences[exon_3_end_pull:exon_3_start_pull]
            as1_sequence = as1_exon_1_seq + as1_exon_3_seq
            header_values = [gene, chr_num, str(exon_trio_num),"as1"]
            as1_header = ".".join(header_values)
            exon_trio_seqs_dict.update({as1_header:as1_sequence})
            exon_trio_num += 1
    return exon_trio_seqs_dict

#pull junction kmers
#returns dictionary with key == junction header from previous function and value == kmers for each junction
def pull_kmers_seqs():
    junctions = pull_junctions()
    kmers_dict = {}
    for junction in junctions:
        single_junction_seq = junctions[junction]
        num_kmers = 51
        x = 0
        while x < num_kmers:
            single_kmer_seq = single_junction_seq[x:150+x]
            if junction in kmers_dict:
                kmers_dict[junction].append(single_kmer_seq)
            elif junction not in kmers_dict:
                kmers_dict.update({junction:[single_kmer_seq]})
            x += 1
    return kmers_dict


#write kmers to file
def write_kmers():
    kmers = pull_kmers_seqs()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for junction in kmers:
            kmer_number = 1
            single_junction = kmers[junction]
            for k in single_junction:
                header = ">%s.kmer%s\n" % (str(junction),str(kmer_number))
                single_seq = "".join(k)
                out.write(header)
                out.write(single_seq + "\n")
                kmer_number += 1


write_kmers()

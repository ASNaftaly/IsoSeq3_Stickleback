#script to pull specific region of sequence from fasta file
#mainly to check another script
#will need to read in full fasta file, split by chromosome then read in as standard input (chr number, strand, start position, end position)
#to run script: PullSeq_Check.py <chromosome fasta file> <specify isoform id> <specify chr num> <specify strand> <specify start pos>

import sys

def read_ensembl_fasta():
    fasta_file = sys.argv[1]
    fasta_dict = {}
    final_fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                new_line = line.split(" ")
                full_isoform_id = new_line[0].strip(" ")
                fasta_id = full_isoform_id.strip(">")
                final_fasta_id = fasta_id.strip("\n")
            else:
                new_line = line.strip("\n")
                if final_fasta_id in fasta_dict:
                    fasta_dict[final_fasta_id].append(new_line)
                elif final_fasta_id not in fasta_dict:
                    fasta_dict.update({final_fasta_id:[new_line]})
        for chr in fasta_dict:
            final_seq = []
            single_seq = fasta_dict[chr]
            for seq in single_seq:
                final_seq += seq
            final_fasta_dict.update({chr:final_seq})
    return final_fasta_dict

#pull specific sequence
def pull_seq():
    chromosome_seqs = read_ensembl_fasta()
    isoform = sys.argv[2]
    chr_num = sys.argv[3]
    strand = sys.argv[4]
    start_pos = int(sys.argv[5])
    output_file = sys.argv[6]
    with open(output_file, 'a') as out:
        single_chr = chromosome_seqs[chr_num]
        if strand == "+":
            upstream_pos = start_pos - 40
            header = ">%s\n" % str(isoform)
            pull_seq = single_chr[upstream_pos:start_pos+1]
            seq = "".join(pull_seq)
            final_seq = seq + "\n"
            out.write(header)
            out.write(final_seq)
        elif strand == "-":
            upstream_pos = start_pos + 40
            header = ">%s\n" % str(isoform)
            pull_seq = single_chr[start_pos:upstream_pos+1]
            seq = "".join(pull_seq)
            final_seq = seq + "\n"
            out.write(header)
            out.write(final_seq)


pull_seq()

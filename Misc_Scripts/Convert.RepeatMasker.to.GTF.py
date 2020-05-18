#converting repeat masker output into gtf format
#gtf format: chr.num \t source (repeat_masker) \t feature (type of repeat) \t start (start position) \t end (end position) \t score (.) \t strand (+/-) \t frame (0,1,2) \t attribute (semicolon list of tag-value pairs with additional information)
#repeat masker format: SW score  % substitutions  in matching region  % of bases gap  % of bases inserted  name of query sequence  start position query  end position query  #of bases in query past ending of match match is with the complement of the consensus sequence in database (i.e. strand) name of repeat repeat type
#to run script: python3 Convert.RepeatMasker.to.GTF.py <repeat masker output> <new repeat masker output in gtf format>
#Author: Alice Naftaly, May 2020

import sys


#read in repeat masker file
#nead columns = chr num, query start, query end, strand, name of repeat (for additional info), repeat type (for additional info and feature)
#returns dictionary with key == chr_num and value == [chr_num, query_start, query_end, strand, repeat_name, repeat_type]
def read_rm_file():
    rm_file = sys.argv[1]
    rm_dict = {}
    with open(rm_file, 'r') as repeat_masker:
        for line in repeat_masker:
            new_line = line.split()
            if len(new_line) == 13 or len(new_line) == 14 or len(new_line) == 0:
                continue
            else:
                chr_num = new_line[4]
                query_start = new_line[5]
                query_end = new_line[6]
                strand = new_line[8]
                repeat_name = new_line[9]
                repeat_type = new_line[10]
                dict_value = [chr_num, query_start, query_end, strand, repeat_name, repeat_type]
                if chr_num in rm_dict:
                    rm_dict[chr_num].append(dict_value)
                elif chr_num not in rm_dict:
                    rm_dict.update({chr_num:[dict_value]})
    return rm_dict


#create gtf lines
#return list of gtf lines (with data in proper order)
def create_gtf():
    rm_dict = read_rm_file()
    all_gtf_lines = []
    for chr in rm_dict:
        single_chr = rm_dict[chr]
        for entry in single_chr:
            chr_num = entry[0]
            start = entry[1]
            end = entry[2]
            strand_unchanged = entry[3]
            if strand_unchanged == "+":
                final_strand = strand_unchanged
            elif strand_unchanged == "C":
                final_strand = "-"
            repeat_name = entry[4]
            repeat_type_full = entry[5].split("/")
            if len(repeat_type_full) == 1:
                feature = repeat_type_full[0]
                attributes = repeat_name
            elif len(repeat_type_full) > 1:
                feature = repeat_type_full[0]
                attributes = repeat_type_full[1] + ";" + repeat_name
            gtf_line = [chr_num, "Repeat_Masker", feature, start, end, ".", final_strand, "0", attributes]
            all_gtf_lines.append(gtf_line)
    return all_gtf_lines


#write gtf
def write_gtf():
    gtf_lines = create_gtf()
    output = sys.argv[2]
    with open(output ,'a') as out:
        for line in gtf_lines:
            add_tabs = "\t".join(line)
            final_line = add_tabs + "\n"
            out.write(final_line)

write_gtf()

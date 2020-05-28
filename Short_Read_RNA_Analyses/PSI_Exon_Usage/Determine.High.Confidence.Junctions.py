#determine junctions with high confidence:
#based on read coverage = max(min(#no_as1, #no_as2),#as1) >= 5 AND min(#no_as1,#no_as2)+#as1 >= 10
#based on AS events with imbalances in read counts between exon pairs = absolute value (log2(#C1A/#Ac2)) <= 1 OR max(#C1A, #AC2) < #C1C2
#a high confidence junction will have Junction \t Y \t N
#output format:
#Junction.ID \t Sufficient.Read.Coverage \t Imbalance.in.Read.Counts \n
#to run script: python3 Determine.High.Confidence.Junctions.py <junctions bed file> <normalized read counts for each junction> <output file>
#Author: Alice Naftaly, May 2020

import sys
import math

#read in junctions positions file
#returns list of junctions ids without no_as1, no_as2, as1
def read_junction_positions():
    junction_file = sys.argv[1]
    junction_list = []
    with open(junction_file, 'r') as junctions:
        for line in junctions:
            new_line = line.split()
            junction_id_full = new_line[3].split(".")
            junction_id_new = junction_id_full[0] + "." + junction_id_full[1]
            junction_list.append(junction_id_new)
    return junction_list


#read in normalized read counts for each junction
#retuns a dictionary with key == junction id (gene_strand.exonnumber) and value == junction type (no_as1, no_as2, as1), read counts
def read_normalized_read_counts():
    read_counts_file = sys.argv[2]
    read_counts_dict = {}
    with open(read_counts_file, 'r') as read_counts:
        for line in read_counts:
            new_line = line.split()
            junction_id_full = new_line[0].split(".")
            junction_id_final = junction_id_full[0] + "." + junction_id_full[2]
            junction_type = junction_id_full[3]
            read_counts = new_line[1]
            counts = [junction_type, read_counts]
            if junction_id_final in read_counts_dict:
                read_counts_dict[junction_id_final].append(counts)
            elif junction_id_final not in read_counts_dict:
                read_counts_dict.update({junction_id_final:[counts]})
    return read_counts_dict


#adjust counts for normalized read counts based on if values are no_as1, no_as2, or as1
#will add in zeros where needed
#returns dictionary where key == junction id and value == counts for [no_as1, as1, no_as2]
def adjust_read_counts():
    read_counts_dict = read_normalized_read_counts()
    final_counts = {}
    for junction in read_counts_dict:
        single_junction = read_counts_dict[junction]
        junction_types = []
        junction_read_counts = []
        for a in single_junction:
            junction_types.append(a[0])
            junction_read_counts.append(a[1])
        if "no_as1" in junction_types:
            no_as1_index = junction_types.index("no_as1")
            no_as1_counts = str(junction_read_counts[no_as1_index])
        elif "no_as1" not in junction_types:
            no_as1_counts = "0"
        if "as1" in junction_types:
            as1_index = junction_types.index("as1")
            as1_counts = str(junction_read_counts[as1_index])
        elif "as1" not in junction_types:
            as1_counts = "0"
        if "no_as2" in junction_types:
            no_as2_index = junction_types.index("no_as2")
            no_as2_counts = str(junction_read_counts[no_as2_index])
        elif "no_as2" not in junction_types:
            no_as2_counts = "0"
        final_read_counts = [no_as1_counts, as1_counts, no_as2_counts]
        final_counts.update({junction:final_read_counts})
    return final_counts


#determine high enough read coverage to be sure of junction PSI calculations
#sufficient read coverage: Y == high confidence junction; N == low confidence junction
def sufficient_read_coverage():
    read_counts = adjust_read_counts()
    sufficient_read_cov_dict = {}
    for junction in read_counts:
        single_junction = read_counts[junction]
        no_as1 = float(single_junction[0])
        as1 = float(single_junction[1])
        no_as2 = float(single_junction[2])
        minimum_no_as = min(no_as1, no_as2)
        maximum_no_as_and_as = max(as1, minimum_no_as)
        addition_check = minimum_no_as + as1
        if maximum_no_as_and_as >= 5 and addition_check >= 10:
            coverage = "Y"
        else:
            coverage = "N"
        sufficient_read_cov_dict.update({junction:coverage})
    return sufficient_read_cov_dict


#determine if AS events are created by imbalance of read counts
#if imbalance == "Y" this means the reads are imbalanced across the junction; if imbalance == "N", the reads are good
def check_imbalance_read_counts():
    read_counts = adjust_read_counts()
    imbalance_read_counts_dict = {}
    for junction in read_counts:
        single_junction = read_counts[junction]
        no_as1 = float(single_junction[0])
        as1 = float(single_junction[1])
        no_as2 = float(single_junction[2])
        max_no_as_vs_as = max(no_as1, no_as2)
        if no_as2 == 0 and no_as1 == 0:
            imbalance_read_counts_dict.update({junction:"N"})
        elif no_as2 == 0 and no_as1 != 0:
            if max_no_as_vs_as < as1:
                imbalance = "N"
            else:
                imbalance = "Y"
        elif no_as2 != 0 and no_as1 != 0:
            log2_no_as = abs(math.log2(no_as1/no_as2))
            if max_no_as_vs_as < as1 or log2_no_as <= 1:
                imbalance = "N"
            else:
                imbalance = "Y"
        imbalance_read_counts_dict.update({junction:imbalance})
    return imbalance_read_counts_dict



#write output
def write():
    read_coverage = sufficient_read_coverage()
    imbalances = check_imbalance_read_counts()
    output = sys.argv[3]
    with open(output ,'a') as out:
        header = "Junction.ID\tSufficient.Read.Coverage\tImbalance.in.Read.Counts\n"
        out.write(header)
        for junction in read_coverage:
            single_read_cov = read_coverage[junction]
            single_imbalance = imbalances[junction]
            final = "%s\t%s\t%s\n" % (str(junction), str(single_read_cov), str(single_imbalance))
            out.write(final)


write()

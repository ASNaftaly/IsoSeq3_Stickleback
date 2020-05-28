#calculating PSI for exon junctions
#will calculate separately for each exon junction
#will read in normalized read counts for each junction and  junction positions
#junctions position file can be any subset of junctions as long as it is in bed format
#to run script: python3 Calc.PSI.per.Junction.py < junction positions bed file> <normalized read counts for junctions> <output file for PSI calcs>
#Author: Alice Naftaly, May 2020

import sys

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


#calculating PSI (before *100)
def calc_psi():
    read_counts = adjust_read_counts()
    psi_values = {}
    for junction in read_counts:
        single_junction = read_counts[junction]
        no_as1 = float(single_junction[0])
        as1 = float(single_junction[1])
        no_as2 = float(single_junction[2])
        average_no_alternative_splicing = (no_as1 + no_as2) / 2
        psi_calc = average_no_alternative_splicing / (as1 + average_no_alternative_splicing)
        percentage_psi_calc = round(100 * psi_calc,2)
        psi_values.update({junction:percentage_psi_calc})
    return psi_values


#write PSI values to a file
def write():
    psi_values = calc_psi()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Junction.ID\tPSI.Percentage\n"
        for junction in psi_values:
            final = "%s\t%s\n" % (str(junction), str(psi_values[junction]))
            out.write(final)


write()


calc_psi()

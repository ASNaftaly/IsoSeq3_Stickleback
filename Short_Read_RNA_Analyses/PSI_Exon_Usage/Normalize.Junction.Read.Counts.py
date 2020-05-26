#normalize junction read counts
#need to read in normalization factor for each junction and then read counts for each junction
#there will be junctions with no reads
#to run script: Normalize.Junction.Read.Counts.py <normalization factor file> <raw read counts for each junction> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read in normalization factors
#returns dictionary with key == junction_id and value == normalization factor
def read_normalization_factors():
    norm_factor_file = sys.argv[1]
    norm_factor_dict = {}
    with open(norm_factor_file, 'r') as nf_info:
        for line in nf_info:
            new_line = line.split()
            junction_id = new_line[0]
            norm_factor = float(new_line[1])
            norm_factor_dict.update({junction_id:norm_factor})
    return norm_factor_dict

#read raw read counts file
#returns dictionary with key == junction id and value == read counts at junction
def read_raw_read_counts():
    read_counts_file = sys.argv[2]
    read_counts_dict = {}
    with open(read_counts_file, 'r') as read_counts:
        for line in read_counts:
            new_line = line.split()
            junction_id = new_line[0]
            junction_counts = int(new_line[1])
            read_counts_dict.update({junction_id:junction_counts})
    return read_counts_dict

#normalizaing read counts by norm factor
#returns dictionary with key == junction id and value == normalized read count
def normalize():
    norm_factors = read_normalization_factors()
    raw_read_counts = read_raw_read_counts()
    normalized_read_counts = {}
    for junction in raw_read_counts:
        if junction in norm_factors:
            single_raw_read_count = raw_read_counts[junction]
            single_norm_factor = norm_factors[junction]
            normalized_read_count = round(single_raw_read_count * single_norm_factor, 6q)
            normalized_read_counts.update({junction:normalized_read_count})
    return normalized_read_counts

#write normalized read counts to new file
def write():
    normalized_read_counts = normalize()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for junction in normalized_read_counts:
            single_junction = normalized_read_counts[junction]
            final = "%s\t%s\n" % (str(junction), str(single_junction))
            out.write(final)

write()

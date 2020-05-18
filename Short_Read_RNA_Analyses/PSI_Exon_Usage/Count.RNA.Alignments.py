#This is very similar to Count.Kmers.Alignments.py, but for RNA reads, the unique alignment is for the read, so each read should only align once
#Will also put all of the duplicate reads into a separate file to examine if need be
#reads in bed file created from mapq 30 filtered and sorted bam file and outputs a file with the following format:
#Junction.Name \t #.Reads
#also need to filter out read that do not overlap the boundary (92-108bp) as 100bp from each exon was pulled
#to run script: python3 Count.RNA.Alignments.py <bed file> <output file; unique reads that span exon boundary> <output file; duplicate reads per junction> <output file; reads that do not span the boundary>
#Author: Alice Naftaly, May 2020


import sys

#read in bed file
#returns dictionary with key == junction identifier (chr num + strand with gene id and exon number as well as no_as1 (exon 1-2 junction), no_as2 (exon 2-3 junction), or as1 (exon 1-3 junction))
def read_bed():
    bed_file = sys.argv[1]
    bed_dict = {}
    with open(bed_file, 'r') as bed_info:
        for line in bed_info:
            new_line = line.split()
            junction_iden = new_line[0]
            start_pos = int(new_line[1])
            end_pos = int(new_line[2])
            read_iden = new_line[3]
            dict_value = [start_pos, end_pos, read_iden]
            if junction_iden in bed_dict:
                bed_dict[junction_iden].append(dict_value)
            elif junction_iden not in bed_dict:
                bed_dict.update({junction_iden:[dict_value]})
    return bed_dict


#filter out reads that don't cover the boundary
#the boundary is 92-108bp for each junction
#returns 2 lists, 1 with reads that do cover the boundary and 1 with reads that do not fully cover the boundary
def filter_cover_boundary():
    bed_dict = read_bed()
    reads_do_not_cover_boundary = []
    reads_cover_boundary = []
    for junction in bed_dict:
        single_junction = bed_dict[junction]
        read_starts = []
        read_ends = []
        read_names = []
        for read in single_junction:
            read_starts.append(int(read[0]))
            read_ends.append(int(read[1]))
            read_names.append(read[2])
        for ind, value in enumerate(read_starts):
            #if reads span full boundary
            if value <= 92 and read_ends[ind] >= 108:
                reads_cover_boundary.append(read_names[ind])
            #if reads do not span full boundary (includes reads that partially cover boundary)
            else:
                reads_do_not_cover_boundary.append(read_names[ind])
    return reads_cover_boundary, reads_do_not_cover_boundary


#remove any potential duplicate read mapping
#returns 3 lists = 1 = final unique reads across the boundary, duplicate reads across the boundary, reads that do not span the boundary
def remove_read_duplicates():
    reads_across_boundary, reads_not_across_boundary = filter_cover_boundary()
    #identify which reads are duplicated
    all_reads = set()
    duplicate_reads = []
    for idx, val in enumerate(reads_across_boundary):
        if val not in all_reads:
            all_reads.add(val)
        else:
            duplicate_reads.append(val)
    duplicate_read_alignments = []
    final_unique_reads = []
    #get indeces for duplicate values
    for overlap in duplicate_reads:
        duplicate_read_alignments += [i for i, x in enumerate(reads_across_boundary) if x == overlap]
    #remove indeces from reads_across_boundary
    for index, k in enumerate(reads_across_boundary):
        if index in duplicate_read_alignments:
            continue
        elif index not in duplicate_read_alignments:
            final_unique_reads.append(k)
    return final_unique_reads, duplicate_reads, reads_not_across_boundary



#count reads across boundary, duplicate reads, and reads not across boundary

#counts non duplicate reads that span boundary
#will only add reads that fit criteria
#returns dictionary with key == junction identifier and value == number of reads that fit criteria
def count_final_reads():
    print("Filtered out duplicate reads")
    final_reads, duplicates, reads_not_spanning_boundary = remove_read_duplicates()
    bed_dict = read_bed()
    final_junction_read_count = {}
    for junction in bed_dict:
        single_junction = bed_dict[junction]
        aligned_junction_read_count = 0
        for read in single_junction:
            read_name = read[2]
            if read_name in final_reads:
                aligned_junction_read_count += 1
        if aligned_junction_read_count > 0:
            final_junction_read_count.update({junction:str(aligned_junction_read_count)})
    print("Counted boundary spanning unique reads")
    return final_junction_read_count


#counts duplicate reads
#will return values that have 0 duplicates as well as those that have more
#returns dictionary with key == junction identifier and value == number of duplicate reads
def count_duplicate_reads():
    final_reads, duplicates, reads_not_spanning_boundary = remove_read_duplicates()
    bed_dict = read_bed()
    duplicate_read_count = {}
    for junction in bed_dict:
        single_junction = bed_dict[junction]
        aligned_junction_read_count = 0
        for read in single_junction:
            read_name = read[2]
            if read_name in duplicates:
                aligned_junction_read_count += 1
        duplicate_read_count.update({junction:str(aligned_junction_read_count)})
    print("Counted duplicate reads")
    return duplicate_read_count


#count reads that do not span junction
#this can also produce 0 reads for a given junction, if no reads do not span the exon pair boundary
def count_nonspanning_boundary_reads():
    final_reads, duplicates, reads_not_spanning_boundary = remove_read_duplicates()
    bed_dict = read_bed()
    nonspanning_read_count = {}
    for junction in bed_dict:
        single_junction = bed_dict[junction]
        aligned_junction_read_count = 0
        for read in single_junction:
            read_name = read[2]
            if read_name in reads_not_spanning_boundary:
                aligned_junction_read_count += 1
        nonspanning_read_count.update({junction:str(aligned_junction_read_count)})
    print("Counted non boundary spanning reads")
    return nonspanning_read_count



#write all files

def write_nonduplicate_boundaryspanning_reads():
    reads_to_count = count_final_reads()
    output = sys.argv[2]
    with open(output, 'a') as out:
        header = "Junction.Identifier\tNo.of.Reads\n"
        out.write(header)
        for junction in reads_to_count:
            single_junction = reads_to_count[junction]
            final = "%s\t%s\n" % (str(junction), str(single_junction))
            out.write(final)


def write_duplicate_boundaryspanning_reads():
    reads_to_count = count_duplicate_reads()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Junction.Identifier\tNo.of.Reads\n"
        out.write(header)
        for junction in reads_to_count:
            single_junction = reads_to_count[junction]
            final = "%s\t%s\n" % (str(junction), str(single_junction))
            out.write(final)


def write_non_boundaryspanning_reads():
    reads_to_count = count_nonspanning_boundary_reads()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Junction.Identifier\tNo.of.Reads\n"
        out.write(header)
        for junction in reads_to_count:
            single_junction = reads_to_count[junction]
            final = "%s\t%s\n" % (str(junction), str(single_junction))
            out.write(final)


#call all functions
def call():
    PSI_counts = write_nonduplicate_boundaryspanning_reads()
    duplicates = write_duplicate_boundaryspanning_reads()
    non_boundaryspanning = write_non_boundaryspanning_reads()

call()

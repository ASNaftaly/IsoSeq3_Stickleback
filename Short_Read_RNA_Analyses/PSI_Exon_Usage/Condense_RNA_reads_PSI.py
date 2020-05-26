#condense RNA unique reads across junctions
#had to split file into 14 smaller files to run script; now I need to combine all of the junctions that have counts in different files
#used cat to combine 14 files and now will read them into a dictionary and condense them
#to run script: python3 Condense_RNA_reads_PSI.py <uncondensed file with all RNA counts per junction> <output file>
#Author: Alice Naftaly, May 2020

import sys

#read in uncondensed file
#returns dictionary with key == junction id and value == junction counts
def read_uncondensed():
    input_file = sys.argv[1]
    junction_dict = {}
    with open(input_file, 'r') as info:
        for line in info:
            if line.startswith("Junction.Identifier"):
                continue
            else:
                new_line = line.split()
                junction_id = new_line[0]
                junction_counts = new_line[1]
                if junction_id in junction_dict:
                    junction_dict[junction_id].append(junction_counts)
                elif junction_id not in junction_dict:
                    junction_dict.update({junction_id:[junction_counts]})
    return junction_dict


#collapse any junction that has more than one set of counts
#returns dictionary where every key(junction) has one total count
def collapse_counts():
    junction_dict = read_uncondensed()
    final_junction_dict = {}
    for junction in junction_dict:
        single_junction = junction_dict[junction]
        if len(single_junction) == 1:
            final_junction_dict.update({junction:single_junction[0]})
        elif len(single_junction) > 1:
            final_counts = int(single_junction[0]) + int(single_junction[1])
            final_junction_dict.update({junction:str(final_counts)})
    return final_junction_dict

#write new junction counts:
def write():
    junction_dict = collapse_counts()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for junction in junction_dict:
            single_junction = junction_dict[junction]
            final = "%s\t%s\n" % (str(junction), str(single_junction))
            out.write(final)

write()

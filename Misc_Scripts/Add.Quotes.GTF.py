#add "" to attributes tags in gtf so other programs will read them correctly
#have to use "" not '' as the single quotes don't work properly
#to run script: python3 Add.Quotes.GTF.py <final gtf> <output formated gtf>
#Author: Alice Naftaly, Sept 2020

import sys

#read in gtf file
def read_gtf():
    gtf_file = sys.argv[1]
    output = sys.argv[2]
    with open(gtf_file, 'r') as gtf, open(output, 'a') as out:
        for line in gtf:
            new_line = line.split("\t")
            attributes = new_line[8].split(";")
            formatted_attributes = []
            for att_tag in attributes:
                split_tag = att_tag.split(" ")
                if len(split_tag) == 2:
                    gene_id = str(split_tag[0]) + " " + '"' + str(split_tag[1]) + '"'
                    formatted_attributes.append(gene_id)
                elif len(split_tag) == 3:
                    other_atts = " " + str(split_tag[1]) + " " + '"' + str(split_tag[2]) + '"'
                    formatted_attributes.append(other_atts)
            joined_attributes = ";".join(formatted_attributes)
            final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0], new_line[1], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], joined_attributes)
            out.write(final)
read_gtf()

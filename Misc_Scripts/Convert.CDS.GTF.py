#CDS regions for - strand isoforms is backwards (need to have smaller number as start)
#to run script: python3 Convert.CDS.GTF.py <input gtf file> <output gtf file>
#Author: Alice Naftaly, Sept 2020

import sys

#read in gtf file
def read_gtf():
    gtf_file = sys.argv[1]
    output = sys.argv[2]
    with open(gtf_file, 'r') as gtf, open(output, 'a') as out:
        for line in gtf:
            new_line = line.split("\t")
            if new_line[4] == ".":
                final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0], new_line[1], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], new_line[8])
                #out.write(final)
            elif new_line[4] != ".":
                start = int(new_line[3])
                stop = int(new_line[4])
                if start < stop:
                    final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0], new_line[1], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], new_line[8])
                    out.write(final)
                elif stop < start:
                    final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (new_line[0], new_line[1], new_line[2], str(stop), str(start), new_line[5], new_line[6], new_line[7], new_line[8])
                    out.write(final)


read_gtf()

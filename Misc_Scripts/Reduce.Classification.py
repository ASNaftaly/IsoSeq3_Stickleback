#need to remove first column from converted classification files and remove extra empty lines
#this is to run other analyses such as using short read RNA-seq to confirm isoforms for each tissue and sex
#will keep the file with both the combined sexes isoform and the individual tissue isoforms separately
#to run script: python3 Reduce.Classification.py <converted classification file> <final classification file output>


import sys

#read in classification file
#returns dictionary with key == "header" or "isoform" and value == full line from file in list format (separated by tabs)
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("combined.sexes.isoform"):
                new_line = line.split()
                header = new_line[1:len(new_line)]
                class_dict.update({"header":header})
            elif line.startswith("PB"):
                new_line = line.split()
                reduced_line = new_line[1:len(new_line)]
                if "isoform" in class_dict:
                    class_dict["isoform"].append(reduced_line)
                elif "isoform" not in class_dict:
                    class_dict.update({"isoform":[reduced_line]})
    return class_dict

#write new output
def write():
    class_dict = read_class()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for key in class_dict:
            if key == "header":
                header_line = "\t".join(class_dict[key])
                final = header_line + "\n"
                out.write(final)
            elif key == "isoform":
                single_key = class_dict[key]
                for isoform in single_key:
                    new_line = "\t".join(isoform)
                    final = new_line + "\n"
                    out.write(final)

write()

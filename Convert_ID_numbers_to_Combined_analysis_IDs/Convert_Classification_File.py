#creating new classification file that has an added column with the combined Sexes isoform id
#need to read in original filtered classification file and then add a column before the other columns for the combined sexes isoform ID
#will need to read in tab delimited file with isoforms from combined sexes and single tissue analysis
#to run script: python3 Convert_Classificatoin_File.py <isoform comparison file from Sort_best_matches_BLAST.py> <individual sexes classification file> <output, new classification file with 1 extra column at beginning> <output, excluded isoforms from single tissue analysis>
#Author: Alice Naftaly, February 2020

import sys

#read in isoform comparison file
#input format = Individual tissue isoform ID \t combined sexes isoform ID
#returns dictionary with key = individual tissue isoform id and value == combined tissues isoform id
def read_isoform_comp():
    isoform_comp_file = sys.argv[1]
    isoform_comp_dict = {}
    with open(isoform_comp_file, 'r') as isoform_comp:
        for line in isoform_comp:
            if line.startswith("PB"):
                new_line = line.split()
                individual_tissue_id = new_line[0]
                combined_tissues_id = new_line[1]
                isoform_comp_dict.update({individual_tissue_id:combined_tissues_id})
    return isoform_comp_dict

#read in individual tissue classification file
#will want to have header line as an entry as well as all of the isoforms
#returns dictionary with key = header or isoform id and value == full line
def read_class():
    class_file = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_input:
        for line in class_input:
            if line.startswith("isoform"):
                new_line = line.split()
                class_dict.update({"header":new_line})
            else:
                new_line = line.split("\t")
                isoform_id = new_line[0]
                class_dict.update({isoform_id:new_line})
    return class_dict

#create new dictionary combining new isoform id from combined sexes to classification data
#returns 2 dictionaries
#dictionary 1 = combined_dict = key = isoform, value = list of line values + combined sexes isoform
#dictionary 2 = excluded dict, same format as combined dict, but only includes single tissue isoforms that were not matched to combined sexes isoforms
def combine():
    isoform_comps = read_isoform_comp()
    class_dict = read_class()
    combined_dict = {}
    excluded_dict = {}
    for key in class_dict:
        if key == "header":
            header_line = class_dict[key]
            excluded_dict.update({"header":header_line})
            header_line.insert(0,"combined.sexes.isoform")
            combined_dict.update({"header":header_line})
        elif key in isoform_comps:
            combined_sexes_isoform = isoform_comps[key]
            single_class_line = class_dict[key]
            single_class_line.insert(0, combined_sexes_isoform)
            combined_dict.update({key:single_class_line})
        elif key not in isoform_comps:
            single_class_line = class_dict[key]
            excluded_dict.update({key:single_class_line})
    print("Number of Isoforms kept after filter 1")
    print(len(combined_dict))
    return combined_dict, excluded_dict

#write combined and excluded dictionaries to files
def write_combined():
    combined_dict, excluded_dict = combine()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for key in combined_dict:
            if key == "header":
                final = "\t".join(combined_dict[key])
                final_out = final + "\n"
                out.write(final_out)
            else:
                final = "\t".join(combined_dict[key])
                final_out = final + "\n"
                out.write(final_out)

def write_excluded():
    combined_dict, excluded_dict = combine()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for key in excluded_dict:
            if key == "header":
                final = "\t".join(excluded_dict[key])
                final_out = final + "\n"
                out.write(final_out)
            else:
                final = "\t".join(excluded_dict[key])
                final_out = final + "\n"
                out.write(final_out)


#calling write functions
def call():
    combined = write_combined()
    excluded = write_excluded()

call()

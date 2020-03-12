#Want to know how many isoforms when BLAST directly to same isoforms (database is created from fasta file used as BLAST input) a full BLAST matches vs partial matches

import sys

#read in blast output
#returns dictionary with key = isoform id and value = alignment
def read_blast():
    blast_output = sys.argv[1]
    blast_dict = {}
    with open(blast_output, 'r') as blast_results:
        for line in blast_results:
            if line.startswith("#"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                if isoform_id in blast_dict:
                    blast_dict[isoform_id].append(new_line)
                elif isoform_id not in blast_dict:
                    blast_dict.update({isoform_id:[new_line]})
    print("Read BLAST Output File")
    return blast_dict


#need the length of combined sexes isoforms
#will create a dictionary with key = isoform id and value = length of isoform
def pull_combined_isoform():
    combined_isoform_file = sys.argv[2]
    combined_isoform_dict = {}
    with open(combined_isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("isoform"):
                continue
            else:
                new_line = line.split()
                isoform_id = new_line[0]
                length = new_line[3]
                combined_isoform_dict.update({isoform_id:length})
    print("Read combined sexes classification file and created isoform length dictionary")
    return combined_isoform_dict

#compare isoform length with BLAST alignment lengths
def compare():
    combined_isoforms = pull_combined_isoform()
    blast_output = read_blast()
    same_isoforms = []
    for key in blast_output:
        single_blast = blast_output[key]
        single_isoform = combined_isoforms[key]
        for value in single_blast:
            alignment_length = value[3]
            matching_isoform = value[1]
            if alignment_length == single_isoform:
                same_isoforms.append(single_isoform)
                '''print(key)
                print(matching_isoform)
                print(alignment_length)
                print(single_isoform)'''
    print(len(combined_isoforms))
    print(len(set(same_isoforms)))
    print(len(same_isoforms))



compare()

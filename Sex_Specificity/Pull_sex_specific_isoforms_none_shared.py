#pull sex specific isoforms that have no shared isoforms
#to run script: python3 Pull_sex_specific_isoforms_none_shared.py <male specific isoforms by splicing file> <female specific isoforms by splicing file> <shared isoforms file> <male single isoforms output file> <female single isoforms output file>
#Author: Alice Naftaly Sept 2020

import sys

#read male specific isoforms
def read_male_isoforms():
    male_file = sys.argv[1]
    male_isoforms = {}
    with open(male_file,'r') as male_isos:
        for line in male_isos:
            male_isoform = line.strip("\n")
            male_iso_split = male_isoform.split(".")
            gene_id = male_iso_split[0] + "." + male_iso_split[1]
            if gene_id in male_isoforms:
                male_isoforms[gene_id].append(male_isoform)
            elif gene_id not in male_isoforms:
                male_isoforms.update({gene_id:[male_isoform]})
    return male_isoforms


#read female specific isoforms
def read_female_isoforms():
    female_file = sys.argv[2]
    female_isoforms = {}
    with open(female_file,'r') as female_isos:
        for line in female_isos:
            female_isoform = line.strip("\n")
            female_iso_split = female_isoform.split(".")
            gene_id = female_iso_split[0] + "." + female_iso_split[1]
            if gene_id in female_isoforms:
                female_isoforms[gene_id].append(female_isoform)
            elif gene_id not in female_isoforms:
                female_isoforms.update({gene_id:[female_isoform]})
    return female_isoforms

#shared isoforms:
#read male specific isoforms
def read_shared_isoforms():
    shared_file = sys.argv[3]
    shared_isoforms = {}
    with open(shared_file,'r') as shared_isos:
        for line in shared_isos:
            shared_isoform = line.strip("\n")
            shared_iso_split = shared_isoform.split(".")
            gene_id = shared_iso_split[0] + "." + shared_iso_split[1]
            if gene_id in shared_isoforms:
                shared_isoforms[gene_id].append(shared_isoform)
            elif gene_id not in shared_isoforms:
                shared_isoforms.update({gene_id:[shared_isoform]})
    return shared_isoforms


#pull isoform pairs where no isoforms are shared
def compare():
    shared_isoforms = read_shared_isoforms()
    male_isoforms = read_male_isoforms()
    female_isoforms = read_female_isoforms()
    final_dict = {}
    for gene in male_isoforms:
        if gene in female_isoforms and gene not in shared_isoforms:
            single_female_isoforms = female_isoforms[gene]
            single_male_isoforms = male_isoforms[gene]
            if len(single_female_isoforms) == 1 and len(single_male_isoforms) == 1:
                final_dict.update({gene:[male_isoforms[gene]]})
                final_dict[gene].append(female_isoforms[gene])
    return final_dict


def write_output():
    final_dict = compare()
    male_output = sys.argv[4]
    female_output = sys.argv[5]
    with open(male_output, 'a') as out, open(female_output, 'a') as out2:
        for gene in final_dict:
            male_isos = final_dict[gene][0]
            female_isos = final_dict[gene][1]
            for iso in male_isos:
                out.write(iso + "\n")
            for iso2 in female_isos:
                out2.write(iso2 + "\n")


write_output()

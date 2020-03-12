#comparing isoforms between single tissue analyses
#this will be easier that before as all of the single tissue analyses have the combined sexes isoform ID added to them, so I can just compare IDs rather than compare sequences, locations, etc.
#also for protein coding genes, there are several isoforms that are predicted to code for the same proteins, so want to take this into account
#to run script: python3 Compare_single_tissues.py <single tissue classification file female liver> <single tissue classification file male liver> <single tissue classification file female brain> <single tissue classification file male brain> <single tissue classification file female pronephros> <single tissue classification file male pronephros> <single tissue classification file ovary> <single tissue classification file testis> <collapsed protein blast results from Remove_duplicates_protein_BLAST.py>
#have 8 single tissue files for this analysis
#Author: Alice Naftaly, March 2020

import sys

#read in classification files for each single tissue
#return lists of isoform IDs
def read_female_liver_class():
    class_file = sys.argv[1]
    female_liver_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_liver_isoforms.append(isoform)
    return female_liver_isoforms

def read_male_liver_class():
    class_file = sys.argv[2]
    male_liver_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_liver_isoforms.append(isoform)
    return male_liver_isoforms


def read_female_brain_class():
    class_file = sys.argv[3]
    female_brain_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_brain_isoforms.append(isoform)
    return female_brain_isoforms

def read_male_brain_class():
    class_file = sys.argv[4]
    male_brain_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_brain_isoforms.append(isoform)
    return male_brain_isoforms


def read_female_pronephros_class():
    class_file = sys.argv[5]
    female_pronephros_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_pronephros_isoforms.append(isoform)
    return female_pronephros_isoforms

def read_male_pronephros_class():
    class_file = sys.argv[6]
    male_pronephros_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_pronephros_isoforms.append(isoform)
    return male_pronephros_isoforms

def read_ovary_class():
    class_file = sys.argv[7]
    ovary_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                ovary_isoforms.append(isoform)
    return ovary_isoforms

def read_testis_class():
    class_file = sys.argv[8]
    testis_isoforms = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                testis_isoforms.append(isoform)
    return testis_isoforms

#read in protein coding genes that are collapsed
#file created from Remove_duplicates_protein_BLAST.py
#pci = protein coding isoforms
def pull_collapsed_pci():
    input_file = sys.argv[9]
    count = 0
    collapsed_dict = {}
    with open(input_file, 'r') as collapsed_proteins:
        for line in collapsed_proteins:
            new_line = line.split()
            collapsed_dict.update({count:new_line})
            count += 1
    return collapsed_dict

#now to compare isoforms from different tissues or sets
#all tissues:
def compare_all_tissues():
    female_liver_isoforms = read_female_liver_class()
    male_liver_isoforms = read_male_liver_class()
    female_brain_isoforms = read_female_brain_class()
    male_brain_isoforms = read_male_brain_class()
    female_pronephros_isoforms = read_female_pronephros_class()
    male_pronephros_isoforms = read_male_pronephros_class()
    ovary_isoforms = read_ovary_class()
    testis_isoforms = read_testis_class()
    collapsed_pci = pull_collapsed_pci()
    all_tissues = 0
    all_somatic_tissues = 0
    all_male_tissues = 0
    all_female_tissues = 0
    liver_tissues = 0
    brain_tissues = 0
    pronephros_tissues = 0
    gonad_tissues = 0
    for key in collapsed_pci:
        single_key = collapsed_pci[key]
        if len(single_key) == 1:
            test_isoform = single_key[0]
            if test_isoform in female_liver_isoforms and test_isoform in male_liver_isoforms and test_isoform in female_brain_isoforms and test_isoform in male_brain_isoforms and test_isoform in female_pronephros_isoforms and test_isoform in male_pronephros_isoforms and test_isoform in ovary_isoforms and test_isoform in testis_isoforms:
                all_tissues += 1
            elif test_isoform in female_liver_isoforms and test_isoform in male_liver_isoforms and test_isoform in female_brain_isoforms and test_isoform in male_brain_isoforms and test_isoform in female_pronephros_isoforms and test_isoform in male_pronephros_isoforms:
                all_somatic_tissues += 1
            elif test_isoform in female_liver_isoforms and test_isoform in female_brain_isoforms and test_isoform in female_pronephros_isoforms and test_isoform in ovary_isoforms:
                all_female_tissues += 1
            elif test_isoform in male_liver_isoforms and test_isoform in male_brain_isoforms and test_isoform in male_pronephros_isoforms and test_isoform in testis_isoforms:
                all_male_tissues += 1
            elif test_isoform in female_liver_isoforms and test_isoform in male_liver_isoforms:
                liver_tissues += 1
            elif test_isoform in female_brain_isoforms and test_isoform in male_brain_isoforms:
                brain_tissues += 1
            elif test_isoform in female_pronephros_isoforms and test_isoform in male_pronephros_isoforms:
                pronephros_tissues += 1
            if test_isoform in ovary_isoforms and test_isoform in testis_isoforms:
                gonad_tissues += 1
        elif len(single_key) > 1:
            temp_list = []
            for test_isoform in single_key:
                if test_isoform in female_liver_isoforms:
                    new_value = "female_liver"
                    temp_list.append(new_value)
                if test_isoform in male_liver_isoforms:
                    new_value = "male_liver"
                    temp_list.append(new_value)
                if test_isoform in female_brain_isoforms:
                    new_value = "female_brain"
                    temp_list.append(new_value)
                if test_isoform in male_brain_isoforms:
                    new_value = "male_brain"
                    temp_list.append(new_value)
                if test_isoform in female_pronephros_isoforms:
                    new_value = "female_pronephros"
                    temp_list.append(new_value)
                if test_isoform in male_pronephros_isoforms:
                    new_value = "male_pronephros"
                    temp_list.append(new_value)
                if test_isoform in ovary_isoforms:
                    new_value = "ovary"
                    temp_list.append(new_value)
                if test_isoform in testis_isoforms:
                    new_value = "testis"
                    temp_list.append(new_value)
            set_temp_list = list(set(temp_list))
            if len(set_temp_list) == 8:
                all_tissues += 1
            elif len(set_temp_list) == 2:
                if "female_liver" in set_temp_list and "male_liver" in set_temp_list:
                    liver_tissues += 1
                elif "female_brain" in set_temp_list and "male_brain" in set_temp_list:
                    brain_tissues += 1
                elif "female_pronephros" in set_temp_list and "male_pronephros" in set_temp_list:
                    pronephros_tissues += 1
                elif "ovary" in set_temp_list and "testis" in set_temp_list:
                    gonad_tissues += 1
            elif len(set_temp_list) == 6:
                if "female_liver" in set_temp_list and "male_liver" in set_temp_list and "female_brain" in set_temp_list and "male_brain" in set_temp_list and "female_pronephros" in set_temp_list and "male_pronephros" in set_temp_list:
                    all_somatic_tissues += 1
            elif len(set_temp_list) == 4:
                if "female_liver" in set_temp_list and "female_brain" in set_temp_list and "female_pronephros" in set_temp_list and "ovary" in set_temp_list:
                    all_female_tissues += 1
                elif "male_liver" in set_temp_list and "male_brain" in set_temp_list and "male_pronephros" in set_temp_list and "testis" in set_temp_list:
                    all_male_tissues += 1
    print(all_tissues)
    print(all_somatic_tissues)
    print(all_female_tissues)
    print(all_male_tissues)
    print(liver_tissues)
    print(brain_tissues)
    print(pronephros_tissues)
    print(gonad_tissues)


compare_all_tissues()

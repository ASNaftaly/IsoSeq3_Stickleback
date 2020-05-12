#splitting male and female classification files into 3 classification files (1 with all shared isoforms, 1 with male specific isoforms and 1 with female specific isoforms)
#
#to run script: python3 Sort.by.Sex.Classification.py <male full classification file> <female full classification file> 

import sys


#read male classification file
#returns dictionary with key == isoform and value == classification line
def read_male_file():
    male_file = sys.argv[1]
    male_dict = {}
    with open(male_file, 'r') as male_class:
        for line in male_class:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_dict.update({isoform:line})
    return male_dict


#read female classification file
#returns dictionary with key == isoform and value == classification line
def read_female_file():
    female_file = sys.argv[2]
    female_dict = {}
    with open(female_file, 'r') as female_class:
        for line in female_class:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_dict.update({isoform:line})
    return female_dict

#split isoforms into 3 lists
def split_isoforms():
    male_dict = read_male_file()
    all_male_isoforms = list(male_dict.keys())
    female_dict = read_female_file()
    all_female_isoforms = list(female_dict.keys())
    print(all_male_isoforms)

split_isoforms()

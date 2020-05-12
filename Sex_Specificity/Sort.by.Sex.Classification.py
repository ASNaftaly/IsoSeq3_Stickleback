#splitting male and female classification files into 3 classification files (1 with all shared isoforms, 1 with male specific isoforms and 1 with female specific isoforms)
#to run script: python3 Sort.by.Sex.Classification.py <male full classification file> <female full classification file> <output all shared isoforms with 1 isoform per line> < output all male specific isoforms> <output all female specific isoforms>
#Author: Alice Naftaly, May 2020

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
#returns 3 lists of isoforms (shared, male specific, female specific)
def split_isoforms():
    male_dict = read_male_file()
    all_male_isoforms = list(male_dict.keys())
    female_dict = read_female_file()
    all_female_isoforms = list(female_dict.keys())
    shared_isoforms = []
    male_specific_isoforms = []
    female_specific_isoforms = []
    for isoform in all_male_isoforms:
        if isoform in all_female_isoforms:
            shared_isoforms.append(isoform)
        elif isoform not in all_female_isoforms:
            male_specific_isoforms.append(isoform)
    for iso in all_female_isoforms:
        if iso not in all_male_isoforms:
            female_specific_isoforms.append(iso)
    return shared_isoforms, male_specific_isoforms, female_specific_isoforms

#write isoforms to files
def write_shared_isoforms():
    shared_isoforms, male_specific_isoforms, female_specific_isoforms = split_isoforms()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for iso in shared_isoforms:
            final = "%s\n" % str(iso)
            out.write(final)


def write_male_specific_isoforms():
    shared_isoforms, male_specific_isoforms, female_specific_isoforms = split_isoforms()
    output = sys.argv[4]
    with open(output, 'a') as out:
        for iso in male_specific_isoforms:
            final = "%s\n" % str(iso)
            out.write(final)


def write_female_specific_isoforms():
    shared_isoforms, male_specific_isoforms, female_specific_isoforms = split_isoforms()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for iso in female_specific_isoforms:
            final = "%s\n" % str(iso)
            out.write(final)


#call all functions
def call():
    shared = write_shared_isoforms()
    male_specific = write_male_specific_isoforms()
    female_specific = write_female_specific_isoforms()

call()

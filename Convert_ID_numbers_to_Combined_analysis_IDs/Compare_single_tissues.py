#comparing isoforms between single tissue analyses
#this will be easier that before as all of the single tissue analyses have the combined sexes isoform ID added to them, so I can just compare IDs rather than compare sequences, locations, etc.
#to run script: python3 Compare_single_tissues.py <single tissue classification file female liver> <single tissue classification file male liver> <single tissue classification file female brain> <single tissue classification file male brain> <single tissue classification file female pronephros> <single tissue classification file male pronephros> <single tissue classification file ovary> <single tissue classification file testis>
#have 8 single tissue files for this analysis

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


#summary read out of total number of isoforms per tissue:
def total_isoform_counts():
    fl_iso = read_female_liver_class()
    ml_iso = read_male_liver_class()
    fb_iso = read_female_brain_class()
    mb_iso = read_male_brain_class()
    fp_iso = read_female_pronephros_class()
    mp_iso = read_male_pronephros_class()
    o_iso = read_ovary_class()
    t_iso = read_testis_class()
    print("Total number of Isoforms in Female Liver")
    print(len(fl_iso))
    print("Total number of Isoforms in Male Liver")
    print(len(ml_iso))
    print("Total number of Isoforms in Female Brain")
    print(len(fb_iso))
    print("Total number of Isoforms in Male Brain")
    print(len(mb_iso))
    print("Total number of Isoforms in Female Pronephros")
    print(len(fp_iso))
    print("Total number of Isoforms in Male Pronephros")
    print(len(mp_iso))
    print("Total number of Isoforms in Ovary")
    print(len(o_iso))
    print("Total number of Isoforms in Testis")
    print(len(t_iso))



#need to compare isoforms
def compare():
    fl_iso = set(read_female_liver_class())
    ml_iso = set(read_male_liver_class())
    fb_iso = set(read_female_brain_class())
    mb_iso = set(read_male_brain_class())
    fp_iso = set(read_female_pronephros_class())
    mp_iso = set(read_male_pronephros_class())
    o_iso = set(read_ovary_class())
    t_iso = set(read_testis_class())
    #set intersection for all tissues/sexes
    shared_between_all_samples = fl_iso.intersection(ml_iso, fb_iso, mb_iso, fp_iso, mp_iso, o_iso, t_iso)
    print(len(shared_between_all_samples))
    #set intersections for all liver samples
    shared_liver = fl_iso.intersection(ml_iso)
    print(len(shared_liver))
    shared_brain = fb_iso.intersection(mb_iso)
    print(len(shared_brain))
    shared_pronephros = fp_iso.intersection(mp_iso)
    print(len(shared_pronephros))
    shared_gonad = o_iso.intersection(t_iso)
    print(len(shared_gonad))
    #set intersection for all tissues excpet gonads
    shared_somatic_tissues = fl_iso.intersection(ml_iso, fb_iso, mb_iso, fp_iso, mp_iso)
    print(len(shared_somatic_tissues))
    #set intersectionf for each sex
    female_shared_tissues = fl_iso.intersection(fb_iso, fp_iso, o_iso)
    print(len(female_shared_tissues))
    male_shared_tissues = ml_iso.intersection(mb_iso, mp_iso, t_iso)
    print(len(male_shared_tissues))



#call all functions
def call():
    isoforms_counts = total_isoform_counts()
    compare_function = compare()

call()

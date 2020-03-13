#comparing isoforms between single tissue analyses
#this will be easier that before as all of the single tissue analyses have the combined sexes isoform ID added to them, so I can just compare IDs rather than compare sequences, locations, etc.
#will use exon counts file for input for single tissues to remove any converted isoform ids that don't match up correctly (Gene_Isoform_Counts.py does this)
#to run script: python3 Compare_single_tissues.py <single tissue exon counts file female liver> <single tissue exon file male liver> <single tissue exon file female brain> <single tissue exon file male brain> <single tissue exon file female pronephros> <single tissue exon file male pronephros> <single tissue exon file ovary> <single tissue classification file testis>
#have 8 single tissue files for this analysis

import sys

#read in exon counts files for each single tissue
#return lists of isoform IDs
def read_female_liver_file():
    input_file = sys.argv[1]
    female_liver_isoforms = []
    with open(input_file, 'r') as fl_info:
        for line in fl_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_liver_isoforms.append(isoform)
    return female_liver_isoforms

def read_male_liver_file():
    input_file = sys.argv[2]
    male_liver_isoforms = []
    with open(input_file, 'r') as ml_info:
        for line in ml_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_liver_isoforms.append(isoform)
    return male_liver_isoforms


def read_female_brain_file():
    input_file = sys.argv[3]
    female_brain_isoforms = []
    with open(input_file, 'r') as fb_info:
        for line in fb_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_brain_isoforms.append(isoform)
    return female_brain_isoforms

def read_male_brain_file():
    input_file = sys.argv[4]
    male_brain_isoforms = []
    with open(input_file, 'r') as mb_info:
        for line in mb_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_brain_isoforms.append(isoform)
    return male_brain_isoforms


def read_female_pronephros_file():
    input_file = sys.argv[5]
    female_pronephros_isoforms = []
    with open(input_file, 'r') as fp_info:
        for line in fp_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_pronephros_isoforms.append(isoform)
    return female_pronephros_isoforms

def read_male_pronephros_file():
    input_file = sys.argv[6]
    male_pronephros_isoforms = []
    with open(input_file, 'r') as mp_info:
        for line in mp_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_pronephros_isoforms.append(isoform)
    return male_pronephros_isoforms


def read_ovary_file():
    input_file = sys.argv[7]
    ovary_isoforms = []
    with open(input_file, 'r') as o_info:
        for line in o_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                ovary_isoforms.append(isoform)
    return ovary_isoforms

def read_testis_file():
    input_file = sys.argv[8]
    testis_isoforms = []
    with open(input_file, 'r') as t_info:
        for line in t_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                testis_isoforms.append(isoform)
    return testis_isoforms


#summary read out of total number of isoforms per tissue:
def total_isoform_counts():
    fl_iso = read_female_liver_file()
    ml_iso = read_male_liver_file()
    fb_iso = read_female_brain_file()
    mb_iso = read_male_brain_file()
    fp_iso = read_female_pronephros_file()
    mp_iso = read_male_pronephros_file()
    o_iso = read_ovary_file()
    t_iso = read_testis_file()
    print("Total number of Isoforms in Female Liver")
    print(len(set(fl_iso)))
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
    fl_iso = set(read_female_liver_file())
    ml_iso = set(read_male_liver_file())
    fb_iso = set(read_female_brain_file())
    mb_iso = set(read_male_brain_file())
    fp_iso = set(read_female_pronephros_file())
    mp_iso = set(read_male_pronephros_file())
    o_iso = set(read_ovary_file())
    t_iso = set(read_testis_file())
    print("set intersection for all tissues/sexes")
    shared_between_all_samples = fl_iso.intersection(ml_iso, fb_iso, mb_iso, fp_iso, mp_iso, o_iso, t_iso)
    print(len(shared_between_all_samples))
    print("set intersections for liver samples")
    shared_liver = fl_iso.intersection(ml_iso)
    print(len(shared_liver))
    print("set intersections for brain samples")
    shared_brain = fb_iso.intersection(mb_iso)
    print(len(shared_brain))
    print("set intersections for pronephros samples")
    shared_pronephros = fp_iso.intersection(mp_iso)
    print(len(shared_pronephros))
    print("set intersections for gonads samples")
    shared_gonad = o_iso.intersection(t_iso)
    print(len(shared_gonad))
    '''print("set intersection for all tissues excpet gonads")
    shared_somatic_tissues = fl_iso.intersection(ml_iso, fb_iso, mb_iso, fp_iso, mp_iso)
    print(len(shared_somatic_tissues))
    print("set intersectionf for each sex")
    female_shared_tissues = fl_iso.intersection(fb_iso, fp_iso, o_iso)
    print(len(female_shared_tissues))
    male_shared_tissues = ml_iso.intersection(mb_iso, mp_iso, t_iso)
    print(len(male_shared_tissues))'''
    print("set intersection for liver and brain")
    shared_liver_brain = fl_iso.intersection(ml_iso, fb_iso, mb_iso)
    print(len(shared_liver_brain))
    print("set intersection for liver and pronephros")
    shared_liver_pronephros = fl_iso.intersection(ml_iso, fp_iso, mp_iso)
    print(len(shared_liver_pronephros))
    print("set intersection for liver and gonads")
    shared_liver_gonads = fl_iso.intersection(ml_iso, o_iso, t_iso)
    print(len(shared_liver_gonads))
    print("set intersection for brain and Pronephros")
    shared_brain_pronephros = fb_iso.intersection(mb_iso, fp_iso, mp_iso)
    print(len(shared_brain_pronephros))
    print("set intersection for brain and gonads")
    shared_brain_gonads = fb_iso.intersection(mb_iso, o_iso, t_iso)
    print(len(shared_brain_gonads))
    print("set intersection for pronephros and gonads")
    shared_pronephros_gonads = fp_iso.intersection(mp_iso, o_iso, t_iso)
    print(len(shared_pronephros_gonads))



#call all functions
def call():
    #isoforms_counts = total_isoform_counts()
    compare_function = compare()

call()

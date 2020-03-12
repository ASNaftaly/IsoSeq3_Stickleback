#Need to examine how many genes are shared between tissues from single tissue analysis
#currently have the the single tissue analyses blasted to combined sexes analysis to have the same isoform ID
#will need to pull the gene from the combined sexes classification file so the novel genes can be compared more easily
#to run script: python3 Compare_single_tissues_genes.py <combined sexes filtered classification file> <female liver converted classification file> <male liver converted_classification file> <female brain converted classification file> <male brain converted classification file> <female pronephros converted classification file> <male pronephros converted classification file> <ovary converted classification file> <testis converted classification file> <output shared all samples genes> <output shared somatic samples genes> <output shared female samples genes> <output shared male samples genes> <output shared liver samples genes> <output shared brain samples genes> <output shared pronephros samples genes> <output shared gonad samples genes>
#Author: Alice Naftaly, March 2020

import sys

#read in classification file and keep genesform id and gene id, transcript id
#returns dictionary with key = isoform id and value == [gene id, transcript id]
def read_combined_sexes_class():
    combined_class_file = sys.argv[1]
    class_dict = {}
    with open(combined_class_file, 'r') as combined_class:
        for line in combined_class:
            if line.startswith('PB'):
                new_line = line.split()
                isoform_id = new_line[0]
                gene_id = new_line[6]
                transcript_id = new_line[7]
                final_value = [gene_id, transcript_id]
                class_dict.update({isoform_id:final_value})
    return class_dict


#read in each converted single tissue classification file
#only need the isoform id from the combined sexes analysis (first value in line)
def read_female_liver_class():
    fl_file = sys.argv[2]
    fl_isoforms = []
    with open(fl_file, 'r') as fl:
        for line in fl:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                fl_isoforms.append(isoform_id)
    return fl_isoforms

def read_male_liver_class():
    ml_file = sys.argv[3]
    ml_isoforms = []
    with open(ml_file, 'r') as ml:
        for line in ml:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                ml_isoforms.append(isoform_id)
    return ml_isoforms

def read_female_brain_class():
    fb_file = sys.argv[4]
    fb_isoforms = []
    with open(fb_file, 'r') as fb:
        for line in fb:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                fb_isoforms.append(isoform_id)
    return fb_isoforms

def read_male_brain_class():
    mb_file = sys.argv[5]
    mb_isoforms = []
    with open(mb_file, 'r') as mb:
        for line in mb:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                mb_isoforms.append(isoform_id)
    return mb_isoforms

def read_female_pronephros_class():
    fp_file = sys.argv[6]
    fp_isoforms = []
    with open(fp_file, 'r') as fp:
        for line in fp:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                fp_isoforms.append(isoform_id)
    return fp_isoforms

def read_male_pronephros_class():
    mp_file = sys.argv[7]
    mp_isoforms = []
    with open(mp_file, 'r') as mp:
        for line in mp:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                mp_isoforms.append(isoform_id)
    return mp_isoforms

def read_ovary_class():
    o_file = sys.argv[8]
    o_isoforms = []
    with open(o_file, 'r') as o:
        for line in o:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                o_isoforms.append(isoform_id)
    return o_isoforms

def read_testis_class():
    t_file = sys.argv[9]
    t_isoforms = []
    with open(t_file, 'r') as t:
        for line in t:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                t_isoforms.append(isoform_id)
    return t_isoforms


#Now going to convert the isoform ids for each tissue into list of gene ids
def convert_fl():
    combined_sexes = read_combined_sexes_class()
    female_liver_isoforms = read_female_liver_class()
    female_liver_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in female_liver_isoforms:
            #adds gene id to new list
            female_liver_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_female_liver_genes = set(female_liver_genes)
    '''print("Total number of female liver isoforms")
    print(len(female_liver_isoforms))
    print("Total number of genes in female liver sample")
    print(len(set_female_liver_genes))'''
    return set_female_liver_genes

def convert_ml():
    combined_sexes = read_combined_sexes_class()
    male_liver_isoforms = read_male_liver_class()
    male_liver_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in male_liver_isoforms:
            #adds gene id to new list
            male_liver_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_male_liver_genes = set(male_liver_genes)
    '''print("Total number of male liver isoforms")
    print(len(male_liver_isoforms))
    print("Total number of genes in male liver sample")
    print(len(set_male_liver_genes))'''
    return set_male_liver_genes

def convert_fb():
    combined_sexes = read_combined_sexes_class()
    female_brain_isoforms = read_female_brain_class()
    female_brain_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in female_brain_isoforms:
            #adds gene id to new list
            female_brain_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_female_brain_genes = set(female_brain_genes)
    '''print("Total number of female brain isoforms")
    print(len(female_brain_isoforms))
    print("Total number of genes in female brain sample")
    print(len(set_female_brain_genes))'''
    return set_female_brain_genes

def convert_mb():
    combined_sexes = read_combined_sexes_class()
    male_brain_isoforms = read_male_brain_class()
    male_brain_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in male_brain_isoforms:
            #adds gene id to new list
            male_brain_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_male_brain_genes = set(male_brain_genes)
    '''print("Total number of male brain isoforms")
    print(len(male_brain_isoforms))
    print("Total number of genes in male brain sample")
    print(len(set_male_brain_genes))'''
    return set_male_brain_genes

def convert_fp():
    combined_sexes = read_combined_sexes_class()
    female_pronephros_isoforms = read_female_pronephros_class()
    female_pronephros_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in female_pronephros_isoforms:
            #adds gene id to new list
            female_pronephros_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_female_pronephros_genes = set(female_pronephros_genes)
    '''print("Total number of female pronephros isoforms")
    print(len(female_pronephros_isoforms))
    print("Total number of genes in female pronephros sample")
    print(len(set_female_pronephros_genes))'''
    return set_female_pronephros_genes

def convert_mp():
    combined_sexes = read_combined_sexes_class()
    male_pronephros_isoforms = read_male_pronephros_class()
    male_pronephros_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in male_pronephros_isoforms:
            #adds gene id to new list
            male_pronephros_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_male_pronephros_genes = set(male_pronephros_genes)
    '''print("Total number of male pronephros isoforms")
    print(len(male_pronephros_isoforms))
    print("Total number of genes in male pronephros sample")
    print(len(set_male_pronephros_genes))'''
    return set_male_pronephros_genes

def convert_o():
    combined_sexes = read_combined_sexes_class()
    female_ovary_isoforms = read_ovary_class()
    female_ovary_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in female_ovary_isoforms:
            #adds gene id to new list
            female_ovary_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_female_ovary_genes = set(female_ovary_genes)
    '''print("Total number of female ovary isoforms")
    print(len(female_ovary_isoforms))
    print("Total number of genes in ovary sample")
    print(len(set_female_ovary_genes))'''
    return set_female_ovary_genes

def convert_t():
    combined_sexes = read_combined_sexes_class()
    male_testis_isoforms = read_testis_class()
    male_testis_genes = []
    for key in combined_sexes:
        single_combined_sexes = combined_sexes[key]
        if key in male_testis_isoforms:
            #adds gene id to new list
            male_testis_genes.append(single_combined_sexes[0])
    #will set list to remove any duplicates
    set_male_testis_genes = set(male_testis_genes)
    '''print("Total number of male testis isoforms")
    print(len(male_testis_isoforms))
    print("Total number of genes in testis sample")
    print(len(set_male_testis_genes))'''
    return set_male_testis_genes


#now to compare genes between tissues
def compare_tissues():
    fl_genes = convert_fl()
    ml_genes = convert_ml()
    fb_genes = convert_fb()
    mb_genes = convert_mb()
    fp_genes = convert_fp()
    mp_genes = convert_mp()
    o_genes = convert_o()
    t_genes = convert_t()
    '''test_gene = "ENSGACG00000015771"
    if test_gene in fl_genes:
        print("Yes")
    if test_gene in ml_genes:
        print("Yes")
    if test_gene in fb_genes:
        print("Yes")
    if test_gene in mb_genes:
        print("Yes")
    if test_gene in fp_genes:
        print("Yes")
    if test_gene in mp_genes:
        print("Yes")
    if test_gene in o_genes:
        print("Yes")
    if test_gene in t_genes:
        print("Yes")'''
    #set intersection for all tissues/sexes
    shared_between_all_samples = fl_genes.intersection(ml_genes, fb_genes, mb_genes, fp_genes, mp_genes, o_genes, t_genes)
    print(len(shared_between_all_samples))
    #set intersection for all tissues excpet gonads
    shared_somatic_tissues = fl_genes.intersection(ml_genes, fb_genes, mb_genes, fp_genes, mp_genes)
    print(len(shared_somatic_tissues))
    #set intersectionf for each sex
    female_shared_tissues = fl_genes.intersection(fb_genes, fp_genes, o_genes)
    print(len(female_shared_tissues))
    male_shared_tissues = ml_genes.intersection(mb_genes, mp_genes, t_genes)
    print(len(male_shared_tissues))
    #set intersections for all liver samples
    shared_liver = fl_genes.intersection(ml_genes)
    print(len(shared_liver))
    shared_brain = fb_genes.intersection(mb_genes)
    print(len(shared_brain))
    shared_pronephros = fp_genes.intersection(mp_genes)
    print(len(shared_pronephros))
    shared_gonad = o_genes.intersection(t_genes)
    print(len(shared_gonad))
    return shared_between_all_samples, shared_somatic_tissues, female_shared_tissues, male_shared_tissues, shared_liver, shared_brain, shared_pronephros, shared_gonad


#write genes to files for each category

def write_shared_all_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[10]
    with open(output, 'a') as out:
        for gene in shared_all_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_somatic_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[11]
    with open(output, 'a') as out:
        for gene in shared_somatic_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_female_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[12]
    with open(output, 'a') as out:
        for gene in shared_female_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_male_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[13]
    with open(output, 'a') as out:
        for gene in shared_male_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_liver_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[14]
    with open(output, 'a') as out:
        for gene in shared_liver_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_brain_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[15]
    with open(output, 'a') as out:
        for gene in shared_brain_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_pronephros_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[16]
    with open(output, 'a') as out:
        for gene in shared_pronephros_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_gonad_samples():
    shared_all_samples, shared_somatic_samples, shared_female_samples, shared_male_samples, shared_liver_samples, shared_brain_samples, shared_pronephros_samples, shared_gonad_samples = compare_tissues()
    output = sys.argv[17]
    with open(output, 'a') as out:
        for gene in shared_gonad_samples:
            final = "%s\n" % str(gene)
            out.write(final)


#call all functions
def call():
    shared_all_samples = write_shared_all_samples()
    shared_somatic = write_shared_somatic_samples()
    shared_female = write_shared_female_samples()
    shared_male = write_shared_male_samples()
    shared_liver = write_shared_liver_samples()
    shared_brain = write_shared_brain_samples()
    shared_pronephros = write_shared_pronephros_samples()
    shared_gonad = write_shared_gonad_samples()

call()

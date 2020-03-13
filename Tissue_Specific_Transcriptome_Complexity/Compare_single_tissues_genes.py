#Need to examine how many genes are shared between tissues from single tissue analysis
#currently have the the single tissue analyses blasted to combined sexes analysis to have the same form ID
#will use isoform counts file for input for single tissues to remove any converted isoform ids that don't match up correctly (Gene_Isoform_Counts.py does this)
#to run script: python3 Compare_single_tissues_genes.py <female liver isoform file> <male liver isoform file> <female brain isoform file> <male brain isoform file> <female pronephros isoform file> <male pronephros isoform file> <ovary isoform file> <testis isoform file> <output shared all samples genes> <output shared liver samples genes> <output shared brain samples genes> <output shared pronephros samples genes> <output shared gonad samples genes> <output shared liver & brain genes> <output shared liver & pronephros genes> <output shared liver & gonad genes> <output shared brain & pronephros genes> <output shared brain & gonad genes> <output shared pronephros & gonad genes>
#Author: Alice Naftaly, March 2020

import sys


#read in each isoform counts file
#returns a list of gene ids
def read_female_liver_file():
    input_file = sys.argv[1]
    female_liver_genes = []
    with open(input_file, 'r') as fl_info:
        for line in fl_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                female_liver_genes.append(gene)
    return female_liver_genes

def read_male_liver_file():
    input_file = sys.argv[2]
    male_liver_genes = []
    with open(input_file, 'r') as ml_info:
        for line in ml_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                male_liver_genes.append(gene)
    return male_liver_genes


def read_female_brain_file():
    input_file = sys.argv[3]
    female_brain_genes = []
    with open(input_file, 'r') as fb_info:
        for line in fb_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                female_brain_genes.append(gene)
    return female_brain_genes

def read_male_brain_file():
    input_file = sys.argv[4]
    male_brain_genes = []
    with open(input_file, 'r') as mb_info:
        for line in mb_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                male_brain_genes.append(gene)
    return male_brain_genes


def read_female_pronephros_file():
    input_file = sys.argv[5]
    female_pronephros_genes = []
    with open(input_file, 'r') as fp_info:
        for line in fp_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                female_pronephros_genes.append(gene)
    return female_pronephros_genes

def read_male_pronephros_file():
    input_file = sys.argv[6]
    male_pronephros_genes = []
    with open(input_file, 'r') as mp_info:
        for line in mp_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                male_pronephros_genes.append(gene)
    return male_pronephros_genes

def read_ovary_file():
    input_file = sys.argv[7]
    ovary_genes = []
    with open(input_file, 'r') as o_info:
        for line in o_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                ovary_genes.append(gene)
    return ovary_genes

def read_testis_file():
    input_file = sys.argv[8]
    testis_genes = []
    with open(input_file, 'r') as t_info:
        for line in t_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                testis_genes.append(gene)
    return testis_genes

def number_of_genes():
    fl_final_genes = set(read_female_liver_file())
    ml_final_genes = set(read_male_liver_file())
    fb_final_genes = set(read_female_brain_file())
    mb_final_genes = set(read_male_brain_file())
    fp_final_genes = set(read_female_pronephros_file())
    mp_final_genes = set(read_male_pronephros_file())
    o_final_genes = set(read_ovary_file())
    t_final_genes = set(read_testis_file())
    print("Number of Female Liver Genes")
    print(len(fl_final_genes))
    print("Number of Male Liver Genes")
    print(len(ml_final_genes))
    print("Number of Female Brain Genes")
    print(len(fb_final_genes))
    print("Number of Male Brain Genes")
    print(len(mb_final_genes))
    print("Number of Female Pronephros Genes")
    print(len(fp_final_genes))
    print("Number of Male Pronephros Genes")
    print(len(mp_final_genes))
    print("Number of Ovary Genes")
    print(len(o_final_genes))
    print("Number of Testis Genes")
    print(len(t_final_genes))
    return fl_final_genes, ml_final_genes, fb_final_genes, mb_final_genes, fp_final_genes, mp_final_genes, o_final_genes, t_final_genes


#now to compare genes between tissues
def compare_tissues():
    fl_genes, ml_genes, fb_genes, mb_genes, fp_genes, mp_genes, o_genes, t_genes = number_of_genes()
    print("set intersection for all tissues/sexes")
    shared_between_all_samples = fl_genes.intersection(ml_genes, fb_genes, mb_genes, fp_genes, mp_genes, o_genes, t_genes)
    print(len(shared_between_all_samples))
    print("set intersections for liver samples")
    shared_liver = fl_genes.intersection(ml_genes)
    print(len(shared_liver))
    print("set intersections for brain samples")
    shared_brain = fb_genes.intersection(mb_genes)
    print(len(shared_brain))
    print("set intersections for pronephros samples")
    shared_pronephros = fp_genes.intersection(mp_genes)
    print(len(shared_pronephros))
    print("set intersections for gonads samples")
    shared_gonad = o_genes.intersection(t_genes)
    print(len(shared_gonad))
    print("set intersection for liver and brain")
    shared_liver_brain = fl_genes.intersection(ml_genes, fb_genes, mb_genes)
    print(len(shared_liver_brain))
    print("set intersection for liver and pronephros")
    shared_liver_pronephros = fl_genes.intersection(ml_genes, fp_genes, mp_genes)
    print(len(shared_liver_pronephros))
    print("set intersection for liver and gonads")
    shared_liver_gonads = fl_genes.intersection(ml_genes, o_genes, t_genes)
    print(len(shared_liver_gonads))
    print("set intersection for brain and Pronephros")
    shared_brain_pronephros = fb_genes.intersection(mb_genes, fp_genes, mp_genes)
    print(len(shared_brain_pronephros))
    print("set intersection for brain and gonads")
    shared_brain_gonads = fb_genes.intersection(mb_genes, o_genes, t_genes)
    print(len(shared_brain_gonads))
    print("set intersection for pronephros and gonads")
    shared_pronephros_gonads = fp_genes.intersection(mp_genes, o_genes, t_genes)
    print(len(shared_pronephros_gonads))
    return shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads

#write genes to files for each category

def write_shared_all_samples():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[9]
    with open(output, 'a') as out:
        for gene in shared_between_all_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_liver_samples():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[10]
    with open(output, 'a') as out:
        for gene in shared_liver:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_brain_samples():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[11]
    with open(output, 'a') as out:
        for gene in shared_brain:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_pronephros_samples():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[12]
    with open(output, 'a') as out:
        for gene in shared_pronephros:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_gonad_samples():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[13]
    with open(output, 'a') as out:
        for gene in shared_gonad:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_liver_brain():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[14]
    with open(output, 'a') as out:
        for gene in shared_liver_brain:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_liver_pronephros():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[15]
    with open(output, 'a') as out:
        for gene in shared_liver_pronephros:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_liver_gonads():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[16]
    with open(output, 'a') as out:
        for gene in shared_liver_gonads:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_brain_pronephros():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[17]
    with open(output, 'a') as out:
        for gene in shared_brain_pronephros:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_brain_gonads():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[18]
    with open(output, 'a') as out:
        for gene in shared_brain_gonads:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_pronephros_gonads():
    shared_between_all_samples, shared_liver, shared_brain, shared_pronephros, shared_gonad, shared_liver_brain, shared_liver_pronephros, shared_liver_gonads, shared_brain_pronephros, shared_brain_gonads, shared_pronephros_gonads = compare_tissues()
    output = sys.argv[19]
    with open(output, 'a') as out:
        for gene in shared_pronephros_gonads:
            final = "%s\n" % str(gene)
            out.write(final)

#call all functions
def call():
    shared_all_samples = write_shared_all_samples()
    shared_liver = write_shared_liver_samples()
    shared_brain = write_shared_brain_samples()
    shared_pronephros = write_shared_pronephros_samples()
    shared_gonad = write_shared_gonad_samples()
    shared_liver_brain = write_shared_liver_brain()
    shared_liver_pronephros = write_shared_liver_pronephros()
    shared_liver_gonads = write_shared_liver_gonads()
    shared_brain_pronephros = write_shared_brain_pronephros()
    shared_brain_gonads = write_shared_brain_gonads()
    shared_pronephros_gonads = write_shared_pronephros_gonads()

call()

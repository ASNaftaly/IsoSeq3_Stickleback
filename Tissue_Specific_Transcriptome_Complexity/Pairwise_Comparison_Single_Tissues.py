#Pairwise comparison of shared genes and isoforms between each of the single tissues
#this script is a combination of Compare_single_tissues.py and Compare_single_tissues_genes.py; but only looks at pairwise comparison of all single tissues
#reads in all single tissue exon counts files and all single tissue isoform counts
#return counts of shared genes and shared isoforms
#will not create any output for now
#to run script: python3 Pairwise_Comparison_Single_Tissues.py < <single tissue exon counts file female liver> <single tissue exon file male liver> <single tissue exon file female brain> <single tissue exon file male brain> <single tissue exon file female pronephros> <single tissue exon file male pronephros> <single tissue exon file ovary> <single tissue classification file testis> <female liver isoform file> <male liver isoform file> <female brain isoform file> <male brain isoform file> <female pronephros isoform file> <male pronephros isoform file> <ovary isoform file> <testis isoform file>
#Author: Alice Naftaly, March 2020

import sys

#read in exon counts files for each single tissue
#return lists of isoform IDs
def read_female_liver_exons_file():
    input_file = sys.argv[1]
    female_liver_isoforms = []
    with open(input_file, 'r') as fl_info:
        for line in fl_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_liver_isoforms.append(isoform)
    return female_liver_isoforms

def read_male_liver_exons_file():
    input_file = sys.argv[2]
    male_liver_isoforms = []
    with open(input_file, 'r') as ml_info:
        for line in ml_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_liver_isoforms.append(isoform)
    return male_liver_isoforms


def read_female_brain_exons_file():
    input_file = sys.argv[3]
    female_brain_isoforms = []
    with open(input_file, 'r') as fb_info:
        for line in fb_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_brain_isoforms.append(isoform)
    return female_brain_isoforms

def read_male_brain_exons_file():
    input_file = sys.argv[4]
    male_brain_isoforms = []
    with open(input_file, 'r') as mb_info:
        for line in mb_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_brain_isoforms.append(isoform)
    return male_brain_isoforms


def read_female_pronephros_exons_file():
    input_file = sys.argv[5]
    female_pronephros_isoforms = []
    with open(input_file, 'r') as fp_info:
        for line in fp_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_pronephros_isoforms.append(isoform)
    return female_pronephros_isoforms

def read_male_pronephros_exons_file():
    input_file = sys.argv[6]
    male_pronephros_isoforms = []
    with open(input_file, 'r') as mp_info:
        for line in mp_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_pronephros_isoforms.append(isoform)
    return male_pronephros_isoforms


def read_ovary_exons_file():
    input_file = sys.argv[7]
    ovary_isoforms = []
    with open(input_file, 'r') as o_info:
        for line in o_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                ovary_isoforms.append(isoform)
    return ovary_isoforms

def read_testis_exons_file():
    input_file = sys.argv[8]
    testis_isoforms = []
    with open(input_file, 'r') as t_info:
        for line in t_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                testis_isoforms.append(isoform)
    return testis_isoforms

#read in each isoform counts file
#returns a list of gene ids
def read_female_liver_isoforms_file():
    input_file = sys.argv[9]
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

def read_male_liver_isoforms_file():
    input_file = sys.argv[10]
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


def read_female_brain_isoforms_file():
    input_file = sys.argv[11]
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

def read_male_brain_isoforms_file():
    input_file = sys.argv[12]
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


def read_female_pronephros_isoforms_file():
    input_file = sys.argv[13]
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

def read_male_pronephros_isoforms_file():
    input_file = sys.argv[14]
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

def read_ovary_isoforms_file():
    input_file = sys.argv[15]
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

def read_testis_isoforms_file():
    input_file = sys.argv[16]
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

#genes breakdown
def pairwise_gene_comparison():
    fl_final_genes = set(read_female_liver_isoforms_file())
    ml_final_genes = set(read_male_liver_isoforms_file())
    fb_final_genes = set(read_female_brain_isoforms_file())
    mb_final_genes = set(read_male_brain_isoforms_file())
    fp_final_genes = set(read_female_pronephros_isoforms_file())
    mp_final_genes = set(read_male_pronephros_isoforms_file())
    o_final_genes = set(read_ovary_isoforms_file())
    t_final_genes = set(read_testis_isoforms_file())
    #liver comparison
    fl_ml_total_genes = len(fl_final_genes) + len(ml_final_genes)
    print("Total number of female and male liver genes, including shared genes")
    print(fl_ml_total_genes)
    print("Number of shared liver genes")
    fl_ml_shared_genes = fl_final_genes.intersection(ml_final_genes)
    print(len(fl_ml_shared_genes))
    print("\n")
    #brain comparison
    fb_mb_total_genes = len(fb_final_genes) + len(mb_final_genes)
    print("Total number of female and male brain genes, including shared genes")
    print(fb_mb_total_genes)
    print("Number of shared brain genes")
    fb_mb_shared_genes = fb_final_genes.intersection(mb_final_genes)
    print(len(fb_mb_shared_genes))
    print("\n")
    #pronephros comparison
    print("Total number of female and male pronephros genes, including shared genes")
    fp_mp_total_genes = len(fp_final_genes) + len(mp_final_genes)
    print(fp_mp_total_genes)
    print("Number of shared pronephros genes")
    fp_mp_shared_genes = fp_final_genes.intersection(mp_final_genes)
    print(len(fp_mp_shared_genes))
    print("\n")
    #gonads comparison
    print("Total number of female and male gonads genes, including shared genes")
    o_t_total_genes = len(o_final_genes) + len(t_final_genes)
    print(o_t_total_genes)
    print("Number of shared gonads genes")
    o_t_shared_genes = o_final_genes.intersection(t_final_genes)
    print(len(o_t_shared_genes))
    print("\n")
    #female liver vs female brain
    print("Total number of female liver and female brain genes, including shared genes")
    fl_fb_total_genes = len(fl_final_genes) + len(fb_final_genes)
    print(fl_fb_total_genes)
    print("Number of shared female liver & female brain genes")
    fl_fb_shared_genes = fl_final_genes.intersection(fb_final_genes)
    print(len(fl_fb_shared_genes))
    print("\n")
    #female liver vs female pronephros
    print("Total number of female liver and female pronephros genes, including shared genes")
    fl_fp_total_genes = len(fl_final_genes) + len(fp_final_genes)
    print(fl_fp_total_genes)
    print("Number of shared female liver & female pronephros genes")
    fl_fp_shared_genes = fl_final_genes.intersection(fp_final_genes)
    print(len(fl_fp_shared_genes))
    print("\n")
    #female liver vs ovary
    print("Total number of female liver and ovary genes, including shared genes")
    fl_o_total_genes = len(fl_final_genes) + len(o_final_genes)
    print(fl_o_total_genes)
    print("Number of shared female liver & ovary genes")
    fl_o_shared_genes = fl_final_genes.intersection(o_final_genes)
    print(len(fl_o_shared_genes))
    print("\n")
    #female liver vs male brain
    print("Total number of female liver and male brain genes, including shared genes")
    fl_mb_total_genes = len(fl_final_genes) + len(mb_final_genes)
    print(fl_mb_total_genes)
    print("Number of shared female liver & male brain genes")
    fl_mb_shared_genes = fl_final_genes.intersection(mb_final_genes)
    print(len(fl_mb_shared_genes))
    print("\n")
    #female liver vs male pronephros
    print("Total number of female liver and male pronephros genes, including shared genes")
    fl_mp_total_genes = len(fl_final_genes) + len(mp_final_genes)
    print(fl_mp_total_genes)
    print("Number of shared female liver & male pronephros genes")
    fl_mp_shared_genes = fl_final_genes.intersection(mp_final_genes)
    print(len(fl_mp_shared_genes))
    print("\n")
    #female liver vs testis
    print("Total number of female liver and testis genes, including shared genes")
    fl_t_total_genes = len(fl_final_genes) + len(t_final_genes)
    print(fl_t_total_genes)
    print("Number of shared female liver & testis genes")
    fl_t_shared_genes = fl_final_genes.intersection(t_final_genes)
    print(len(fl_t_shared_genes))
    print("\n")
    #male liver vs male brain
    print("Total number of male liver and male brain genes, including shared genes")
    ml_mb_total_genes = len(ml_final_genes) + len(mb_final_genes)
    print(ml_mb_total_genes)
    print("Number of shared male liver & male brain genes")
    ml_mb_shared_genes = ml_final_genes.intersection(mb_final_genes)
    print(len(ml_mb_shared_genes))
    print("\n")
    #male liver vs male pronephros
    print("Total number of male liver and male pronephros genes, including shared genes")
    ml_mp_total_genes = len(ml_final_genes) + len(mp_final_genes)
    print(ml_mp_total_genes)
    print("Number of shared male liver & male pronephros genes")
    ml_mp_shared_genes = ml_final_genes.intersection(mp_final_genes)
    print(len(ml_mp_shared_genes))
    print("\n")
    #male liver vs testis
    print("Total number of male liver and testis genes, including shared genes")
    ml_t_total_genes = len(ml_final_genes) + len(t_final_genes)
    print(ml_t_total_genes)
    print("Number of shared male liver & testis genes")
    ml_t_shared_genes = ml_final_genes.intersection(t_final_genes)
    print(len(ml_t_shared_genes))
    print("\n")
    #male liver vs female brain
    print("Total number of male liver and female brain genes, including shared genes")
    ml_fb_total_genes = len(ml_final_genes) + len(fb_final_genes)
    print(ml_fb_total_genes)
    print("Number of shared male liver & female brain genes")
    ml_fb_shared_genes = ml_final_genes.intersection(fb_final_genes)
    print(len(ml_fb_shared_genes))
    print("\n")
    #male liver vs female pronephros
    print("Total number of male liver and female pronephros genes, including shared genes")
    ml_fp_total_genes = len(ml_final_genes) + len(fp_final_genes)
    print(ml_fp_total_genes)
    print("Number of shared male liver & female pronephros genes")
    ml_fp_shared_genes = ml_final_genes.intersection(fp_final_genes)
    print(len(ml_fp_shared_genes))
    print("\n")
    #male liver vs ovary
    print("Total number of male liver and ovary genes, including shared genes")
    ml_o_total_genes = len(ml_final_genes) + len(o_final_genes)
    print(ml_o_total_genes)
    print("Number of shared male liver & ovary genes")
    ml_o_shared_genes = ml_final_genes.intersection(o_final_genes)
    print(len(ml_o_shared_genes))
    print("\n")
    #female brain vs female pronephros
    print("Total number of female brain and female pronephros genes, including shared genes")
    fb_fp_total_genes = len(fb_final_genes) + len(fp_final_genes)
    print(fb_fp_total_genes)
    print("Number of shared female brain & female pronephros genes")
    fb_fp_shared_genes = fb_final_genes.intersection(fp_final_genes)
    print(len(fb_fp_shared_genes))
    print("\n")
    #female brain vs ovary
    print("Total number of female brain and ovary genes, including shared genes")
    fb_o_total_genes = len(fb_final_genes) + len(o_final_genes)
    print(fb_o_total_genes)
    print("Number of shared female brain & ovary genes")
    fb_o_shared_genes = fb_final_genes.intersection(o_final_genes)
    print(len(fb_o_shared_genes))
    print("\n")
    #female brain vs male pronephros
    print("Total number of female brain and male pronephros genes, including shared genes")
    fb_mp_total_genes = len(fb_final_genes) + len(mp_final_genes)
    print(fb_mp_total_genes)
    print("Number of shared female brain & male pronephros genes")
    fb_mp_shared_genes = fb_final_genes.intersection(mp_final_genes)
    print(len(fb_mp_shared_genes))
    print("\n")
    #female brain vs testis
    print("Total number of female brain and testis genes, including shared genes")
    fb_t_total_genes = len(fb_final_genes) + len(t_final_genes)
    print(fb_t_total_genes)
    print("Number of shared female brain & testis genes")
    fb_t_shared_genes = fb_final_genes.intersection(t_final_genes)
    print(len(fb_t_shared_genes))
    print("\n")
    #male brain vs male pronephros
    print("Total number of male brain and male pronephros, including shared genes")
    mb_mp_total_genes = len(mb_final_genes) + len(mp_final_genes)
    print(mb_mp_total_genes)
    print("Number of shared male brain & male pronephros genes")
    mb_mp_shared_genes = mb_final_genes.intersection(mp_final_genes)
    print(len(mb_mp_shared_genes))
    print("\n")
    #male brain vs testis
    print("Total number of male brain and testis, including shared genes")
    mb_t_total_genes = len(mb_final_genes) + len(t_final_genes)
    print(mb_t_total_genes)
    print("Number of shared male brain & testis genes")
    mb_t_shared_genes = mb_final_genes.intersection(t_final_genes)
    print(len(mb_t_shared_genes))
    print("\n")
    #male brain vs female pronephros
    print("Total number of male brain and female pronephros, including shared genes")
    mb_fp_total_genes = len(mb_final_genes) + len(fp_final_genes)
    print(mb_fp_total_genes)
    print("Number of shared male brain & female pronephros genes")
    mb_fp_shared_genes = mb_final_genes.intersection(fp_final_genes)
    print(len(mb_fp_shared_genes))
    print("\n")
    #male brain vs ovary
    print("Total number of male brain and ovary, including shared genes")
    mb_o_total_genes = len(mb_final_genes) + len(o_final_genes)
    print(mb_o_total_genes)
    print("Number of shared male brain & ovary genes")
    mb_o_shared_genes = mb_final_genes.intersection(o_final_genes)
    print(len(mb_o_shared_genes))
    print("\n")
    #female pronephros vs ovary
    print("Total number of female pronephros and ovary, including shared genes")
    fp_o_total_genes = len(fp_final_genes) + len(o_final_genes)
    print(fp_o_total_genes)
    print("Number of shared female pronephros & ovary genes")
    fp_o_shared_genes = fp_final_genes.intersection(o_final_genes)
    print(len(fp_o_shared_genes))
    print("\n")
    #female pronephros vs testis
    print("Total number of female pronephros and testis, including shared genes")
    fp_t_total_genes = len(fp_final_genes) + len(t_final_genes)
    print(fp_t_total_genes)
    print("Number of shared female pronephros & testis genes")
    fp_t_shared_genes = fp_final_genes.intersection(t_final_genes)
    print(len(fp_t_shared_genes))
    print("\n")
    #male pronephros vs testis
    print("Total number of male pronephros and testis, including shared genes")
    mp_t_total_genes = len(mp_final_genes) + len(t_final_genes)
    print(mp_t_total_genes)
    print("Number of shared male pronephros & testis genes")
    mp_t_shared_genes = mp_final_genes.intersection(t_final_genes)
    print(len(mp_t_shared_genes))
    print("\n")
    #male pronephros vs ovary
    print("Total number of male pronephros and ovary, including shared genes")
    mp_o_total_genes = len(mp_final_genes) + len(o_final_genes)
    print(mp_o_total_genes)
    print("Number of shared male pronephros & ovary genes")
    mp_o_shared_genes = mp_final_genes.intersection(o_final_genes)
    print(len(mp_o_shared_genes))
    print("\n")



def pairwise_isoform_comparison():
    fl_final_isoforms = set(read_female_liver_exons_file())
    ml_final_isoforms = set(read_male_liver_exons_file())
    fb_final_isoforms = set(read_female_brain_exons_file())
    mb_final_isoforms = set(read_male_brain_exons_file())
    fp_final_isoforms = set(read_female_pronephros_exons_file())
    mp_final_isoforms = set(read_male_pronephros_exons_file())
    o_final_isoforms = set(read_ovary_exons_file())
    t_final_isoforms = set(read_testis_exons_file())
    #liver comparison
    fl_ml_total_isoforms = len(fl_final_isoforms) + len(ml_final_isoforms)
    print("Total number of female and male liver isoforms, including shared isoforms")
    print(fl_ml_total_isoforms)
    print("Number of shared liver isoforms")
    fl_ml_shared_isoforms = fl_final_isoforms.intersection(ml_final_isoforms)
    print(len(fl_ml_shared_isoforms))
    print("\n")
    #brain comparison
    fb_mb_total_isoforms = len(fb_final_isoforms) + len(mb_final_isoforms)
    print("Total number of female and male brain isoforms, including shared isoforms")
    print(fb_mb_total_isoforms)
    print("Number of shared brain isoforms")
    fb_mb_shared_isoforms = fb_final_isoforms.intersection(mb_final_isoforms)
    print(len(fb_mb_shared_isoforms))
    print("\n")
    #pronephros comparison
    print("Total number of female and male pronephros isoforms, including shared isoforms")
    fp_mp_total_isoforms = len(fp_final_isoforms) + len(mp_final_isoforms)
    print(fp_mp_total_isoforms)
    print("Number of shared pronephros isoforms")
    fp_mp_shared_isoforms = fp_final_isoforms.intersection(mp_final_isoforms)
    print(len(fp_mp_shared_isoforms))
    print("\n")
    #gonads comparison
    print("Total number of female and male gonads isoforms, including shared isoforms")
    o_t_total_isoforms = len(o_final_isoforms) + len(t_final_isoforms)
    print(o_t_total_isoforms)
    print("Number of shared gonads isoforms")
    o_t_shared_isoforms = o_final_isoforms.intersection(t_final_isoforms)
    print(len(o_t_shared_isoforms))
    print("\n")
    #female liver vs female brain
    print("Total number of female liver and female brain isoforms, including shared isoforms")
    fl_fb_total_isoforms = len(fl_final_isoforms) + len(fb_final_isoforms)
    print(fl_fb_total_isoforms)
    print("Number of shared female liver & female brain isoforms")
    fl_fb_shared_isoforms = fl_final_isoforms.intersection(fb_final_isoforms)
    print(len(fl_fb_shared_isoforms))
    print("\n")
    #female liver vs female pronephros
    print("Total number of female liver and female pronephros isoforms, including shared isoforms")
    fl_fp_total_isoforms = len(fl_final_isoforms) + len(fp_final_isoforms)
    print(fl_fp_total_isoforms)
    print("Number of shared female liver & female pronephros isoforms")
    fl_fp_shared_isoforms = fl_final_isoforms.intersection(fp_final_isoforms)
    print(len(fl_fp_shared_isoforms))
    print("\n")
    #female liver vs ovary
    print("Total number of female liver and ovary isoforms, including shared isoforms")
    fl_o_total_isoforms = len(fl_final_isoforms) + len(o_final_isoforms)
    print(fl_o_total_isoforms)
    print("Number of shared female liver & ovary isoforms")
    fl_o_shared_isoforms = fl_final_isoforms.intersection(o_final_isoforms)
    print(len(fl_o_shared_isoforms))
    print("\n")
    #female liver vs male brain
    print("Total number of female liver and male brain isoforms, including shared isoforms")
    fl_mb_total_isoforms = len(fl_final_isoforms) + len(mb_final_isoforms)
    print(fl_mb_total_isoforms)
    print("Number of shared female liver & male brain isoforms")
    fl_mb_shared_isoforms = fl_final_isoforms.intersection(mb_final_isoforms)
    print(len(fl_mb_shared_isoforms))
    print("\n")
    #female liver vs male pronephros
    print("Total number of female liver and male pronephros isoforms, including shared isoforms")
    fl_mp_total_isoforms = len(fl_final_isoforms) + len(mp_final_isoforms)
    print(fl_mp_total_isoforms)
    print("Number of shared female liver & male pronephros isoforms")
    fl_mp_shared_isoforms = fl_final_isoforms.intersection(mp_final_isoforms)
    print(len(fl_mp_shared_isoforms))
    print("\n")
    #female liver vs testis
    print("Total number of female liver and testis isoforms, including shared isoforms")
    fl_t_total_isoforms = len(fl_final_isoforms) + len(t_final_isoforms)
    print(fl_t_total_isoforms)
    print("Number of shared female liver & testis isoforms")
    fl_t_shared_isoforms = fl_final_isoforms.intersection(t_final_isoforms)
    print(len(fl_t_shared_isoforms))
    print("\n")
    #male liver vs male brain
    print("Total number of male liver and male brain isoforms, including shared isoforms")
    ml_mb_total_isoforms = len(ml_final_isoforms) + len(mb_final_isoforms)
    print(ml_mb_total_isoforms)
    print("Number of shared male liver & male brain isoforms")
    ml_mb_shared_isoforms = ml_final_isoforms.intersection(mb_final_isoforms)
    print(len(ml_mb_shared_isoforms))
    print("\n")
    #male liver vs male pronephros
    print("Total number of male liver and male pronephros isoforms, including shared isoforms")
    ml_mp_total_isoforms = len(ml_final_isoforms) + len(mp_final_isoforms)
    print(ml_mp_total_isoforms)
    print("Number of shared male liver & male pronephros isoforms")
    ml_mp_shared_isoforms = ml_final_isoforms.intersection(mp_final_isoforms)
    print(len(ml_mp_shared_isoforms))
    print("\n")
    #male liver vs testis
    print("Total number of male liver and testis isoforms, including shared isoforms")
    ml_t_total_isoforms = len(ml_final_isoforms) + len(t_final_isoforms)
    print(ml_t_total_isoforms)
    print("Number of shared male liver & testis isoforms")
    ml_t_shared_isoforms = ml_final_isoforms.intersection(t_final_isoforms)
    print(len(ml_t_shared_isoforms))
    print("\n")
    #male liver vs female brain
    print("Total number of male liver and female brain isoforms, including shared isoforms")
    ml_fb_total_isoforms = len(ml_final_isoforms) + len(fb_final_isoforms)
    print(ml_fb_total_isoforms)
    print("Number of shared male liver & female brain isoforms")
    ml_fb_shared_isoforms = ml_final_isoforms.intersection(fb_final_isoforms)
    print(len(ml_fb_shared_isoforms))
    print("\n")
    #male liver vs female pronephros
    print("Total number of male liver and female pronephros isoforms, including shared isoforms")
    ml_fp_total_isoforms = len(ml_final_isoforms) + len(fp_final_isoforms)
    print(ml_fp_total_isoforms)
    print("Number of shared male liver & female pronephros isoforms")
    ml_fp_shared_isoforms = ml_final_isoforms.intersection(fp_final_isoforms)
    print(len(ml_fp_shared_isoforms))
    print("\n")
    #male liver vs ovary
    print("Total number of male liver and ovary isoforms, including shared isoforms")
    ml_o_total_isoforms = len(ml_final_isoforms) + len(o_final_isoforms)
    print(ml_o_total_isoforms)
    print("Number of shared male liver & ovary isoforms")
    ml_o_shared_isoforms = ml_final_isoforms.intersection(o_final_isoforms)
    print(len(ml_o_shared_isoforms))
    print("\n")
    #female brain vs female pronephros
    print("Total number of female brain and female pronephros isoforms, including shared isoforms")
    fb_fp_total_isoforms = len(fb_final_isoforms) + len(fp_final_isoforms)
    print(fb_fp_total_isoforms)
    print("Number of shared female brain & female pronephros isoforms")
    fb_fp_shared_isoforms = fb_final_isoforms.intersection(fp_final_isoforms)
    print(len(fb_fp_shared_isoforms))
    print("\n")
    #female brain vs ovary
    print("Total number of female brain and ovary isoforms, including shared isoforms")
    fb_o_total_isoforms = len(fb_final_isoforms) + len(o_final_isoforms)
    print(fb_o_total_isoforms)
    print("Number of shared female brain & ovary isoforms")
    fb_o_shared_isoforms = fb_final_isoforms.intersection(o_final_isoforms)
    print(len(fb_o_shared_isoforms))
    print("\n")
    #female brain vs male pronephros
    print("Total number of female brain and male pronephros isoforms, including shared isoforms")
    fb_mp_total_isoforms = len(fb_final_isoforms) + len(mp_final_isoforms)
    print(fb_mp_total_isoforms)
    print("Number of shared female brain & male pronephros isoforms")
    fb_mp_shared_isoforms = fb_final_isoforms.intersection(mp_final_isoforms)
    print(len(fb_mp_shared_isoforms))
    print("\n")
    #female brain vs testis
    print("Total number of female brain and testis isoforms, including shared isoforms")
    fb_t_total_isoforms = len(fb_final_isoforms) + len(t_final_isoforms)
    print(fb_t_total_isoforms)
    print("Number of shared female brain & testis isoforms")
    fb_t_shared_isoforms = fb_final_isoforms.intersection(t_final_isoforms)
    print(len(fb_t_shared_isoforms))
    print("\n")
    #male brain vs male pronephros
    print("Total number of male brain and male pronephros, including shared isoforms")
    mb_mp_total_isoforms = len(mb_final_isoforms) + len(mp_final_isoforms)
    print(mb_mp_total_isoforms)
    print("Number of shared male brain & male pronephros isoforms")
    mb_mp_shared_isoforms = mb_final_isoforms.intersection(mp_final_isoforms)
    print(len(mb_mp_shared_isoforms))
    print("\n")
    #male brain vs testis
    print("Total number of male brain and testis, including shared isoforms")
    mb_t_total_isoforms = len(mb_final_isoforms) + len(t_final_isoforms)
    print(mb_t_total_isoforms)
    print("Number of shared male brain & testis isoforms")
    mb_t_shared_isoforms = mb_final_isoforms.intersection(t_final_isoforms)
    print(len(mb_t_shared_isoforms))
    print("\n")
    #male brain vs female pronephros
    print("Total number of male brain and female pronephros, including shared isoforms")
    mb_fp_total_isoforms = len(mb_final_isoforms) + len(fp_final_isoforms)
    print(mb_fp_total_isoforms)
    print("Number of shared male brain & female pronephros isoforms")
    mb_fp_shared_isoforms = mb_final_isoforms.intersection(fp_final_isoforms)
    print(len(mb_fp_shared_isoforms))
    print("\n")
    #male brain vs ovary
    print("Total number of male brain and ovary, including shared isoforms")
    mb_o_total_isoforms = len(mb_final_isoforms) + len(o_final_isoforms)
    print(mb_o_total_isoforms)
    print("Number of shared male brain & ovary isoforms")
    mb_o_shared_isoforms = mb_final_isoforms.intersection(o_final_isoforms)
    print(len(mb_o_shared_isoforms))
    print("\n")
    #female pronephros vs ovary
    print("Total number of female pronephros and ovary, including shared isoforms")
    fp_o_total_isoforms = len(fp_final_isoforms) + len(o_final_isoforms)
    print(fp_o_total_isoforms)
    print("Number of shared female pronephros & ovary isoforms")
    fp_o_shared_isoforms = fp_final_isoforms.intersection(o_final_isoforms)
    print(len(fp_o_shared_isoforms))
    print("\n")
    #female pronephros vs testis
    print("Total number of female pronephros and testis, including shared isoforms")
    fp_t_total_isoforms = len(fp_final_isoforms) + len(t_final_isoforms)
    print(fp_t_total_isoforms)
    print("Number of shared female pronephros & testis isoforms")
    fp_t_shared_isoforms = fp_final_isoforms.intersection(t_final_isoforms)
    print(len(fp_t_shared_isoforms))
    print("\n")
    #male pronephros vs testis
    print("Total number of male pronephros and testis, including shared isoforms")
    mp_t_total_isoforms = len(mp_final_isoforms) + len(t_final_isoforms)
    print(mp_t_total_isoforms)
    print("Number of shared male pronephros & testis isoforms")
    mp_t_shared_isoforms = mp_final_isoforms.intersection(t_final_isoforms)
    print(len(mp_t_shared_isoforms))
    print("\n")
    #male pronephros vs ovary
    print("Total number of male pronephros and ovary, including shared isoforms")
    mp_o_total_isoforms = len(mp_final_isoforms) + len(o_final_isoforms)
    print(mp_o_total_isoforms)
    print("Number of shared male pronephros & ovary isoforms")
    mp_o_shared_isoforms = mp_final_isoforms.intersection(o_final_isoforms)
    print(len(mp_o_shared_isoforms))
    print("\n")


#call functions
def call():
    gene_counts = pairwise_gene_comparison()
    isoform_counts = pairwise_isoform_comparison()

call()

#Need to examine how many genes are shared between tissues from single tissue analysis
#will use isoform counts file for input for joined tissues to remove any converted isoform ids that don't match up correctly (Gene_Isoform_Counts.py does this)
#to run script: python3 Compare_joined_tissues_genes.py <isoform counts file all female tissues> <isoform counts file all male tissues> <isoform counts file no gonads females> <isoform counts file no gonads males> <output shared all samples genes> <output shared no gonads samples genes>
#Author: Alice Naftaly, March 2020

import sys

#read in each isoform counts file
#returns a list of gene ids
def read_all_female_file():
    input_file = sys.argv[1]
    all_female_genes = []
    with open(input_file, 'r') as af_info:
        for line in af_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                all_female_genes.append(gene)
    return all_female_genes

def read_all_male_file():
    input_file = sys.argv[2]
    all_male_genes = []
    with open(input_file, 'r') as mf_info:
        for line in mf_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                all_male_genes.append(gene)
    return all_male_genes

def read_nogonads_female_file():
    input_file = sys.argv[3]
    nogonads_female_genes = []
    with open(input_file, 'r') as fng_info:
        for line in fng_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                nogonads_female_genes.append(gene)
    return nogonads_female_genes

def read_nogonads_male_file():
    input_file = sys.argv[4]
    nogonads_male_genes = []
    with open(input_file, 'r') as mng_info:
        for line in mng_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene = new_line[0]
                nogonads_male_genes.append(gene)
    return nogonads_male_genes

def compare_tissues():
    af_final_genes = set(read_all_female_file())
    am_final_genes = set(read_all_male_file())
    fng_final_genes = set(read_nogonads_female_file())
    mng_final_genes = set(read_nogonads_male_file())
    print("set intersection between the sexes for genes")
    shared_between_all_samples = af_final_genes.intersection(am_final_genes, fng_final_genes, mng_final_genes)
    print(len(shared_between_all_samples))
    print("set intersections for all somatic tissues for genes")
    shared_no_gonads = fng_final_genes.intersection(mng_final_genes)
    print(len(shared_no_gonads))
    return shared_between_all_samples, shared_no_gonads

#write gene ids to output files
def write_shared_all_samples():
    shared_between_all_samples, shared_no_gonads = compare_tissues()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for gene in shared_between_all_samples:
            final = "%s\n" % str(gene)
            out.write(final)

def write_shared_nogonads():
    shared_between_all_samples, shared_no_gonads = compare_tissues()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for gene in shared_no_gonads:
            final = "%s\n" % str(gene)
            out.write(final)

#call all functions
def call():
    shared_all_samples = write_shared_all_samples()
    shared_nogonads = write_shared_nogonads()

call()

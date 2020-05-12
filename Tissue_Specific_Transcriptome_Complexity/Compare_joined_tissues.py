#comparing isoforms between all tissue & no gonads tissue analyses analyses
#will use exon counts file for input for joined tissues to remove any converted isoform ids that don't match up correctly (Gene_Isoform_Counts.py does this)
#to run script: python3 Compare_joined_tissues.py <exon counts file all female tissues> <exon counts file all male tissues> <exon counts file no gonads females> <exon counts file no gonads males> <output shared all samples isoforms> <output shared no gonads samples isoforms>
#author: Alice Naftaly, March 2020, edited May 2020


import sys

#read in exon counts files for each single tissue
#return lists of isoform IDs
def read_all_female_file():
    input_file = sys.argv[1]
    all_female_isoforms = []
    with open(input_file, 'r') as af_info:
        for line in af_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                all_female_isoforms.append(isoform)
    return all_female_isoforms

def read_all_male_file():
    input_file = sys.argv[2]
    all_male_isoforms = []
    with open(input_file, 'r') as am_info:
        for line in am_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                all_male_isoforms.append(isoform)
    return all_male_isoforms

def read_nogonads_female_file():
    input_file = sys.argv[3]
    female_nogonads_isoforms = []
    with open(input_file, 'r') as fng_info:
        for line in fng_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                female_nogonads_isoforms.append(isoform)
    return female_nogonads_isoforms

def read_nogonads_male_file():
    input_file = sys.argv[4]
    male_nogonads_isoforms = []
    with open(input_file, 'r') as mng_info:
        for line in mng_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                male_nogonads_isoforms.append(isoform)
    return male_nogonads_isoforms

def compare_tissues():
    af_iso = set(read_all_female_file())
    am_iso = set(read_all_male_file())
    fng_iso = set(read_nogonads_female_file())
    mng_iso = set(read_nogonads_male_file())
    print("set intersection between the sexes for isoforms")
    shared_between_all_samples = af_iso.intersection(am_iso)
    print(len(shared_between_all_samples))
    print("set intersections for all somatic tissues for isoforms")
    shared_no_gonads = fng_iso.intersection(mng_iso)
    print(len(shared_no_gonads))
    return shared_between_all_samples, shared_no_gonads

#write gene ids to output files
def write_shared_all_samples():
    shared_between_all_samples, shared_no_gonads = compare_tissues()
    output = sys.argv[5]
    with open(output, 'a') as out:
        for isoform in shared_between_all_samples:
            final = "%s\n" % str(isoform)
            out.write(final)

def write_shared_nogonads():
    shared_between_all_samples, shared_no_gonads = compare_tissues()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for isoform in shared_no_gonads:
            final = "%s\n" % str(isoform)
            out.write(final)

#call all functions
def call():
    shared_all_samples = write_shared_all_samples()
    shared_nogonads = write_shared_nogonads()

call()

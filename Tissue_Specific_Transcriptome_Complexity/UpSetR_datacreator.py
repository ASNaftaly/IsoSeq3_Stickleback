#creating input file for UpSetR
#need identifiers for the first column = gene ids
#need to pull all genes from combined sexes classification file
#Then the next columns will be: Female Liver, Male Liver, Female Brain, Male Brain, Female Pronephros, Male Pronephros, Ovary, Gonad, Joined All Female tissues, Joined All Male Tissues, Joined No gonads Female Tissues, Joined No gonads Male Tissues
#each of these will have a 0 or 1 in the file depending on if that gene is in the set (present = 1, absent = 0)
#to run script: python3 UpSetR_datacreator.py <combined sexes classification file no dups> <isoform counts file female liver> <isoform counts male liver> <isoform counts female brain> <isoform counts male brain> <isoform counts female pronephros> <isoform counts male pronephros> <isoform counts ovary> <isoform counts testis> <isoform counts joined all female tissues> <isoform counts joined all male tissues> <isoform counts joined no gonads female tissues> <isoform counts joined no gonads male tissues>


import sys

#read in combined sexes classification file
#just need gene ids from this file
#returns list of gene ids
def read_combined_sexes_class():
    class_file = sys.argv[1]
    gene_list = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                gene_id = new_line[6]
                gene_list.append(gene_id)
    final_gene_list = list(set(gene_list))
    return final_gene_list


#next need to read in isoform counts files that have gene ids for each of the separate single/joined tissues
def read_fl_genes():
    fl_genes_file = sys.argv[2]
    fl_genes = []
    with open(fl_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                fl_genes.append(gene_id)
    return fl_genes

def read_ml_genes():
    ml_genes_file = sys.argv[3]
    ml_genes = []
    with open(ml_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                ml_genes.append(gene_id)
    return ml_genes

def read_fb_genes():
    fb_genes_file = sys.argv[4]
    fb_genes = []
    with open(fb_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                fb_genes.append(gene_id)
    return fb_genes

def read_mb_genes():
    mb_genes_file = sys.argv[5]
    mb_genes = []
    with open(mb_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                mb_genes.append(gene_id)
    return mb_genes

def read_fp_genes():
    fp_genes_file = sys.argv[6]
    fp_genes = []
    with open(fp_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                fp_genes.append(gene_id)
    return fp_genes

def read_mp_genes():
    mp_genes_file = sys.argv[7]
    mp_genes = []
    with open(mp_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                mp_genes.append(gene_id)
    return mp_genes

def read_o_genes():
    o_genes_file = sys.argv[8]
    o_genes = []
    with open(o_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                o_genes.append(gene_id)
    return o_genes

def read_t_genes():
    t_genes_file = sys.argv[9]
    t_genes = []
    with open(t_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                t_genes.append(gene_id)
    return t_genes

#jaf = joined all female
def read_jaf_genes():
    jaf_genes_file = sys.argv[10]
    jaf_genes = []
    with open(jaf_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                jaf_genes.append(gene_id)
    return jaf_genes

#jam = joined all male
def read_jam_genes():
    jam_genes_file = sys.argv[11]
    jam_genes = []
    with open(jam_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                jam_genes.append(gene_id)
    return jam_genes

#jngf = joined no gonads female
def read_jngf_genes():
    jngf_genes_file = sys.argv[12]
    jngf_genes = []
    with open(jngf_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                jngf_genes.append(gene_id)
    return jngf_genes

#jngm = joined no gonads male
def read_jngm_genes():
    jngm_genes_file = sys.argv[13]
    jngm_genes = []
    with open(jngm_genes_file, 'r') as gene_info:
        for line in gene_info:
            if line.startswith("Gene"):
                continue
            else:
                new_line = line.split()
                gene_id = new_line[0]
                jngm_genes.append(gene_id)
    return jngm_genes


#after reading all of the genes for the single/joined tissues, need to compare these to the combined sexes genes to create list of 0,1
#order is important here!
#returns dictionary with key = gene and value == 0 or 1 depending on if the gene is in the single/joined tissue
def convert_fl_genes():
    combined_sexes_genes = read_combined_sexes_class()
    fl_genes = read_fl_genes()
    final_fl_genes = {}
    for gene in combined_sexes_genes:
        if gene in fl_genes:
            new_value = "1"
        elif gene not in fl_genes:
            new_value = "0"
        final_fl_genes.update({gene:new_value})
    return final_fl_genes

def convert_ml_genes():
    combined_sexes_genes = read_combined_sexes_class()
    ml_genes = read_ml_genes()
    final_ml_genes = {}
    for gene in combined_sexes_genes:
        if gene in ml_genes:
            new_value = "1"
        elif gene not in ml_genes:
            new_value = "0"
        final_ml_genes.update({gene:new_value})
    return final_ml_genes

def convert_fb_genes():
    combined_sexes_genes = read_combined_sexes_class()
    fb_genes = read_fb_genes()
    final_fb_genes = {}
    for gene in combined_sexes_genes:
        if gene in fb_genes:
            new_value = "1"
        elif gene not in fb_genes:
            new_value = "0"
        final_fb_genes.update({gene:new_value})
    return final_fb_genes

def convert_mb_genes():
    combined_sexes_genes = read_combined_sexes_class()
    mb_genes = read_mb_genes()
    final_mb_genes = {}
    for gene in combined_sexes_genes:
        if gene in mb_genes:
            new_value = "1"
        elif gene not in mb_genes:
            new_value = "0"
        final_mb_genes.update({gene:new_value})
    return final_mb_genes

def convert_fp_genes():
    combined_sexes_genes = read_combined_sexes_class()
    fp_genes = read_fp_genes()
    final_fp_genes = {}
    for gene in combined_sexes_genes:
        if gene in fp_genes:
            new_value = "1"
        elif gene not in fp_genes:
            new_value = "0"
        final_fp_genes.update({gene:new_value})
    return final_fp_genes

def convert_mp_genes():
    combined_sexes_genes = read_combined_sexes_class()
    mp_genes = read_mp_genes()
    final_mp_genes = {}
    for gene in combined_sexes_genes:
        if gene in mp_genes:
            new_value = "1"
        elif gene not in mp_genes:
            new_value = "0"
        final_mp_genes.update({gene:new_value})
    return final_mp_genes


def convert_o_genes():
    combined_sexes_genes = read_combined_sexes_class()
    o_genes = read_o_genes()
    final_o_genes = {}
    for gene in combined_sexes_genes:
        if gene in o_genes:
            new_value = "1"
        elif gene not in o_genes:
            new_value = "0"
        final_o_genes.update({gene:new_value})
    return final_o_genes

def convert_t_genes():
    combined_sexes_genes = read_combined_sexes_class()
    t_genes = read_t_genes()
    final_t_genes = {}
    for gene in combined_sexes_genes:
        if gene in t_genes:
            new_value = "1"
        elif gene not in t_genes:
            new_value = "0"
        final_t_genes.update({gene:new_value})
    return final_t_genes

def convert_jaf_genes():
    combined_sexes_genes = read_combined_sexes_class()
    jaf_genes = read_jaf_genes()
    final_jaf_genes = {}
    for gene in combined_sexes_genes:
        if gene in jaf_genes:
            new_value = "1"
        elif gene not in jaf_genes:
            new_value = "0"
        final_jaf_genes.update({gene:new_value})
    return final_jaf_genes

def convert_jam_genes():
    combined_sexes_genes = read_combined_sexes_class()
    jam_genes = read_jam_genes()
    final_jam_genes = {}
    for gene in combined_sexes_genes:
        if gene in jam_genes:
            new_value = "1"
        elif gene not in jam_genes:
            new_value = "0"
        final_jam_genes.update({gene:new_value})
    return final_jam_genes

def convert_jngf_genes():
    combined_sexes_genes = read_combined_sexes_class()
    jngf_genes = read_jngf_genes()
    final_jngf_genes = {}
    for gene in combined_sexes_genes:
        if gene in jngf_genes:
            new_value = "1"
        elif gene not in jngf_genes:
            new_value = "0"
        final_jngf_genes.update({gene:new_value})
    return final_jngf_genes

def convert_jngm_genes():
    combined_sexes_genes = read_combined_sexes_class()
    jngm_genes = read_jngm_genes()
    final_jngm_genes = {}
    for gene in combined_sexes_genes:
        if gene in jngm_genes:
            new_value = "1"
        elif gene not in jngm_genes:
            new_value = "0"
        final_jngm_genes.update({gene:new_value})
    return final_jngm_genes


#now need to create final file format
def write():
    combined_sexes_genes = read_combined_sexes_class()
    fl_genes = convert_fl_genes()
    ml_genes = convert_ml_genes()
    fb_genes = convert_fb_genes()
    mb_genes = convert_mb_genes()
    fp_genes = convert_fp_genes()
    mp_genes = convert_mp_genes()
    o_genes = convert_o_genes()
    t_genes = convert_t_genes()
    jaf_genes = convert_jaf_genes()
    jam_genes = convert_jam_genes()
    jngf_genes = convert_jngf_genes()
    jngm_genes = convert_jngm_genes()
    output = sys.argv[14]
    with open(output, 'a') as out:
        header = "Identifier\tFemale.Liver\tMale.Liver\tFemale.Brain\tMale.Brain\tFemale.Pronephros\tMale.Pronephros\tOvary\tTestis\tJoined.All.Female\tJoined.All.Male\tJoined.No.Gonads.Female\tJoined.No.Gonads.Male\n"
        out.write(header)
        for gene in combined_sexes_genes:
            single_fl = fl_genes[gene]
            single_ml = ml_genes[gene]
            single_fb = fb_genes[gene]
            single_mb = mb_genes[gene]
            single_fp = fp_genes[gene]
            single_mp = mp_genes[gene]
            single_o = o_genes[gene]
            single_t = t_genes[gene]
            single_jaf = jaf_genes[gene]
            single_jam = jam_genes[gene]
            single_jngf = jngf_genes[gene]
            single_jngm = jngm_genes[gene]
            final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(gene), str(single_fl), str(single_ml), str(single_fb), str(single_mb), str(single_fp), str(single_mp), str(single_o), str(single_t), str(single_jaf), str(single_jam), str(single_jngf), str(single_jngm))
            out.write(final)

write()

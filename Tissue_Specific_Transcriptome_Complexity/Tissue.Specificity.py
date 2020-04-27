#need to identify tissue specific isoforms i.e. isoforms that are only present in one tissue (either present in both males and females or only one sex)
#will read isoform ids from exon counts files and then pull other information from the individual tissue classification file

#to run script: python3 Tissue.Specificity.py <individual tissue exon counts file, female liver> <individual tissue exon counts file, male liver> <individual tissue exon counts file, female brain> <individual tissue exon counts file, male brain> <individual tissue exon counts file, female pronephros> <individual tissue exon counts file, male pronephros> <individual tissue exon counts file, ovary> <individual tissue exon counts file, testis> <filtered classification file, female liver> <filtered classification file, male liver> <filtered classification file, female brain> <filtered classification file, male brain> <filtered classification file, female pronephros> <filtered classification file, male pronephros> <filtered classification file, ovary> <filtered classification file, testis>

import sys

#reads in exon counts files for each single tissue
#returns list of isoform ids

def read_fl_exon_counts():
    isoform_file = sys.argv[1]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_ml_exon_counts():
    isoform_file = sys.argv[2]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_fb_exon_counts():
    isoform_file = sys.argv[3]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_mb_exon_counts():
    isoform_file = sys.argv[4]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_fp_exon_counts():
    isoform_file = sys.argv[5]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_mp_exon_counts():
    isoform_file = sys.argv[6]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_o_exon_counts():
    isoform_file = sys.argv[7]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

def read_t_exon_counts():
    isoform_file = sys.argv[8]
    isoform_list = []
    with open(isoform_file, 'r') as isoforms:
        for line in isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_list.append(new_line[0])
    return isoform_list

#read in classification file
#only read in isoforms from isoform lists
#returns dictionary with key = isoform id and value == [gene id, transcript id, number of exons, splice type, structural category]
def read_fl_class():
    class_file = sys.argv[9]
    fl_isoforms = read_fl_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in fl_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_ml_class():
    class_file = sys.argv[10]
    ml_isoforms = read_ml_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in ml_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_fb_class():
    class_file = sys.argv[11]
    fb_isoforms = read_fb_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in fb_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_mb_class():
    class_file = sys.argv[12]
    mb_isoforms = read_mb_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in mb_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_fp_class():
    class_file = sys.argv[13]
    fp_isoforms = read_fp_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in fp_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_mp_class():
    class_file = sys.argv[14]
    mp_isoforms = read_mp_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in mp_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_o_class():
    class_file = sys.argv[15]
    o_isoforms = read_o_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in o_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

def read_t_class():
    class_file = sys.argv[16]
    t_isoforms = read_t_exon_counts()
    final_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                if isoform in t_isoforms:
                    exon_num = new_line[5]
                    splice_type = new_line[6]
                    gene_id = new_line[7]
                    transcript_id = new_line[8]
                    structural_type = new_line[13]
                    final = [gene_id, transcript_id, exon_num, splice_type, structural_type]
                    final_dict.update({isoform:final})
    return final_dict

#combine sexes
#will also look at gonads separately and together
#need to combine all isoforms between the sexes and remove duplicates
#returns dictionaries for liver, brain, pronephros, shared gonads, and individual gonad tissues
def combine_sexes():
    fl_isoforms = read_fl_class()
    ml_isoforms = read_ml_class()
    fb_isoforms = read_fb_class()
    mb_isoforms = read_mb_class()
    fp_isoforms = read_fp_class()
    mp_isoforms = read_mp_class()
    o_isoforms = read_o_class()
    t_isoforms = read_t_class()
    #merge liver dictionaries
    final_liver_dictionary = {}
    final_liver_dictionary.update(fl_isoforms)
    final_liver_dictionary.update(ml_isoforms)
    #merge brain dictionaries
    final_brain_dictionary = {}
    final_brain_dictionary.update(fb_isoforms)
    final_brain_dictionary.update(mb_isoforms)
    #merge pronephros dictionaries
    final_pronephros_dictionary = {}
    final_pronephros_dictionary.update(fp_isoforms)
    final_pronephros_dictionary.update(mp_isoforms)
    #shared gonad isoforms; only keeping shared because ovary and testis are different tissues
    combined_gonad_dictionary = {}
    for key in o_isoforms:
        if key in t_isoforms:
            combined_gonad_dictionary.update({key:o_isoforms[key]})
    return final_liver_dictionary, final_brain_dictionary, final_pronephros_dictionary, combined_gonad_dictionary, o_isoforms, t_isoforms


#Sort out tissue specific isoforms
#uses set difference to get unique isoforms for each tissue
#returns dictionaries for each tissue set
def sort_unique_isoforms():
    liver_dict, brain_dict, pronephros_dict, combined_gonad_dict, o_dict, t_dict = combine_sexes()
    liver_isoforms = set(liver_dict.keys())
    brain_isoforms = set(brain_dict.keys())
    pronephros_isoforms = set(pronephros_dict.keys())
    combined_gonad_isoforms = set(combined_gonad_dict.keys())
    o_isoforms = set(o_dict.keys())
    t_isoforms = set(t_dict.keys())
    final_liver_dict = {}
    final_brain_dict = {}
    final_pronephros_dict = {}
    final_combined_gonad_dict = {}
    final_o_dict = {}
    final_t_dict = {}
    #liver unique isoforms:
    unique_liver_isoforms = liver_isoforms.difference(brain_isoforms, pronephros_isoforms, combined_gonad_isoforms)
    #brain unique isoforms:
    unique_brain_isoforms = brain_isoforms.difference(liver_isoforms, pronephros_isoforms, combined_gonad_isoforms)
    #pronephros unique isoforms:
    unique_pronephros_isoforms = pronephros_isoforms.difference(liver_isoforms, brain_isoforms, combined_gonad_isoforms)
    #unique combined gonad isoforms:
    unique_combined_gonad_isoforms = combined_gonad_isoforms.difference(liver_isoforms, brain_isoforms, pronephros_isoforms)
    #unique ovary isoforms:
    unique_ovary_isoforms = o_isoforms.difference(liver_isoforms, brain_isoforms, pronephros_isoforms, t_isoforms)
    #unique testis_isoforms:
    unique_testis_isoforms = t_isoforms.difference(liver_isoforms, brain_isoforms, pronephros_isoforms, o_isoforms)
    #now create new dictionaries with only unique keys
    #final_liver dictionary
    for l_key in unique_liver_isoforms:
        single_key = liver_dict[l_key]
        final_liver_dict.update({l_key:single_key})
    #final brain dictionary
    for b_key in unique_brain_isoforms:
        single_key = brain_dict[b_key]
        final_brain_dict.update({b_key:single_key})
    #final pronephros dictionary
    for p_key in unique_pronephros_isoforms:
        single_key = pronephros_dict[p_key]
        final_pronephros_dict.update({p_key:single_key})
    #final combined gonad dictionary
    for g_key in unique_combined_gonad_isoforms:
        single_key = combined_gonad_dict[g_key]
        final_combined_gonad_dict.update({g_key:single_key})
    #final ovary dictionary
    for o_key in unique_ovary_isoforms:
        single_key = o_dict[o_key]
        final_o_dict.update({o_key:single_key})
    #final testis dictionary
    for t_key in unique_testis_isoforms:
        single_key = t_dict[t_key]
        final_t_dict.update({t_key:single_key})
    return final_liver_dict, final_brain_dict, final_pronephros_dict, final_combined_gonad_dict, final_o_dict, final_t_dict


#write output files
def write_liver():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[17]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in liver_dict:
            single_isoform = liver_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

def write_brain():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[18]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in brain_dict:
            single_isoform = brain_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

def write_pronephros():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[19]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in pronephros_dict:
            single_isoform = pronephros_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

def write_gonads():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[20]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in gonad_dict:
            single_isoform = gonad_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

def write_ovary():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[21]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in ovary_dict:
            single_isoform = ovary_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

def write_testis():
    liver_dict, brain_dict, pronephros_dict, gonad_dict, ovary_dict, testis_dict = sort_unique_isoforms()
    output = sys.argv[22]
    with open(output, 'a') as out:
        header="Isoform.ID\tGene.ID\tTranscript.ID\tExon.Number\tSplice.Type\tSub.Type\n"
        out.write(header)
        for isoform in testis_dict:
            single_isoform = testis_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(single_isoform[0]), str(single_isoform[1]), str(single_isoform[2]), str(single_isoform[3]), str(single_isoform[4]))
            out.write(final)

#call all functions
def call():
    liver = write_liver()
    brain = write_brain()
    pronephros = write_pronephros()
    gonads = write_gonads()
    ovary = write_ovary()
    testis = write_testis()

call()

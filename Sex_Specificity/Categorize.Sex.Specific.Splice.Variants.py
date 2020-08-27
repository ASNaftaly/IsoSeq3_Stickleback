#further classifying splicing variants from sex-specific isoforms by Splicing
#will read in sex specific isoforms by splicing and create dictionaries where key == gene
#note: realized some of the sex specific isoforms share genes that are only present in one sex in the sex specfic isoforms while others have the genes in the shared file; need to read in the shared isoforms too.
#will also need to read in GTF file with exons, bed files with TSSs and TTSs (combined full bed file with weird isoforms)
#to run script: python3 Categorize.Sex.Specific.Splice.Variants.py <female specific splice isoforms> <male specific splice isoforms> <shared isoforms> <isoseq gtf file> <corrected bed file> <output file 1; female specific splicing file> <output file 2;male specific splicing file>
#Author: Alice Naftaly, Aug 2020

import sys

#read female specific isoforms
#returns dictionary with key == gene and value == isoform id
def read_female_isoforms():
    female_file = sys.argv[1]
    female_dict = {}
    with open(female_file, 'r') as female_info:
        for line in female_info:
            new_line = line.split()
            isoform = new_line[1]
            gene = new_line[0]
            if gene in female_dict:
                female_dict[gene].append(isoform)
            elif gene not in female_dict:
                female_dict.update({gene:[isoform]})
    return female_dict

#read male specific isoforms
#returns dictionary with key == gene and value == isoform id
def read_male_isoforms():
    male_file = sys.argv[2]
    male_dict = {}
    with open(male_file, 'r') as male_info:
        for line in male_info:
            new_line = line.split()
            isoform = new_line[1]
            gene = new_line[0]
            if gene in male_dict:
                male_dict[gene].append(isoform)
            elif gene not in male_dict:
                male_dict.update({gene:[isoform]})
    return male_dict


#read shared isoforms
#returns dictionary with key == gene and value == isoform id
def read_shared_isoforms():
    shared_file = sys.argv[3]
    shared_dict = {}
    with open(shared_file, 'r') as shared_info:
        for line in shared_info:
            new_line = line.split()
            isoform = new_line[1]
            gene = new_line[0]
            if gene in shared_dict:
                shared_dict[gene].append(isoform)
            elif gene not in shared_dict:
                shared_dict.update({gene:[isoform]})
    return shared_dict



#read in isoseq gtf file
#returns a dictionary with key == isoform id and value == list of exons with each exon == [strand, start pos, end pos]
#minus strand isoforms are in the correct order for exon number but the start pos is the 2 value in the list for each individual exon
def read_isoseq_gtf():
    gtf_file = sys.argv[4]
    exon_dict = {}
    final_exon_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            exon_start = new_line[3]
            exon_end = new_line[4]
            strand = new_line[6]
            isoform_id_full = new_line[9].strip(";")
            isoform_id = isoform_id_full.strip("\"")
            dict_value = [strand, exon_start, exon_end]
            if isoform_id in exon_dict:
                exon_dict[isoform_id].append(dict_value)
            elif isoform_id not in exon_dict:
                exon_dict.update({isoform_id:[dict_value]})
    for isoform in exon_dict:
        single_isoform = exon_dict[isoform]
        if len(single_isoform) == 1:
            final_exon_dict.update({isoform:single_isoform})
        elif len(single_isoform) > 1:
            strand = single_isoform[0][0]
            if strand == "+":
                final_exon_dict.update({isoform:single_isoform})
            elif strand == "-":
                single_isoform.reverse()
                final_exon_dict.update({isoform:single_isoform})
    return final_exon_dict


#read in bed file
#returns dictionary with key == isoform and value == [chr, strand, start pos, end pos]
def read_bed():
    bed_file = sys.argv[5]
    bed_dict = {}
    with open(bed_file, 'r') as bed_info:
        for line in bed_info:
            new_line = line.split()
            chr = new_line[0]
            pos_1 = new_line[1]
            pos_2 = new_line[2]
            isoform = new_line[3]
            strand = new_line[5]
            if strand == "+":
                final = [chr, strand, pos_1, pos_2]
            elif strand == "-":
                final = [chr, strand, pos_2, pos_1]
            bed_dict.update({isoform:final})
    return bed_dict


#combine shared isoforms with sex_specific isoforms
#begins by adding isoforms that are shared to male and female specific isoforms
#then removes all genes that have the exact same isoforms in both sexes (this means for a given gene there were no isoform differences between all of the isoforms in either sex; had to be the same number of isoforms as well)
#returns dictionaries with key == gene and value == all isoforms
def combine_isoforms():
    shared_isoforms = read_shared_isoforms()
    female_isoforms = read_female_isoforms()
    male_isoforms = read_male_isoforms()
    for gene in shared_isoforms:
        if gene in female_isoforms:
            single_gene_shared = shared_isoforms[gene]
            for iso in single_gene_shared:
                female_isoforms[gene].append(iso)
        elif gene not in female_isoforms:
            single_gene_shared = shared_isoforms[gene]
            female_isoforms.update({gene:single_gene_shared})
        if gene in male_isoforms:
            single_gene_shared = shared_isoforms[gene]
            for iso in single_gene_shared:
                male_isoforms[gene].append(iso)
        elif gene not in male_isoforms:
            single_gene_shared = shared_isoforms[gene]
            male_isoforms.update({gene:single_gene_shared})
    exact_same_genes = []
    for g in female_isoforms:
        if g in male_isoforms:
            single_female = female_isoforms[g]
            single_male = male_isoforms[g]
            if single_female == single_male:
                exact_same_genes.append(g)
    for exact in exact_same_genes:
        del female_isoforms[exact]
        del male_isoforms[exact]
    return female_isoforms, male_isoforms

#identify sex specific isoforms that have differences in splicing
#differences in exons
#first check for differences in exon numbers, then in exon positions for isoforms with the same number of exons
#final dictionary entry ==
#isoform \t gene \t Female/Male \t Alternative.Splicing/TSS/TTS.UTR.diffs \t num exons/diff.TSS or diff.TTS \n
def characterize_differences():
    female_isoforms, male_isoforms = combine_isoforms()
    gtf_dict = read_isoseq_gtf()
    characterization_dict = {}
    for gene in female_isoforms:
        single_female_gene = female_isoforms[gene]
        single_male_gene = male_isoforms[gene]
        single_female_gene_exon_counts = {}
        single_male_gene_exon_counts = {}
        #creating dictionaries with key == exon_count and value == isoforms
        #this way, if the keys are different between males and females, then the isoforms have different exons which means alternative splicing
        for f_iso in single_female_gene:
            single_f_iso_exons = str(len(gtf_dict[f_iso]))
            if single_f_iso_exons in single_female_gene_exon_counts:
                single_female_gene_exon_counts[single_f_iso_exons].append(f_iso)
            elif single_f_iso_exons not in single_female_gene_exon_counts:
                single_female_gene_exon_counts.update({single_f_iso_exons:[f_iso]})
        for m_iso in single_male_gene:
            single_m_iso_exons = str(len(gtf_dict[m_iso]))
            if single_m_iso_exons in single_male_gene_exon_counts:
                single_male_gene_exon_counts[single_m_iso_exons].append(m_iso)
            elif single_m_iso_exons not in single_male_gene_exon_counts:
                single_male_gene_exon_counts.update({single_m_iso_exons:[m_iso]})
        #will now compare exon counts among male and female dictionaries
        for f_key in single_female_gene_exon_counts:
            #if the female count is in the male counts; this means that there is a difference in TSS or TTS/UTR regions or perhaps the exons are different but still the same number of exons
            if f_key in single_male_gene_exon_counts:
                single_female_gene_isoforms = single_female_gene_exon_counts[f_key]
                single_male_gene_isoforms = single_male_gene_exon_counts[f_key]
                if single_female_gene_isoforms == single_male_gene_isoforms:
                    #exact same isoform; can skip this as these are shared isoforms
                    continue
                elif len(single_female_gene_isoforms) == 1 and len(single_male_gene_isoforms) == 1 and single_female_gene_isoforms != single_male_gene_isoforms:
                    female_isoform = single_female_gene_isoforms[0]
                    male_isoform = single_male_gene_isoforms[0]
                    single_female_exons = gtf_dict[single_female_gene_isoforms[0]]
                    single_male_exons = gtf_dict[single_male_gene_isoforms[0]]
                    strand = single_female_exons[0][0]
                    female_exons = []
                    male_exons = []
                    for f_exon in single_female_exons:
                        f_pos_1 = f_exon[1]
                        f_pos_2 = f_exon[2]
                        if strand == '+':
                            f_exon_final = [f_pos_1, f_pos_2]
                            female_exons.append(f_exon_final)
                        elif strand == "-":
                            f_exon_final = [f_pos_2, f_pos_1]
                            female_exons.append(f_exon_final)
                    for m_exon in single_male_exons:
                        m_pos_1 = m_exon[1]
                        m_pos_2 = m_exon[2]
                        if strand == '+':
                            m_exon_final = [m_pos_1, m_pos_2]
                            male_exons.append(m_exon_final)
                        elif strand == "-":
                            m_exon_final = [m_pos_2, m_pos_1]
                            male_exons.append(m_exon_final)
                    exon_descriptors = []
                    exon_descriptor_names = []
                    for index, exon in enumerate(female_exons):
                        f_exon_start = exon[0]
                        f_exon_end = exon[1]
                        m_exon_start = male_exons[index][0]
                        m_exon_end = male_exons[index][1]
                        if f_exon_start == m_exon_start and f_exon_end == m_exon_end:
                            exon_descriptors.append([index,"same.exon"])
                            exon_descriptor_names.append("same.exon")
                        elif f_exon_start != m_exon_start and f_exon_end == m_exon_end:
                            exon_descriptors.append([index,"diff.start"])
                            exon_descriptor_names.append("diff.start")
                        elif f_exon_start == m_exon_start and f_exon_end != m_exon_end:
                            exon_descriptors.append([index, "diff.end"])
                            exon_descriptor_names.append("diff.end")
                        elif f_exon_start != m_exon_start and f_exon_end != m_exon_end:
                            exon_descriptors.append([index, "diff.exon"])
                            exon_descriptor_names.append("diff.exon")
                    if exon_descriptors == [[0, "diff.end"]]:
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing\tSingle.Exon.TTS.Different\n" % (female_isoform,gene)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing\tSingle.Exon.TTS.Different\n" % (male_isoform,gene)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors == [[0, "diff.exon"]]:
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.Different\n" % (female_isoform,gene)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.Different\n" % (male_isoform,gene)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif len(exon_descriptors) == 2 and exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[1] == [1, "diff.end"]:
                        f_final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (female_isoform,gene)
                        m_final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (male_isoform,gene)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif len(exon_descriptors) == 2 and exon_descriptors[0] == [0, "diff.end"] and exon_descriptors[1] == [1, "diff.exon"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[1:len(exon_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif len(exon_descriptors) == 2 and exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[1] == [1, "diff.exon"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif len(exon_descriptors) == 2 and exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[1] == [1, "diff.end"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"] and list(set(exon_descriptor_names[1:len(exon_descriptors)-1])) == ["same.exon"]:
                        f_final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (female_isoform,gene)
                        m_final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (male_isoform,gene)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"] and len(list(set(exon_descriptor_names[1:len(exon_descriptors)-1]))) >1:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.exon"]:
                         descriptor_count = 0
                         for descriptor in exon_descriptor_names[1:len(exon_descriptors)]:
                             if descriptor != "same.exon":
                                 descriptor_count += 1
                         f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                         m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                         if female_isoform in characterization_dict:
                             characterization_dict[female_isoform].append(f_final)
                         elif female_isoform not in characterization_dict:
                             characterization_dict.update({female_isoform:[f_final]})
                         if male_isoform in characterization_dict:
                             characterization_dict[male_isoform].append(m_final)
                         elif male_isoform not in characterization_dict:
                             characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.start"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[1:len(exon_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "same.exon"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[1:len(exon_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.start"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[1:len(exon_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.end"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.start"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.end"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[0:len(exon_descriptors)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced;Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[0:len(exon_descriptors)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.start"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.exon"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "same.exon"] and list(set(exon_descriptor_names[1:len(exon_descriptor_names)])) == ["same.exon"]:
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,"1")
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,"1")
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "diff.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "same.exon"] and len(list(set(exon_descriptor_names[1:len(exon_descriptor_names)]))) > 1:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "same.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"] and list(set(exon_descriptor_names[0:len(exon_descriptor_names)-1])) == ["same.exon"]:
                        f_final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (female_isoform,gene)
                        m_final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (male_isoform,gene)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "same.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.end"] and len(list(set(exon_descriptor_names[0:len(exon_descriptor_names)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names[0:len(exon_descriptors)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "same.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "diff.start"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                    elif exon_descriptors[0] == [0, "same.exon"] and exon_descriptors[len(exon_descriptors)-1] == [int(f_key)-1, "same.exon"]:
                        descriptor_count = 0
                        for descriptor in exon_descriptor_names:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        f_final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (female_isoform,gene,descriptor_count)
                        m_final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (male_isoform,gene,descriptor_count)
                        if female_isoform in characterization_dict:
                            characterization_dict[female_isoform].append(f_final)
                        elif female_isoform not in characterization_dict:
                            characterization_dict.update({female_isoform:[f_final]})
                        if male_isoform in characterization_dict:
                            characterization_dict[male_isoform].append(m_final)
                        elif male_isoform not in characterization_dict:
                            characterization_dict.update({male_isoform:[m_final]})
                elif len(single_female_gene_isoforms) == 1 and single_female_gene_isoforms[0] in single_male_gene_isoforms and len(single_male_gene_isoforms) == 2:
                    m_single_exons = []
                    male_exon_numbers = []
                    for m_single in single_male_gene_isoforms:
                        m_single_exons += gtf_dict[m_single]
                        male_exon_numbers.append(len(gtf_dict[m_single]))
                    #goes through exons and looks for differences between two isoforms
                    non_matching_exons = []
                    non_matching_exons_descriptors = []
                    x = 0
                    while x < int(f_key):
                        exon_1 = m_single_exons[x]
                        exon_2 = m_single_exons[x + int(f_key)]
                        exon_1_strand = exon_1[0]
                        if exon_1_strand == "+":
                            exon_1_start = int(exon_1[1])
                            exon_1_end = int(exon_1[2])
                            exon_2_start = int(exon_2[1])
                            exon_2_end = int(exon_2[2])
                        elif exon_1_strand == "-":
                            exon_1_start = int(exon_1[2])
                            exon_1_end = int(exon_1[1])
                            exon_2_start = int(exon_2[2])
                            exon_2_end = int(exon_2[1])
                        #if exon completely matches
                        if exon_1_start == exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "same.exon"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("same.exon")
                            x += 1
                        #if the exon starts do not equal but the ends do:
                        elif exon_1_start != exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "start.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("start.diff")
                            x += 1
                        #if the exon ends are different and the starts are the same
                        elif exon_1_start == exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "end.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("end.diff")
                            x += 1
                        #if both the ends and the starts are different
                        elif exon_1_start != exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "both.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("both.diff")
                            x += 1
                    #if everything is the same except the first exon and the first exon is different in both start and end position, this is alterantive splicing and alternative TSS
                    if non_matching_exons == [[0, "both.diff"]]:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing\tSingle.Exon.Both.Start.Stop.Different\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1])) == ['same.exon'] and [int(f_key)-1, "end.diff"]:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0,"start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons and len(non_matching_exons) == 2:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif len(non_matching_exons) == 2 and [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif len(non_matching_exons) == 2 and [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons and len(non_matching_exons) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif non_matching_exons == [[0, "end.diff"]]:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tSingle.Exon.Different.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "start.diff"] in non_matching_exons:
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(2))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_male_gene_isoforms:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                elif len(single_female_gene_isoforms) == 2 and single_male_gene_isoforms[0] in single_male_gene_isoforms and len(single_male_gene_isoforms) == 1:
                    f_single_exons = []
                    female_exon_numbers = []
                    for f_single in single_female_gene_isoforms:
                        f_single_exons += gtf_dict[f_single]
                        female_exon_numbers.append(len(gtf_dict[f_single]))
                    #goes through exons and looks for differences between two isoforms
                    non_matching_exons = []
                    non_matching_exons_descriptors = []
                    x = 0
                    while x < int(f_key):
                        exon_1 = f_single_exons[x]
                        exon_2 = f_single_exons[x + int(f_key)]
                        exon_1_strand = exon_1[0]
                        if exon_1_strand == "+":
                            exon_1_start = int(exon_1[1])
                            exon_1_end = int(exon_1[2])
                            exon_2_start = int(exon_2[1])
                            exon_2_end = int(exon_2[2])
                        elif exon_1_strand == "-":
                            exon_1_start = int(exon_1[2])
                            exon_1_end = int(exon_1[1])
                            exon_2_start = int(exon_2[2])
                            exon_2_end = int(exon_2[1])
                        #if exon completely matches
                        if exon_1_start == exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "same.exon"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("same.exon")
                            x += 1
                        #if the exon starts do not equal but the ends do:
                        elif exon_1_start != exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "start.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("start.diff")
                            x += 1
                        #if the exon ends are different and the starts are the same
                        elif exon_1_start == exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "end.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("end.diff")
                            x += 1
                        #if both the ends and the starts are different
                        elif exon_1_start != exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "both.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("both.diff")
                            x += 1
                    #if everything is the same except the first exon and the first exon is different in both start and end position, this is alterantive splicing and alternative TSS
                    if non_matching_exons == [[0, "both.diff"]]:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing\tSingle.Exon.Both.Start.Stop.Different\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1])) == ['same.exon'] and [int(f_key)-1, "end.diff"]:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0,"start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons and female_exon_numbers[0] == 2:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif female_exon_numbers[0] == 2 and [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif non_matching_exons == [[0, "end.diff"]]:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tSingle.Exon.Different.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons and list(set(non_matching_exons_descriptors[1:len(non_matching_exons_descriptors)])) == ["same.exon"]:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(1))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "start.diff"] in non_matching_exons:
                        for isoform in single_female_gene_isoforms:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(2))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                else:
                    all_exons = []
                    all_isoforms = []
                    all_isoforms_dict = {}
                    for f_iso in single_female_gene_isoforms:
                        f_exons = gtf_dict[f_iso]
                        strand = f_exons[0][0]
                        for exon in f_exons:
                            if strand == "+":
                                exon_start = exon[1]
                                exon_end = exon[2]
                            elif strand == "-":
                                exon_start = exon[2]
                                exon_end = exon[1]
                            all_exons.append(exon_start)
                            all_exons.append(exon_end)
                        all_isoforms.append(f_iso)
                    for m_iso in single_male_gene_isoforms:
                        m_exons = gtf_dict[m_iso]
                        strand = m_exons[0][0]
                        for exon in m_exons:
                            if strand == "+":
                                exon_start = exon[1]
                                exon_end = exon[2]
                            elif strand == "-":
                                exon_start = exon[2]
                                exon_end = exon[1]
                            all_exons.append(exon_start)
                            all_exons.append(exon_end)
                        all_isoforms.append(m_iso)
                    x = 0
                    y = 0
                    while x < int(f_key)*2:
                        single_exon_start = all_exons[x::int(f_key)*2]
                        single_exon_end = all_exons[x+1::int(f_key)*2]
                        set_exon_starts = list(set(single_exon_start))
                        set_exon_ends = list(set(single_exon_end))
                        if len(set_exon_starts) == 1 and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("same.exon")
                                a += 1
                            all_isoforms_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.start")
                                a += 1
                            all_isoforms_dict.update({y:value})
                        elif len(set_exon_starts) == 1 and len(set_exon_ends) == len(single_exon_end):
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.end")
                                a += 1
                            all_isoforms_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == len(single_exon_end):
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.exon")
                                a += 1
                            all_isoforms_dict.update({y:value})
                        elif len(set_exon_starts) == 1 and 1 <= len(set_exon_ends) <= len(single_exon_end):
                            value = []
                            for v in set_exon_ends:
                                count = single_exon_end.count(v)
                                if count == 1:
                                    start_index = single_exon_end.index(v)
                                    single_value = [start_index, "diff.end"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_end):
                                    start_indeces = [i for i, x in enumerate(single_exon_end) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                                else:
                                    start_indeces = [i for i, x in enumerate(single_exon_end) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "diff.start"]
                                        value.append(single_value)
                            all_isoforms_dict.update({y:value})
                        elif 1 <= len(set_exon_starts) <= len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            for v in set_exon_starts:
                                count = single_exon_start.count(v)
                                if count == 1:
                                    start_index = single_exon_start.index(v)
                                    single_value = [start_index, "diff.start"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_start):
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                                else:
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "diff.start"]
                                        value.append(single_value)
                            all_isoforms_dict.update({y:value})
                        elif 1 <= len(set_exon_starts) <= len(single_exon_start) and 1 <= len(set_exon_ends) <= len(single_exon_end):
                            value = []
                            dup_start_list = []
                            dup_start_indeces = []
                            dup_stop_list = []
                            for v in set_exon_starts:
                                start_count = single_exon_start.count(v)
                                if start_count == 1:
                                    start_index = single_exon_start.index(v)
                                    stop_pos = single_exon_end[start_index]
                                    stop_count = single_exon_end.count(stop_pos)
                                    if stop_count == 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.exon"]
                                        value.append(single_value)
                                    elif stop_count > 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.start"]
                                        value.append(single_value)
                                    all_isoforms_dict.update({y:value})
                                elif start_count > 1:
                                    dup_start_list = [a for i,a in enumerate(single_exon_start) if a == v]
                                    dup_start_indeces = [i for i,a in enumerate(single_exon_start) if a == v]
                                    dup_stop_list = [single_exon_end[i] for i in dup_start_indeces]
                                    if len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) == 1:
                                        for index in dup_start_indeces:
                                            single_value = [index, "same.exon"]
                                            value.append(single_value)
                                    elif len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) > 1:
                                        for stop in list(set(dup_stop_list)):
                                            stop_count = single_exon_end.count(stop)
                                            if stop_count == 1:
                                                stop_index = single_exon_end.index(stop)
                                                single_value = [stop_index, "diff.end"]
                                                value.append(single_value)
                                            elif stop_count > 1:
                                                dup_end_list = [a for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_end_indeces = [i for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_begin_list = [single_exon_start[i] for i in dup_end_indeces]
                                                if len(set(dup_end_list)) == 1 and len(set(dup_begin_list)) == 1:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "same.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                                else:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "diff.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                    all_isoforms_dict.update({y:value})
                        x += 2
                        y += 1
                    converted_all_isoforms_dict = {}
                    for k in all_isoforms_dict:
                        single_k = all_isoforms_dict[k]
                        if len(single_k) == len(all_isoforms) and isinstance(single_k[0], list) == False:
                            if k in converted_all_isoforms_dict:
                                converted_all_isoforms_dict[k].append(single_k)
                            elif k not in converted_all_isoforms_dict:
                                converted_all_isoforms_dict.update({k:[single_k]})
                        else:
                            sorted_k = sorted(single_k, key=lambda x: x[0])
                            if len(sorted_k) == len(all_isoforms):
                                final_value = []
                                for val in sorted_k:
                                    final_value.append(val[1])
                                if k in converted_all_isoforms_dict:
                                    converted_all_isoforms_dict[k].append(final_values)
                                elif k not in converted_all_isoforms_dict:
                                    converted_all_isoforms_dict.update({k:[final_value]})
                            else:
                                num_isoforms = len(all_isoforms)
                                exons = [c[0] for c in sorted_k]
                                exon_values = [c[1] for c in sorted_k]
                                set_exons = list(set(exons))
                                c_more_than_1 = []
                                for c in set_exons:
                                    c_count = exons.count(c)
                                    if c_count > 1:
                                        c_more_than_1.append(c)
                                if len(c_more_than_1) == 1:
                                    exon_index = [i for i,a in enumerate(exons) if a == c_more_than_1[0]]
                                    sorted_k_copy = [a[1] for i,a in enumerate(sorted_k) if i not in exon_index]
                                    dups = []
                                    final_index = exon_index[0]
                                    for ind in exon_index:
                                        dups.append(exon_values[ind])
                                    if "diff.exon" in dups:
                                        final_exon_value = "diff.exon"
                                    final_value = sorted_k_copy
                                    final_value.insert(final_index, final_exon_value)
                                    if k in converted_all_isoforms_dict:
                                        converted_all_isoforms_dict[k].append(final_value)
                                    elif k not in converted_all_isoforms_dict:
                                        converted_all_isoforms_dict.update({k:[final_value]})
                                elif len(c_more_than_1) > 1:
                                    exon_index = []
                                    for val in c_more_than_1:
                                        for i, a in enumerate(exons):
                                            if a == val:
                                                exon_index.append(i)
                                    sorted_k_copy = [a[1] for i,a in enumerate(sorted_k) if i not in exon_index]
                                    dups = []
                                    dup_indeces = []
                                    for v in exon_index:
                                        dups.append(exon_values[v])
                                        dup_indeces.append(exons[v])
                                    final_indeces = dup_indeces[0::2]
                                    h = 0
                                    g = 0
                                    while h < len(dups):
                                        single_exon_dups = dups[h:h+2]
                                        if "diff.exon" in single_exon_dups:
                                            final_exon_value = "diff.exon"
                                            final_value = sorted_k_copy
                                            final_value.insert(final_indeces[g], final_exon_value)
                                        h += 2
                                        g += 1
                                    if k in converted_all_isoforms_dict:
                                        converted_all_isoforms_dict[k].append(final_value)
                                    elif k not in converted_all_isoforms_dict:
                                        converted_all_isoforms_dict.update({k:[final_value]})
                    for index, isoform in enumerate(all_isoforms):
                        single_isoform_values = []
                        if isoform in single_female_gene_isoforms:
                            final_isoform = isoform
                            sex = "Female"
                        elif isoform in single_male_gene_isoforms:
                            final_isoform = isoform
                            sex = "Male"
                        for exon in converted_all_isoforms_dict:
                            single_exon = converted_all_isoforms_dict[exon][0]
                            single_isoform_values.append(single_exon[index])
                        if len(single_isoform_values) == 1 and single_isoform_values == ["same.exon"]:
                            continue
                        elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.end"]:
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.Alternate.TTS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\tSingle.Exon.Both.TSS.TTS.different\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "diff.end":
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "diff.start":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "diff.exon":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "same.exon":
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.end" and single_isoform_values[1] == "diff.end":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.end":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.start":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(2))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.exon":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(2))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "same.exon" and single_isoform_values[1] == "same.exon":
                            continue
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "same.exon" and single_isoform_values[1] == "diff.start":
                            final = "%s\t%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "same.exon" and single_isoform_values[1] == "diff.end":
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif len(single_isoform_values) == 2 and single_isoform_values[0] == "same.exon" and single_isoform_values[1] == "diff.exon":
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(1))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif list(set(single_isoform_values)) == ["same.exon"]:
                            continue
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[1:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[0:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (isoform,gene,sex)
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[1:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon" and len(list(set(single_isoform_values))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[0:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[0:len(single_isoform_values)-1]))[0].startswith("diff"):
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene,sex, str(descriptor_count))
                            if final_isoform in characterization_dict:
                                characterization_dict[final_isoform].append(final)
                            elif final_isoform not in characterization_dict:
                                characterization_dict.update({final_isoform:[final]})
            #if the female key is not in the male dictionary, this means these isoforms likely represent alternative splicing events
            #however if there is more than one isoform per key, need to figure out the differences between the isoforms
            elif f_key not in single_male_gene_exon_counts:
                single_f_counts = single_female_gene_exon_counts[f_key]
                #only 1 isoform == alternative splicing
                #cannot easily get info as to which exons are shared with other isoforms
                if len(single_f_counts) == 1:
                    final = "%s\t%s\tFemale\tAlternative.Splicing\tNo.Other.Isoforms.From.This.Gene.with.this.number.of.exons\n" % (single_f_counts[0],gene)
                    if single_f_counts[0] in characterization_dict:
                        characterization_dict[single_f_counts[0]].append(final)
                    elif single_f_counts[0] not in characterization_dict:
                        characterization_dict.update({single_f_counts[0]:[final]})
                elif len(single_f_counts) == 2:
                    f_single_exons = []
                    female_exon_numbers = []
                    for f_single in single_f_counts:
                        f_single_exons += gtf_dict[f_single]
                        female_exon_numbers.append(len(gtf_dict[f_single]))
                    #goes through exons and looks for differences between two isoforms
                    non_matching_exons = []
                    non_matching_exons_descriptors = []
                    x = 0
                    while x < int(f_key):
                        exon_1 = f_single_exons[x]
                        exon_2 = f_single_exons[x + int(f_key)]
                        exon_1_strand = exon_1[0]
                        if exon_1_strand == "+":
                            exon_1_start = int(exon_1[1])
                            exon_1_end = int(exon_1[2])
                            exon_2_start = int(exon_2[1])
                            exon_2_end = int(exon_2[2])
                        elif exon_1_strand == "-":
                            exon_1_start = int(exon_1[2])
                            exon_1_end = int(exon_1[1])
                            exon_2_start = int(exon_2[2])
                            exon_2_end = int(exon_2[1])
                        #if exon completely matches
                        if exon_1_start == exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "same.exon"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("same.exon")
                            x += 1
                        #if the exon starts do not equal but the ends do:
                        elif exon_1_start != exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "start.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("start.diff")
                            x += 1
                        #if the exon ends are different and the starts are the same
                        elif exon_1_start == exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "end.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("end.diff")
                            x += 1
                        #if both the ends and the starts are different
                        elif exon_1_start != exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "both.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("both.diff")
                            x += 1
                    #if everything is the same except the first exon and the first exon is different in both start and end position, this is alterantive splicing and alternative TSS
                    if non_matching_exons == [[0, "both.diff"]]:
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing\tSingle.Exon.Both.Start.Stop.Different\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1])) == ['same.exon'] and [int(f_key)-1, "end.diff"]:
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0,"start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons and female_exon_numbers[0] == 2:
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif female_exon_numbers[0] == 2 and [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif female_exon_numbers[0] == 2 and [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons:
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(f_key)-1, "same.exon"] in non_matching_exons and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_f_counts:
                            final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                else:
                    f_single_exons = []
                    f_single_exons_dict = {}
                    female_exon_numbers = []
                    for f_single in single_f_counts:
                        f_single_exon = gtf_dict[f_single]
                        female_exon_numbers.append(len(f_single_exon))
                        for f_exon in f_single_exon:
                            strand = f_exon[0]
                            if strand == "+":
                                f_exon_start = f_exon[1]
                                f_exon_end = f_exon[2]
                            elif strand == "-":
                                f_exon_start = f_exon[2]
                                f_exon_end = f_exon[1]
                            final_exon = [f_exon_start, f_exon_end]
                            f_single_exons += final_exon
                    x = 0
                    y = 0
                    while x < int(f_key)*2:
                        single_exon_start = f_single_exons[x::int(f_key)*2]
                        single_exon_end = f_single_exons[x+1::int(f_key)*2]
                        set_exon_starts = list(set(single_exon_start))
                        set_exon_ends = list(set(single_exon_end))
                        if len(set_exon_starts) == 1 and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("same.exon")
                                a += 1
                            f_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.start")
                                a += 1
                            f_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == len(single_exon_end):
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.exon")
                                a += 1
                            f_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == 1 and 1 < len(set_exon_ends) < len(single_exon_end):
                            value = []
                            for v in set_exon_ends:
                                count = single_exon_end.count(v)
                                if count == 1:
                                    start_index = single_exon_end.index(v)
                                    single_value = [start_index, "diff.end"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_end):
                                    start_indeces =  [i for i, x in enumerate(single_exon_end) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                            f_single_exons_dict.update({y:value})
                        elif 1 < len(set_exon_starts) < len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            for v in set_exon_starts:
                                count = single_exon_start.count(v)
                                if count == 1:
                                    start_index = single_exon_start.index(v)
                                    single_value = [start_index, "diff.start"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_start):
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                                else:
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "diff.start"]
                                        value.append(single_value)
                            f_single_exons_dict.update({y:value})
                        elif 1 <= len(set_exon_starts) <= len(single_exon_start) and 1 <= len(set_exon_ends) <= len(single_exon_end):
                            value = []
                            dup_start_list = []
                            dup_start_indeces = []
                            dup_stop_list = []
                            for v in set_exon_starts:
                                start_count = single_exon_start.count(v)
                                if start_count == 1:
                                    start_index = single_exon_start.index(v)
                                    stop_pos = single_exon_end[start_index]
                                    stop_count = single_exon_end.count(stop_pos)
                                    if stop_count == 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.exon"]
                                        value.append(single_value)
                                    elif stop_count > 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.start"]
                                        value.append(single_value)
                                    f_single_exons_dict.update({y:value})
                                elif start_count > 1:
                                    dup_start_list = [a for i,a in enumerate(single_exon_start) if a == v]
                                    dup_start_indeces = [i for i,a in enumerate(single_exon_start) if a == v]
                                    dup_stop_list = [single_exon_end[i] for i in dup_start_indeces]
                                    if len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) == 1:
                                        for index in dup_start_indeces:
                                            single_value = [index, "same.exon"]
                                            value.append(single_value)
                                    elif len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) > 1:
                                        for stop in list(set(dup_stop_list)):
                                            stop_count = single_exon_end.count(stop)
                                            if stop_count == 1:
                                                stop_index = single_exon_end.index(stop)
                                                single_value = [stop_index, "diff.end"]
                                                value.append(single_value)
                                            elif stop_count > 1:
                                                dup_end_list = [a for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_end_indeces = [i for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_begin_list = [single_exon_start[i] for i in dup_end_indeces]
                                                if len(set(dup_end_list)) == 1 and len(set(dup_begin_list)) == 1:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "same.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                                else:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "diff.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                    f_single_exons_dict.update({y:value})
                        y += 1
                        x += 2
                    converted_f_dict_exons = {}
                    for k in f_single_exons_dict:
                        single_k = f_single_exons_dict[k]
                        if len(single_k) == len(single_f_counts) and isinstance(single_k[0], list) == False:
                            if k in converted_f_dict_exons:
                                converted_f_dict_exons[k].append(single_k)
                            elif k not in converted_f_dict_exons:
                                converted_f_dict_exons.update({k:[single_k]})
                        else:
                            sorted_k = sorted(single_k, key=lambda x: x[0])
                            if len(sorted_k) == len(single_f_counts):
                                final_value = []
                                for val in sorted_k:
                                    final_value.append(val[1])
                                if k in converted_f_dict_exons:
                                    converted_f_dict_exons[k].append(final_values)
                                elif k not in converted_f_dict_exons:
                                    converted_f_dict_exons.update({k:[final_value]})
                            else:
                                exons = [c[0] for c in sorted_k]
                                exon_values = [c[1] for c in sorted_k]
                                set_exons = list(set(exons))
                                c_more_than_1 = []
                                for c in set_exons:
                                    c_count = exons.count(c)
                                    if c_count > 1:
                                        c_more_than_1.append(c)
                                if len(c_more_than_1) == 1:
                                    if c_more_than_1[0] == 0:
                                        final_value = ["diff.exon", exon_values[1], exon_values[2]]
                                        if k in converted_f_dict_exons:
                                            converted_f_dict_exons[k].append(final_value)
                                        elif k not in converted_f_dict_exons:
                                            converted_f_dict_exons.update({k:[final_value]})
                                    elif c_more_than_1[0] == 1:
                                        final_value = [exon_values[0], "diff.exon", exon_values[2]]
                                        if k in converted_f_dict_exons:
                                            converted_f_dict_exons[k].append(final_value)
                                        elif k not in converted_f_dict_exons:
                                            converted_f_dict_exons.update({k:[final_value]})
                                    elif c_more_than_1[0] == 2:
                                        final_value = [exon_values[0], exon_values[1], "diff.exon"]
                                        if k in converted_f_dict_exons:
                                            converted_f_dict_exons[k].append(final_value)
                                        elif k not in converted_f_dict_exons:
                                            converted_f_dict_exons.update({k:[final_value]})
                                elif len(c_more_than_1) > 1:
                                    final_value = [exon_values[0], exon_values[1],"diff.exon", "diff.exon", exon_values[4]]
                                    if k in converted_f_dict_exons:
                                        converted_f_dict_exons[k].append(final_value)
                                    elif k not in converted_f_dict_exons:
                                        converted_f_dict_exons.update({k:[final_value]})
                                else:
                                    final_value = ["same.exon", exon_values[0], exon_values[1],"same.exon",  exon_values[2]]
                                    if k in converted_f_dict_exons:
                                        converted_f_dict_exons[k].append(final_value)
                                    elif k not in converted_f_dict_exons:
                                        converted_f_dict_exons.update({k:[final_value]})
                    for index, isoform in enumerate(single_f_counts):
                        single_isoform_values = []
                        for exon in converted_f_dict_exons:
                            single_exon = converted_f_dict_exons[exon][0]
                            single_isoform_values.append(single_exon[index])
                        if single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[1:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[1:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[0:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\tFemale\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[0:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tFemale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
        for m_key in single_male_gene_exon_counts:
            if m_key not in single_female_gene_exon_counts:
                single_m_counts = single_male_gene_exon_counts[m_key]
                #only 1 isoform == alternative splicing
                #cannot easily get info as to which exons are shared with other isoforms
                if len(single_m_counts) == 1:
                    final = "%s\t%s\tMale\tAlternative.Splicing\tNo.Other.Isoforms.From.This.Gene.with.this.number.of.exons\n" % (single_m_counts[0],gene)
                    if single_m_counts[0] in characterization_dict:
                        characterization_dict[single_m_counts[0]].append(final)
                    elif single_m_counts[0] not in characterization_dict:
                        characterization_dict.update({single_m_counts[0]:[final]})
                elif len(single_m_counts) == 2:
                    m_single_exons = []
                    male_exon_numbers = []
                    for m_single in single_m_counts:
                        m_single_exons += gtf_dict[m_single]
                        male_exon_numbers.append(len(gtf_dict[m_single]))
                    #goes through exons and looks for differences between two isoforms
                    non_matching_exons = []
                    non_matching_exons_descriptors = []
                    x = 0
                    while x < int(m_key):
                        exon_1 = m_single_exons[x]
                        exon_2 = m_single_exons[x + int(m_key)]
                        exon_1_strand = exon_1[0]
                        if exon_1_strand == "+":
                            exon_1_start = int(exon_1[1])
                            exon_1_end = int(exon_1[2])
                            exon_2_start = int(exon_2[1])
                            exon_2_end = int(exon_2[2])
                        elif exon_1_strand == "-":
                            exon_1_start = int(exon_1[2])
                            exon_1_end = int(exon_1[1])
                            exon_2_start = int(exon_2[2])
                            exon_2_end = int(exon_2[1])
                        #if exon completely matches
                        if exon_1_start == exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "same.exon"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("same.exon")
                            x += 1
                        #if the exon starts do not equal but the ends do:
                        elif exon_1_start != exon_2_start and exon_1_end == exon_2_end:
                            value = [x, "start.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("start.diff")
                            x += 1
                        #if the exon ends are different and the starts are the same
                        elif exon_1_start == exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "end.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("end.diff")
                            x += 1
                        #if both the ends and the starts are different
                        elif exon_1_start != exon_2_start and exon_1_end != exon_2_end:
                            value = [x, "both.diff"]
                            non_matching_exons.append(value)
                            non_matching_exons_descriptors.append("both.diff")
                            x += 1
                    #if everything is the same except the first exon and the first exon is different in both start and end position, this is alterantive splicing and alternative TSS
                    if non_matching_exons == [[0, "both.diff"]]:
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing\tSingle.Exon.Both.Start.Stop.Different\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1])) == ['same.exon'] and [int(m_key)-1, "end.diff"]:
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0,"start.diff"] in non_matching_exons and [int(m_key)-1, "end.diff"] in non_matching_exons and male_exon_numbers[0] == 2:
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif male_exon_numbers[0] == 2 and [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "both.diff"] in non_matching_exons:
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif male_exon_numbers[0] == 2 and [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "same.exon"] in non_matching_exons:
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\tBoth.Exons.Different,Alternate.TSS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(m_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "end.diff"] in non_matching_exons and [int(m_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "end.diff"] and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(m_key)-1, "end.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[0:len(non_matching_exons)-1]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(m_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(m_key)-1, "same.exon"] in non_matching_exons and len(list(set(non_matching_exons_descriptors[1:len(non_matching_exons)-1]))) > 1:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(m_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "start.diff"] in non_matching_exons and [int(m_key)-1, "both.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors[1:len(non_matching_exons_descriptors)]:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "same.exon"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "both.diff"] in non_matching_exons and [int(m_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                    elif [0, "same.exon"] in non_matching_exons and [int(m_key)-1, "start.diff"] in non_matching_exons:
                        descriptor_count = 0
                        for descriptor in non_matching_exons_descriptors:
                            if descriptor != "same.exon":
                                descriptor_count += 1
                        for isoform in single_m_counts:
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                else:
                    m_single_exons = []
                    m_single_exons_dict = {}
                    male_exon_numbers = []
                    for m_single in single_m_counts:
                        m_single_exon = gtf_dict[m_single]
                        male_exon_numbers.append(len(m_single_exon))
                        for m_exon in m_single_exon:
                            strand = m_exon[0]
                            if strand == "+":
                                m_exon_start = m_exon[1]
                                m_exon_end = m_exon[2]
                            elif strand == "-":
                                m_exon_start = m_exon[2]
                                m_exon_end = m_exon[1]
                            final_exon = [m_exon_start, m_exon_end]
                            m_single_exons += final_exon
                    x = 0
                    y = 0
                    while x < int(m_key)*2:
                        single_exon_start = m_single_exons[x::int(m_key)*2]
                        single_exon_end = m_single_exons[x+1::int(m_key)*2]
                        set_exon_starts = list(set(single_exon_start))
                        set_exon_ends = list(set(single_exon_end))
                        if len(set_exon_starts) == 1 and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("same.exon")
                                a += 1
                            m_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.start")
                                a += 1
                            m_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == len(single_exon_start) and len(set_exon_ends) == len(single_exon_end):
                            value = []
                            a = 0
                            while a < len(single_exon_start):
                                value.append("diff.exon")
                                a += 1
                            m_single_exons_dict.update({y:value})
                        elif len(set_exon_starts) == 1 and 1 < len(set_exon_ends) < len(single_exon_end):
                            value = []
                            for v in set_exon_ends:
                                count = single_exon_end.count(v)
                                if count == 1:
                                    start_index = single_exon_end.index(v)
                                    single_value = [start_index, "diff.end"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_end):
                                    start_indeces =  [i for i, x in enumerate(single_exon_end) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                            m_single_exons_dict.update({y:value})
                        elif 1 < len(set_exon_starts) < len(single_exon_start) and len(set_exon_ends) == 1:
                            value = []
                            for v in set_exon_starts:
                                count = single_exon_start.count(v)
                                if count == 1:
                                    start_index = single_exon_start.index(v)
                                    single_value = [start_index, "diff.start"]
                                    value.append(single_value)
                                elif count + 1 == len(single_exon_start):
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "same.exon"]
                                        value.append(single_value)
                                else:
                                    start_indeces =  [i for i, x in enumerate(single_exon_start) if x == v]
                                    for ind in start_indeces:
                                        single_value = [ind, "diff.start"]
                                        value.append(single_value)
                            m_single_exons_dict.update({y:value})
                        elif 1 <= len(set_exon_starts) <= len(single_exon_start) and 1 <= len(set_exon_ends) <= len(single_exon_end):
                            value = []
                            dup_start_list = []
                            dup_start_indeces = []
                            dup_stop_list = []
                            for v in set_exon_starts:
                                start_count = single_exon_start.count(v)
                                if start_count == 1:
                                    start_index = single_exon_start.index(v)
                                    stop_pos = single_exon_end[start_index]
                                    stop_count = single_exon_end.count(stop_pos)
                                    if stop_count == 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.exon"]
                                        value.append(single_value)
                                    elif stop_count > 1:
                                        start_index = single_exon_start.index(v)
                                        single_value = [start_index, "diff.start"]
                                        value.append(single_value)
                                    m_single_exons_dict.update({y:value})
                                elif start_count > 1:
                                    dup_start_list = [a for i,a in enumerate(single_exon_start) if a == v]
                                    dup_start_indeces = [i for i,a in enumerate(single_exon_start) if a == v]
                                    dup_stop_list = [single_exon_end[i] for i in dup_start_indeces]
                                    if len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) == 1:
                                        for index in dup_start_indeces:
                                            single_value = [index, "same.exon"]
                                            value.append(single_value)
                                    elif len(set(dup_start_list)) == 1 and len(set(dup_stop_list)) > 1:
                                        for stop in list(set(dup_stop_list)):
                                            stop_count = single_exon_end.count(stop)
                                            if stop_count == 1:
                                                stop_index = single_exon_end.index(stop)
                                                single_value = [stop_index, "diff.end"]
                                                value.append(single_value)
                                            elif stop_count > 1:
                                                dup_end_list = [a for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_end_indeces = [i for i,a in enumerate(single_exon_end) if a == stop]
                                                dup_begin_list = [single_exon_start[i] for i in dup_end_indeces]
                                                if len(set(dup_end_list)) == 1 and len(set(dup_begin_list)) == 1:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "same.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                                else:
                                                    b = 0
                                                    while b < len(dup_end_indeces):
                                                        single_value = [dup_end_indeces[b], "diff.exon"]
                                                        value.append(single_value)
                                                        b += 1
                                    m_single_exons_dict.update({y:value})
                        y += 1
                        x += 2
                    converted_m_dict_exons = {}
                    for k in m_single_exons_dict:
                        single_k = m_single_exons_dict[k]
                        if len(single_k) == len(single_m_counts) and isinstance(single_k[0], list) == False:
                            if k in converted_m_dict_exons:
                                converted_m_dict_exons[k].append(single_k)
                            elif k not in converted_m_dict_exons:
                                converted_m_dict_exons.update({k:[single_k]})
                        else:
                            sorted_k = sorted(single_k, key=lambda x: x[0])
                            if len(sorted_k) == len(single_m_counts):
                                final_value = []
                                for val in sorted_k:
                                    final_value.append(val[1])
                                if k in converted_m_dict_exons:
                                    converted_m_dict_exons[k].append(final_values)
                                elif k not in converted_m_dict_exons:
                                    converted_m_dict_exons.update({k:[final_value]})
                            else:
                                exons = [c[0] for c in sorted_k]
                                exon_values = [c[1] for c in sorted_k]
                                set_exons = list(set(exons))
                                c_more_than_1 = []
                                for c in set_exons:
                                    c_count = exons.count(c)
                                    if c_count > 1:
                                        c_more_than_1.append(c)
                                if len(c_more_than_1) == 1:
                                    if c_more_than_1[0] == 0:
                                        final_value = ["diff.exon", exon_values[1], exon_values[2]]
                                        if k in converted_m_dict_exons:
                                            converted_m_dict_exons[k].append(final_value)
                                        elif k not in converted_m_dict_exons:
                                            converted_m_dict_exons.update({k:[final_value]})
                                    elif c_more_than_1[0] == 1:
                                        final_value = [exon_values[0], "diff.exon", exon_values[2]]
                                        if k in converted_m_dict_exons:
                                            converted_m_dict_exons[k].append(final_value)
                                        elif k not in converted_m_dict_exons:
                                            converted_m_dict_exons.update({k:[final_value]})
                                    elif c_more_than_1[0] == 2:
                                        final_value = [exon_values[0], exon_values[1], "diff.exon"]
                                        if k in converted_m_dict_exons:
                                            converted_m_dict_exons[k].append(final_value)
                                        elif k not in converted_m_dict_exons:
                                            converted_m_dict_exons.update({k:[final_value]})
                                else:
                                    if exons[0] == 2:
                                        w = 0
                                        final_value = []
                                        while w <= exons[len(exons)-1]:
                                            if w in exons:
                                                w_index = exons.index(w)
                                                final_value.append(exon_values[w_index])
                                            else:
                                                final_value.append("same.exon")
                                            w += 1
                                        if k in converted_m_dict_exons:
                                            converted_m_dict_exons[k].append(final_value)
                                        elif k not in converted_m_dict_exons:
                                            converted_m_dict_exons.update({k:[final_value]})
                                    else:
                                        final_value = ["diff.end", "same.exon", "diff.end", "same.exon", "same.exon"]
                                        if k in converted_m_dict_exons:
                                            converted_m_dict_exons[k].append(final_value)
                                        elif k not in converted_m_dict_exons:
                                            converted_m_dict_exons.update({k:[final_value]})
                    for index, isoform in enumerate(single_m_counts):
                        single_isoform_values = []
                        for exon in converted_m_dict_exons:
                            single_exon = converted_m_dict_exons[exon][0]
                            single_isoform_values.append(single_exon[index])
                        if single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[1:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[1:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and list(set(single_isoform_values[0:len(single_isoform_values)-1])) == ["same.exon"]:
                            final = "%s\t%s\tMale\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (isoform,gene)
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values[0:len(single_isoform_values)-1]))) > 1:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(single_isoform_values) == 3:
                            descriptor_count = 0
                            for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values[0:len(single_isoform_values)]:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
                        elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                            descriptor_count = 0
                            for descriptor in single_isoform_values:
                                if descriptor != "same.exon":
                                    descriptor_count += 1
                            final = "%s\t%s\tMale\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (isoform,gene, str(descriptor_count))
                            if isoform in characterization_dict:
                                characterization_dict[isoform].append(final)
                            elif isoform not in characterization_dict:
                                characterization_dict.update({isoform:[final]})
    return characterization_dict


#remove shared isoforms from characterization dict
#removes isoforms that are shared between the sexes
#this produced 5999 isoforms in the characterization dict which is higher than it should be; need to filter based on female and male specific isoforms
def remove_shared():
    characterization_dict = characterize_differences()
    shared_isoforms = read_shared_isoforms()
    for key in shared_isoforms:
        single_key = shared_isoforms[key]
        for iso in single_key:
            if iso in characterization_dict:
                del characterization_dict[iso]
    return characterization_dict

#pull female specific isoforms
#PB.13104.5 was not present in the characterization dictionary
#I'll do that one by hand
def pull_characterized_female_specific_isoforms():
    female_isoforms = read_female_isoforms()
    characterization_dict = remove_shared()
    female_characterized_dict = {}
    for key in female_isoforms:
        single_key = female_isoforms[key]
        for iso in single_key:
            if iso in characterization_dict:
                female_characterized_dict.update({iso:characterization_dict[iso]})
    final = "PB.13104.4\tENSGACG00000002172\tFemale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TSS\n"
    female_characterized_dict.update({"PB.13104.5":[final]})
    return female_characterized_dict

#27 isoforms are not present in the characterization dict
# PB.5797.1
#will also need to do these by hand
def pull_characterized_male_specific_isoforms():
    male_isoforms = read_male_isoforms()
    characterization_dict = remove_shared()
    male_characterized_dict = {}
    for key in male_isoforms:
        single_key = male_isoforms[key]
        for iso in single_key:
            if iso in characterization_dict:
                male_characterized_dict.update({iso:characterization_dict[iso]})
    male_characterized_dict.update({"PB.5788.4":["PB.5788.4\tENSGACG00000020291\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.4490.3":["PB.4490.3\tENSGACG00000002358\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.9109.4":["PB.9109.4\tENSGACG00000007457\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS\n"]})
    male_characterized_dict.update({"PB.5105.5":["PB.5105.5\tENSGACG00000011659\tMale\tAltenative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.11733.3":["PB.11733.3\tENSGACG00000016040\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.9739.10":["PB.9739.10\tENSGACG00000014523\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t2.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.1058.22":["PB.1058.22\tENSGACG00000015537\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.17765.6":["PB.17765.6\tENSGACG00000015033\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.12948.6":["PB.12948.6\tENSGACG00000012974\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.14327.7":["PB.14327.7\tENSGACG00000005382\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t2.exons.alternatively.spliced,Alternate.TSS\n"]})
    male_characterized_dict.update({"PB.5083.11":["PB.5083.11\tENSGACG00000011467\tMale\tTSS.TTS.UTR.diffs\tAlternate.TTS\n"]})
    male_characterized_dict.update({"PB.5922.4":["PB.5922.4\tENSGACG00000020441\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TSS\n"]})
    male_characterized_dict.update({"PB.4359.5":["PB.4359.5\tnovelGene_3423\tMale\tTSS.TTS.UTR.diffs\tSingle.Exon.Both.TSS.TTS.different\n"]})
    male_characterized_dict.update({"PB.284.4":["PB.284.4\tENSGACG00000007985\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t2.exons.alternatively.spliced,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.6893.1":["PB.6893.1\tENSGACG00000011040\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.14575.7":["PB.14575.7\tENSGACG00000009681\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.10275.4":["PB.10275.4\tENSGACG00000008944\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t4.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.11706.2":["PB.11706.2\tENSGACG00000015928\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.1004.1":["PB.1004.1\tENSGACG00000015295\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t3.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.14487.2":["PB.14487.2\tENSGACG00000008281\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.2711.11":["PB.2711.11\tENSGACG00000016220\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.10226.3":["PB.10226.3\tENSGACG00000008035\tMale\tAlternative.Splicing\t3.exons.alternatively.spliced\n"]})
    male_characterized_dict.update({"PB.2797.2":["PB.2797.2\tENSGACG00000016706\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.1692.7":["PB.1692.7\tENSGACG00000016702\tMale\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.15475.4":["PB.15475.4\tENSGACG00000010346\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t3.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.9546.1":["PB.9546.1\tENSGACG00000012896\tMale\tAlternative.Splicing,TSS.TTS.UTR.diffs\t1.exon.alternatively.spliced,Alternate.TSS,Alternate.TTS\n"]})
    male_characterized_dict.update({"PB.5797.1":["PB.5797.1\tnovelGene_4599\tMale\tTSS.TTS.UTR.diffs\tSingle.Exon.Both.TSS.TTS.different\n"]})
    return male_characterized_dict

#write output
def write_female_isoforms():
    female_isoforms = pull_characterized_female_specific_isoforms()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for key in female_isoforms:
            single_key = female_isoforms[key]
            out.write(single_key[0])

def write_male_isoforms():
    male_isoforms = pull_characterized_male_specific_isoforms()
    output = sys.argv[7]
    with open(output, 'a') as out:
        for key in male_isoforms:
            single_key = male_isoforms[key]
            if len(single_key) == 1:
                out.write(single_key[0])



#call functions
def call():
    females = write_female_isoforms()
    males = write_male_isoforms()

call()

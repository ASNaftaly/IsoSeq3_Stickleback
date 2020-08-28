#will further categorize isoforms that are tissue specific due to splicing
#will need to read in several files
#to run script: python3 Categorize.Tissue.Specific.Splice.Variants.py <brain filtered splicing file> <brain unfiltered splicing file> <liver filtered splicing file> <liver unfiltered splicing file> <pronephros filtered splicing file> <pronephros unfiltered splicing file> <testis filtered splicing file> <testis unfiltered splicing file> <ovary filtered splicing file> <ovary unfiltered splicing file> <full isoseq gtf file with exon positions> <output 1; brain characterization> <output 2; liver characterization> <output 3; pronephros characterization> <output 4; testis characterization> <output 5; ovary characterization>
#Author: Alice Naftaly, Aug 2020

import sys

#read in brain specific isoforms by splicing
#returns dictionary with key == gene and value == isoforms
def pull_brain_specific_isoforms():
    brain_specific_file = sys.argv[1]
    brain_genes_file = sys.argv[2]
    brain_isoforms_list = []
    brain_specific_dict_withdups = {}
    final_brain_specific_dict = {}
    with open(brain_specific_file, 'r') as isoforms_file, open(brain_genes_file, 'r') as genes_file:
        for line in isoforms_file:
            isoform = line.strip("\n")
            brain_isoforms_list.append(isoform)
        for line_2 in genes_file:
            new_line = line_2.split()
            gene = new_line[0]
            isoform_id = new_line[1]
            if isoform_id in brain_isoforms_list:
                if gene in brain_specific_dict_withdups:
                    brain_specific_dict_withdups[gene].append(isoform_id)
                elif gene not in brain_specific_dict_withdups:
                    brain_specific_dict_withdups.update({gene:[isoform_id]})
    for key in brain_specific_dict_withdups:
        single_key = list(set(brain_specific_dict_withdups[key]))
        final_brain_specific_dict.update({key:single_key})
    return final_brain_specific_dict


#read in liver specific isoforms by splicing
#returns dictionary with key == gene and value == isoforms
def pull_liver_specific_isoforms():
    liver_specific_file = sys.argv[3]
    liver_genes_file = sys.argv[4]
    liver_isoforms_list = []
    liver_specific_dict_withdups = {}
    final_liver_specific_dict = {}
    with open(liver_specific_file, 'r') as isoforms_file, open(liver_genes_file, 'r') as genes_file:
        for line in isoforms_file:
            isoform = line.strip("\n")
            liver_isoforms_list.append(isoform)
        for line_2 in genes_file:
            new_line = line_2.split()
            gene = new_line[0]
            isoform_id = new_line[1]
            if isoform_id in liver_isoforms_list:
                if gene in liver_specific_dict_withdups:
                    liver_specific_dict_withdups[gene].append(isoform_id)
                elif gene not in liver_specific_dict_withdups:
                    liver_specific_dict_withdups.update({gene:[isoform_id]})
    for key in liver_specific_dict_withdups:
        single_key = list(set(liver_specific_dict_withdups[key]))
        final_liver_specific_dict.update({key:single_key})
    return final_liver_specific_dict


#read in pronephros specific isoforms by splicing
#returns dictionary with key == gene and value == isoforms
def pull_pronephros_specific_isoforms():
    pronephros_specific_file = sys.argv[5]
    pronephros_genes_file = sys.argv[6]
    pronephros_isoforms_list = []
    pronephros_specific_dict_withdups = {}
    final_pronephros_specific_dict = {}
    with open(pronephros_specific_file, 'r') as isoforms_file, open(pronephros_genes_file, 'r') as genes_file:
        for line in isoforms_file:
            isoform = line.strip("\n")
            pronephros_isoforms_list.append(isoform)
        for line_2 in genes_file:
            new_line = line_2.split()
            gene = new_line[0]
            isoform_id = new_line[1]
            if isoform_id in pronephros_isoforms_list:
                if gene in pronephros_specific_dict_withdups:
                    pronephros_specific_dict_withdups[gene].append(isoform_id)
                elif gene not in pronephros_specific_dict_withdups:
                    pronephros_specific_dict_withdups.update({gene:[isoform_id]})
    for key in pronephros_specific_dict_withdups:
        single_key = list(set(pronephros_specific_dict_withdups[key]))
        final_pronephros_specific_dict.update({key:single_key})
    return final_pronephros_specific_dict

#read in testis specific isoforms by splicing
#returns dictionary with key == gene and value == isoforms
def pull_testis_specific_isoforms():
    testis_specific_file = sys.argv[7]
    testis_genes_file = sys.argv[8]
    testis_isoforms_list = []
    testis_specific_dict_withdups = {}
    final_testis_specific_dict = {}
    with open(testis_specific_file, 'r') as isoforms_file, open(testis_genes_file, 'r') as genes_file:
        for line in isoforms_file:
            isoform = line.strip("\n")
            testis_isoforms_list.append(isoform)
        for line_2 in genes_file:
            new_line = line_2.split()
            gene = new_line[0]
            isoform_id = new_line[1]
            if isoform_id in testis_isoforms_list:
                if gene in testis_specific_dict_withdups:
                    testis_specific_dict_withdups[gene].append(isoform_id)
                elif gene not in testis_specific_dict_withdups:
                    testis_specific_dict_withdups.update({gene:[isoform_id]})
    for key in testis_specific_dict_withdups:
        single_key = list(set(testis_specific_dict_withdups[key]))
        final_testis_specific_dict.update({key:single_key})
    return final_testis_specific_dict

#read in ovary specific isoforms by splicing
#returns dictionary with key == gene and value == isoforms
def pull_ovary_specific_isoforms():
    ovary_specific_file = sys.argv[9]
    ovary_genes_file = sys.argv[10]
    ovary_isoforms_list = []
    ovary_specific_dict_withdups = {}
    final_ovary_specific_dict = {}
    with open(ovary_specific_file, 'r') as isoforms_file, open(ovary_genes_file, 'r') as genes_file:
        for line in isoforms_file:
            isoform = line.strip("\n")
            ovary_isoforms_list.append(isoform)
        for line_2 in genes_file:
            new_line = line_2.split()
            gene = new_line[0]
            isoform_id = new_line[1]
            if isoform_id in ovary_isoforms_list:
                if gene in ovary_specific_dict_withdups:
                    ovary_specific_dict_withdups[gene].append(isoform_id)
                elif gene not in ovary_specific_dict_withdups:
                    ovary_specific_dict_withdups.update({gene:[isoform_id]})
    for key in ovary_specific_dict_withdups:
        single_key = list(set(ovary_specific_dict_withdups[key]))
        final_ovary_specific_dict.update({key:single_key})
    return final_ovary_specific_dict

#read in isoseq gtf file
#returns a dictionary with key == isoform id and value == list of exons with each exon == [strand, start pos, end pos]
#minus strand isoforms are in the correct order for exon number but the start pos is the 2 value in the list for each individual exon
def read_isoseq_gtf():
    gtf_file = sys.argv[11]
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


#pull isoforms to compare brain specific isoforms to
#returns dictionary with key == gene and value == isoforms
def pull_brain_comparisons():
    brain_specific_isoforms = pull_brain_specific_isoforms()
    liver_specific_isoforms = pull_liver_specific_isoforms()
    pronephros_specific_isoforms = pull_pronephros_specific_isoforms()
    testis_specific_isoforms = pull_testis_specific_isoforms()
    ovary_specific_isoforms = pull_ovary_specific_isoforms()
    brain_comparisons_dict = {}
    for gene in brain_specific_isoforms:
        if gene in liver_specific_isoforms:
            single_liver_gene = liver_specific_isoforms[gene]
            for l_iso in single_liver_gene:
                if gene in brain_comparisons_dict:
                    brain_comparisons_dict[gene].append(l_iso)
                elif gene not in brain_comparisons_dict:
                    brain_comparisons_dict.update({gene:[l_iso]})
        if gene in pronephros_specific_isoforms:
            single_pronephros_gene = pronephros_specific_isoforms[gene]
            for p_iso in single_pronephros_gene:
                if gene in brain_comparisons_dict:
                    brain_comparisons_dict[gene].append(p_iso)
                elif gene not in brain_comparisons_dict:
                    brain_comparisons_dict.update({gene:[p_iso]})
        if gene in testis_specific_isoforms:
            single_testis_gene = testis_specific_isoforms[gene]
            for t_iso in single_testis_gene:
                if gene in brain_comparisons_dict:
                    brain_comparisons_dict[gene].append(t_iso)
                elif gene not in brain_comparisons_dict:
                    brain_comparisons_dict.update({gene:[t_iso]})
        if gene in ovary_specific_isoforms:
            single_ovary_gene = ovary_specific_isoforms[gene]
            for o_iso in single_ovary_gene:
                if gene in brain_comparisons_dict:
                    brain_comparisons_dict[gene].append(o_iso)
                elif gene not in brain_comparisons_dict:
                    brain_comparisons_dict.update({gene:[o_iso]})
    return brain_comparisons_dict


#pull isoforms to compare liver specific isoforms to
#returns dictionary with key == gene and value == isoforms
def pull_liver_comparisons():
    brain_specific_isoforms = pull_brain_specific_isoforms()
    liver_specific_isoforms = pull_liver_specific_isoforms()
    pronephros_specific_isoforms = pull_pronephros_specific_isoforms()
    testis_specific_isoforms = pull_testis_specific_isoforms()
    ovary_specific_isoforms = pull_ovary_specific_isoforms()
    liver_comparisons_dict = {}
    for gene in liver_specific_isoforms:
        if gene in brain_specific_isoforms:
            single_brain_gene = brain_specific_isoforms[gene]
            for b_iso in single_brain_gene:
                if gene in liver_comparisons_dict:
                    liver_comparisons_dict[gene].append(b_iso)
                elif gene not in liver_comparisons_dict:
                    liver_comparisons_dict.update({gene:[b_iso]})
        if gene in pronephros_specific_isoforms:
            single_pronephros_gene = pronephros_specific_isoforms[gene]
            for p_iso in single_pronephros_gene:
                if gene in liver_comparisons_dict:
                    liver_comparisons_dict[gene].append(p_iso)
                elif gene not in liver_comparisons_dict:
                    liver_comparisons_dict.update({gene:[p_iso]})
        if gene in testis_specific_isoforms:
            single_testis_gene = testis_specific_isoforms[gene]
            for t_iso in single_testis_gene:
                if gene in liver_comparisons_dict:
                    liver_comparisons_dict[gene].append(t_iso)
                elif gene not in liver_comparisons_dict:
                    liver_comparisons_dict.update({gene:[t_iso]})
        if gene in ovary_specific_isoforms:
            single_ovary_gene = ovary_specific_isoforms[gene]
            for o_iso in single_ovary_gene:
                if gene in liver_comparisons_dict:
                    liver_comparisons_dict[gene].append(o_iso)
                elif gene not in liver_comparisons_dict:
                    liver_comparisons_dict.update({gene:[o_iso]})
    return liver_comparisons_dict

#pull isoforms to compare pronephros specific isoforms to
#returns dictionary with key == gene and value == isoforms
def pull_pronephros_comparisons():
    brain_specific_isoforms = pull_brain_specific_isoforms()
    liver_specific_isoforms = pull_liver_specific_isoforms()
    pronephros_specific_isoforms = pull_pronephros_specific_isoforms()
    testis_specific_isoforms = pull_testis_specific_isoforms()
    ovary_specific_isoforms = pull_ovary_specific_isoforms()
    pronephros_comparisons_dict = {}
    for gene in pronephros_specific_isoforms:
        if gene in brain_specific_isoforms:
            single_brain_gene = brain_specific_isoforms[gene]
            for b_iso in single_brain_gene:
                if gene in pronephros_comparisons_dict:
                    pronephros_comparisons_dict[gene].append(b_iso)
                elif gene not in pronephros_comparisons_dict:
                    pronephros_comparisons_dict.update({gene:[b_iso]})
        if gene in liver_specific_isoforms:
            single_liver_gene = liver_specific_isoforms[gene]
            for l_iso in single_liver_gene:
                if gene in pronephros_comparisons_dict:
                    pronephros_comparisons_dict[gene].append(l_iso)
                elif gene not in pronephros_comparisons_dict:
                    pronephros_comparisons_dict.update({gene:[l_iso]})
        if gene in testis_specific_isoforms:
            single_testis_gene = testis_specific_isoforms[gene]
            for t_iso in single_testis_gene:
                if gene in pronephros_comparisons_dict:
                    pronephros_comparisons_dict[gene].append(t_iso)
                elif gene not in pronephros_comparisons_dict:
                    pronephros_comparisons_dict.update({gene:[t_iso]})
        if gene in ovary_specific_isoforms:
            single_ovary_gene = ovary_specific_isoforms[gene]
            for o_iso in single_ovary_gene:
                if gene in pronephros_comparisons_dict:
                    pronephros_comparisons_dict[gene].append(o_iso)
                elif gene not in pronephros_comparisons_dict:
                    pronephros_comparisons_dict.update({gene:[o_iso]})
    return pronephros_comparisons_dict

#pull isoforms to compare testis specific isoforms to
#returns dictionary with key == gene and value == isoforms
def pull_testis_comparisons():
    brain_specific_isoforms = pull_brain_specific_isoforms()
    liver_specific_isoforms = pull_liver_specific_isoforms()
    pronephros_specific_isoforms = pull_pronephros_specific_isoforms()
    testis_specific_isoforms = pull_testis_specific_isoforms()
    ovary_specific_isoforms = pull_ovary_specific_isoforms()
    testis_comparisons_dict = {}
    for gene in testis_specific_isoforms:
        if gene in brain_specific_isoforms:
            single_brain_gene = brain_specific_isoforms[gene]
            for b_iso in single_brain_gene:
                if gene in testis_comparisons_dict:
                    testis_comparisons_dict[gene].append(b_iso)
                elif gene not in testis_comparisons_dict:
                    testis_comparisons_dict.update({gene:[b_iso]})
        if gene in liver_specific_isoforms:
            single_liver_gene = liver_specific_isoforms[gene]
            for l_iso in single_liver_gene:
                if gene in testis_comparisons_dict:
                    testis_comparisons_dict[gene].append(l_iso)
                elif gene not in testis_comparisons_dict:
                    testis_comparisons_dict.update({gene:[l_iso]})
        if gene in pronephros_specific_isoforms:
            single_pronephros_gene = pronephros_specific_isoforms[gene]
            for p_iso in single_pronephros_gene:
                if gene in testis_comparisons_dict:
                    testis_comparisons_dict[gene].append(p_iso)
                elif gene not in testis_comparisons_dict:
                    testis_comparisons_dict.update({gene:[p_iso]})
        if gene in ovary_specific_isoforms:
            single_ovary_gene = ovary_specific_isoforms[gene]
            for o_iso in single_ovary_gene:
                if gene in testis_comparisons_dict:
                    testis_comparisons_dict[gene].append(o_iso)
                elif gene not in testis_comparisons_dict:
                    testis_comparisons_dict.update({gene:[o_iso]})
    return testis_comparisons_dict

#pull isoforms to compare ovary specific isoforms to
#returns dictionary with key == gene and value == isoforms
def pull_ovary_comparisons():
    brain_specific_isoforms = pull_brain_specific_isoforms()
    liver_specific_isoforms = pull_liver_specific_isoforms()
    pronephros_specific_isoforms = pull_pronephros_specific_isoforms()
    testis_specific_isoforms = pull_testis_specific_isoforms()
    ovary_specific_isoforms = pull_ovary_specific_isoforms()
    ovary_comparisons_dict = {}
    for gene in ovary_specific_isoforms:
        if gene in brain_specific_isoforms:
            single_brain_gene = brain_specific_isoforms[gene]
            for b_iso in single_brain_gene:
                if gene in ovary_comparisons_dict:
                    ovary_comparisons_dict[gene].append(b_iso)
                elif gene not in ovary_comparisons_dict:
                    ovary_comparisons_dict.update({gene:[b_iso]})
        if gene in liver_specific_isoforms:
            single_liver_gene = liver_specific_isoforms[gene]
            for l_iso in single_liver_gene:
                if gene in ovary_comparisons_dict:
                    ovary_comparisons_dict[gene].append(l_iso)
                elif gene not in ovary_comparisons_dict:
                    ovary_comparisons_dict.update({gene:[l_iso]})
        if gene in pronephros_specific_isoforms:
            single_pronephros_gene = pronephros_specific_isoforms[gene]
            for p_iso in single_pronephros_gene:
                if gene in ovary_comparisons_dict:
                    ovary_comparisons_dict[gene].append(p_iso)
                elif gene not in ovary_comparisons_dict:
                    ovary_comparisons_dict.update({gene:[p_iso]})
        if gene in testis_specific_isoforms:
            single_testis_gene = testis_specific_isoforms[gene]
            for t_iso in single_testis_gene:
                if gene in ovary_comparisons_dict:
                    ovary_comparisons_dict[gene].append(t_iso)
                elif gene not in ovary_comparisons_dict:
                    ovary_comparisons_dict.update({gene:[t_iso]})
    return ovary_comparisons_dict

#characterizing brain isoforms:
def characterize_brain():
    brain_isoforms = pull_brain_specific_isoforms()
    isoforms_to_compare = pull_brain_comparisons()
    exon_dict = read_isoseq_gtf()
    brain_characterization_dict = {}
    for gene in brain_isoforms:
        single_gene_brain = brain_isoforms[gene]
        single_gene_comparisons = isoforms_to_compare[gene]
        comp_exons = []
        #comp dict will have key == brain isoform and value == isoforms to be compared to brain with same number of exons
        comp_dict = {}
        for comp_iso in single_gene_comparisons:
            comp_exons.append([comp_iso,len(exon_dict[comp_iso])])
        for b_iso in single_gene_brain:
            single_b_iso_num_exons = len(exon_dict[b_iso])
            for val in comp_exons:
                if single_b_iso_num_exons == val[1]:
                    if b_iso in comp_dict:
                        comp_dict[b_iso].append(val[0])
                    elif b_iso not in comp_dict:
                        comp_dict.update({b_iso:[val[0]]})
        for b_key in comp_dict:
            brain_exons = exon_dict[b_key]
            single_key_isos = comp_dict[b_key]
            b_exons_all = []
            for b_exon in brain_exons:
                b_strand = brain_exons[0][0]
                if b_strand == "+":
                    b_exon_start = b_exon[1]
                    b_exon_end = b_exon[2]
                elif b_strand == "-":
                    b_exon_start = b_exon[2]
                    b_exon_end = b_exon[1]
                b_exons_all.append(b_exon_start)
                b_exons_all.append(b_exon_end)
            comp_exons = []
            for comp_key in single_key_isos:
                all_comp_exons = exon_dict[comp_key]
                comp_strand = all_comp_exons[0][0]
                for comp_exon in all_comp_exons:
                    if comp_strand == "+":
                        comp_exon_start = comp_exon[1]
                        comp_exon_end = comp_exon[2]
                    elif comp_strand == "-":
                        comp_exon_start = comp_exon[2]
                        comp_exon_end = comp_exon[1]
                    comp_exons.append(comp_exon_start)
                    comp_exons.append(comp_exon_end)
            single_b_iso_dict = {}
            x = 0
            y = 0
            while x < len(brain_exons)*2:
                single_comp_exon_start = comp_exons[x::len(brain_exons)*2]
                single_comp_exon_end = comp_exons[x+1::len(brain_exons)*2]
                single_brain_start = b_exons_all[x]
                single_brain_end = b_exons_all[x+1]
                if single_brain_start in single_comp_exon_start and single_brain_end in single_comp_exon_end:
                    comp_start_index = single_comp_exon_start.index(single_brain_start)
                    comp_end_index = single_comp_exon_end.index(single_brain_end)
                    if comp_start_index == comp_end_index or len(list(set(single_comp_exon_start))) == 1 or len(list(set(single_comp_exon_end))) == 1:
                        if y in single_b_iso_dict:
                            single_b_iso_dict[y].append("same.exon")
                        elif y not in single_b_iso_dict:
                            single_b_iso_dict.update({y:["same.exon"]})
                    elif comp_start_index != comp_end_index:
                        if y in single_b_iso_dict:
                            single_b_iso_dict[y].append("diff.exon")
                        elif y not in single_b_iso_dict:
                            single_b_iso_dict.update({y:["diff.exon"]})
                elif single_brain_start in single_comp_exon_start and single_brain_end not in single_comp_exon_end:
                    if y in single_b_iso_dict:
                        single_b_iso_dict[y].append("diff.end")
                    elif y not in single_b_iso_dict:
                        single_b_iso_dict.update({y:["diff.end"]})
                elif single_brain_start not in single_comp_exon_start and single_brain_end in single_comp_exon_end:
                    if y in single_b_iso_dict:
                        single_b_iso_dict[y].append("diff.start")
                    elif y not in single_b_iso_dict:
                        single_b_iso_dict.update({y:["diff.start"]})
                elif single_brain_start not in single_comp_exon_start and single_brain_end not in single_comp_exon_end:
                    if y in single_b_iso_dict:
                        single_b_iso_dict[y].append("diff.exon")
                    elif y not in single_b_iso_dict:
                        single_b_iso_dict.update({y:["diff.exon"]})
                x += 2
                y += 1
            single_isoform_values = []
            exon_numbers = list(range(0,len(brain_exons)-1))
            if len(exon_numbers) == 0:
                single_isoform_values = single_b_iso_dict[0]
            else:
                for val in exon_numbers:
                    single_isoform_values += single_b_iso_dict[val]
            if len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                final = "%s\t%s\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.different\n" % (b_key,gene)
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.start"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TSS.different\n" % (b_key,gene)
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "same.exon":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (b_key,gene)
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (b_key,gene)
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif list(set(single_isoform_values)) == ["same.exon"]:
                final = "%s\t%s\tCheck.Manually\tShows.up.as.all.same.exons\n" % (b_key,gene)
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (b_key,gene,str(descriptor_count))
                if b_key in brain_characterization_dict:
                    brain_characterization_dict[b_key].append(final)
                elif b_key not in brain_characterization_dict:
                    brain_characterization_dict.update({b_key:[final]})
        for final_b_iso in single_gene_brain:
            if final_b_iso not in comp_dict:
                final = "%s\t%s\tAlternative.Splicing\tNo.other.isoforms.with.same.number.of.exons\n" % (final_b_iso,gene)
                if final_b_iso in brain_characterization_dict:
                    brain_characterization_dict[final_b_iso].append(final)
                elif final_b_iso not in brain_characterization_dict:
                    brain_characterization_dict.update({final_b_iso:[final]})
    return brain_characterization_dict


#characterizing liver isoforms:
def characterize_liver():
    liver_isoforms = pull_liver_specific_isoforms()
    isoforms_to_compare = pull_liver_comparisons()
    exon_dict = read_isoseq_gtf()
    liver_characterization_dict = {}
    for gene in liver_isoforms:
        single_gene_liver = liver_isoforms[gene]
        single_gene_comparisons = isoforms_to_compare[gene]
        comp_exons = []
        #comp dict will have key == brain isoform and value == isoforms to be compared to brain with same number of exons
        comp_dict = {}
        for comp_iso in single_gene_comparisons:
            comp_exons.append([comp_iso,len(exon_dict[comp_iso])])
        for l_iso in single_gene_liver:
            single_l_iso_num_exons = len(exon_dict[l_iso])
            for val in comp_exons:
                if single_l_iso_num_exons == val[1]:
                    if l_iso in comp_dict:
                        comp_dict[l_iso].append(val[0])
                    elif l_iso not in comp_dict:
                        comp_dict.update({l_iso:[val[0]]})
        for l_key in comp_dict:
            liver_exons = exon_dict[l_key]
            single_key_isos = comp_dict[l_key]
            l_exons_all = []
            for l_exon in liver_exons:
                l_strand = liver_exons[0][0]
                if l_strand == "+":
                    l_exon_start = l_exon[1]
                    l_exon_end = l_exon[2]
                elif l_strand == "-":
                    l_exon_start = l_exon[2]
                    l_exon_end = l_exon[1]
                l_exons_all.append(l_exon_start)
                l_exons_all.append(l_exon_end)
            comp_exons = []
            for comp_key in single_key_isos:
                all_comp_exons = exon_dict[comp_key]
                comp_strand = all_comp_exons[0][0]
                for comp_exon in all_comp_exons:
                    if comp_strand == "+":
                        comp_exon_start = comp_exon[1]
                        comp_exon_end = comp_exon[2]
                    elif comp_strand == "-":
                        comp_exon_start = comp_exon[2]
                        comp_exon_end = comp_exon[1]
                    comp_exons.append(comp_exon_start)
                    comp_exons.append(comp_exon_end)
            single_l_iso_dict = {}
            x = 0
            y = 0
            while x < len(liver_exons)*2:
                single_comp_exon_start = comp_exons[x::len(liver_exons)*2]
                single_comp_exon_end = comp_exons[x+1::len(liver_exons)*2]
                single_liver_start = l_exons_all[x]
                single_liver_end = l_exons_all[x+1]
                if single_liver_start in single_comp_exon_start and single_liver_end in single_comp_exon_end:
                    comp_start_index = single_comp_exon_start.index(single_liver_start)
                    comp_end_index = single_comp_exon_end.index(single_liver_end)
                    if comp_start_index == comp_end_index or len(list(set(single_comp_exon_start))) == 1 or len(list(set(single_comp_exon_end))) == 1:
                        if y in single_l_iso_dict:
                            single_l_iso_dict[y].append("same.exon")
                        elif y not in single_l_iso_dict:
                            single_l_iso_dict.update({y:["same.exon"]})
                    elif comp_start_index != comp_end_index:
                        if y in single_l_iso_dict:
                            single_l_iso_dict[y].append("diff.exon")
                        elif y not in single_l_iso_dict:
                            single_l_iso_dict.update({y:["diff.exon"]})
                elif single_liver_start in single_comp_exon_start and single_liver_end not in single_comp_exon_end:
                    if y in single_l_iso_dict:
                        single_l_iso_dict[y].append("diff.end")
                    elif y not in single_l_iso_dict:
                        single_l_iso_dict.update({y:["diff.end"]})
                elif single_liver_start not in single_comp_exon_start and single_liver_end in single_comp_exon_end:
                    if y in single_l_iso_dict:
                        single_l_iso_dict[y].append("diff.start")
                    elif y not in single_l_iso_dict:
                        single_l_iso_dict.update({y:["diff.start"]})
                elif single_liver_start not in single_comp_exon_start and single_liver_end not in single_comp_exon_end:
                    if y in single_l_iso_dict:
                        single_l_iso_dict[y].append("diff.exon")
                    elif y not in single_l_iso_dict:
                        single_l_iso_dict.update({y:["diff.exon"]})
                x += 2
                y += 1
            single_isoform_values = []
            exon_numbers = list(range(0,len(liver_exons)-1))
            if len(exon_numbers) == 0:
                single_isoform_values = single_l_iso_dict[0]
            else:
                for val in exon_numbers:
                    single_isoform_values += single_l_iso_dict[val]
            if len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                final = "%s\t%s\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.different\n" % (l_key,gene)
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.start"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TSS.different\n" % (l_key,gene)
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.end" and single_isoform_values[1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (l_key,gene)
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.exon" and list(set(single_isoform_values[1:len(single_isoform_values)-1])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (l_key,gene)
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif list(set(single_isoform_values)) == ["same.exon"]:
                final = "%s\t%s\tCheck.Manually\tShows.up.as.all.same.exons\n" % (l_key,gene)
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon" and len(list(set(single_isoform_values))) > 1:
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (l_key,gene,str(descriptor_count))
                if l_key in liver_characterization_dict:
                    liver_characterization_dict[l_key].append(final)
                elif l_key not in liver_characterization_dict:
                    liver_characterization_dict.update({l_key:[final]})
        for final_l_iso in single_gene_liver:
            if final_l_iso not in comp_dict:
                final = "%s\t%s\tAlternative.Splicing\tNo.other.isoforms.with.same.number.of.exons\n" % (final_l_iso,gene)
                if final_l_iso in liver_characterization_dict:
                    liver_characterization_dict[final_l_iso].append(final)
                elif final_l_iso not in liver_characterization_dict:
                    liver_characterization_dict.update({final_l_iso:[final]})
    return liver_characterization_dict


#characterizing pronephros isoforms:
def characterize_pronephros():
    pronephros_isoforms = pull_pronephros_specific_isoforms()
    isoforms_to_compare = pull_pronephros_comparisons()
    exon_dict = read_isoseq_gtf()
    pronephros_characterization_dict = {}
    for gene in pronephros_isoforms:
        single_gene_pronephros = pronephros_isoforms[gene]
        single_gene_comparisons = isoforms_to_compare[gene]
        comp_exons = []
        #comp dict will have key == brain isoform and value == isoforms to be compared to brain with same number of exons
        comp_dict = {}
        for comp_iso in single_gene_comparisons:
            comp_exons.append([comp_iso,len(exon_dict[comp_iso])])
        for p_iso in single_gene_pronephros:
            single_p_iso_num_exons = len(exon_dict[p_iso])
            for val in comp_exons:
                if single_p_iso_num_exons == val[1]:
                    if p_iso in comp_dict:
                        comp_dict[p_iso].append(val[0])
                    elif p_iso not in comp_dict:
                        comp_dict.update({p_iso:[val[0]]})
        for p_key in comp_dict:
            pronephros_exons = exon_dict[p_key]
            single_key_isos = comp_dict[p_key]
            p_exons_all = []
            for p_exon in pronephros_exons:
                p_strand = pronephros_exons[0][0]
                if p_strand == "+":
                    p_exon_start = p_exon[1]
                    p_exon_end = p_exon[2]
                elif p_strand == "-":
                    p_exon_start = p_exon[2]
                    p_exon_end = p_exon[1]
                p_exons_all.append(p_exon_start)
                p_exons_all.append(p_exon_end)
            comp_exons = []
            for comp_key in single_key_isos:
                all_comp_exons = exon_dict[comp_key]
                comp_strand = all_comp_exons[0][0]
                for comp_exon in all_comp_exons:
                    if comp_strand == "+":
                        comp_exon_start = comp_exon[1]
                        comp_exon_end = comp_exon[2]
                    elif comp_strand == "-":
                        comp_exon_start = comp_exon[2]
                        comp_exon_end = comp_exon[1]
                    comp_exons.append(comp_exon_start)
                    comp_exons.append(comp_exon_end)
            single_p_iso_dict = {}
            x = 0
            y = 0
            while x < len(pronephros_exons)*2:
                single_comp_exon_start = comp_exons[x::len(pronephros_exons)*2]
                single_comp_exon_end = comp_exons[x+1::len(pronephros_exons)*2]
                single_pronephros_start = p_exons_all[x]
                single_pronephros_end = p_exons_all[x+1]
                if single_pronephros_start in single_comp_exon_start and single_pronephros_end in single_comp_exon_end:
                    comp_start_index = single_comp_exon_start.index(single_pronephros_start)
                    comp_end_index = single_comp_exon_end.index(single_pronephros_end)
                    if comp_start_index == comp_end_index or len(list(set(single_comp_exon_start))) == 1 or len(list(set(single_comp_exon_end))) == 1:
                        if y in single_p_iso_dict:
                            single_p_iso_dict[y].append("same.exon")
                        elif y not in single_p_iso_dict:
                            single_p_iso_dict.update({y:["same.exon"]})
                    elif comp_start_index != comp_end_index:
                        if y in single_p_iso_dict:
                            single_p_iso_dict[y].append("diff.exon")
                        elif y not in single_p_iso_dict:
                            single_p_iso_dict.update({y:["diff.exon"]})
                elif single_pronephros_start in single_comp_exon_start and single_pronephros_end not in single_comp_exon_end:
                    if y in single_p_iso_dict:
                        single_p_iso_dict[y].append("diff.end")
                    elif y not in single_p_iso_dict:
                        single_p_iso_dict.update({y:["diff.end"]})
                elif single_pronephros_start not in single_comp_exon_start and single_pronephros_end in single_comp_exon_end:
                    if y in single_p_iso_dict:
                        single_p_iso_dict[y].append("diff.start")
                    elif y not in single_p_iso_dict:
                        single_p_iso_dict.update({y:["diff.start"]})
                elif single_pronephros_start not in single_comp_exon_start and single_pronephros_end not in single_comp_exon_end:
                    if y in single_p_iso_dict:
                        single_p_iso_dict[y].append("diff.exon")
                    elif y not in single_p_iso_dict:
                        single_p_iso_dict.update({y:["diff.exon"]})
                x += 2
                y += 1
            single_isoform_values = []
            exon_numbers = list(range(0,len(pronephros_exons)-1))
            if len(exon_numbers) == 0:
                single_isoform_values = single_p_iso_dict[0]
            else:
                for val in exon_numbers:
                    single_isoform_values += single_p_iso_dict[val]
            if len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                final = "%s\t%s\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.different\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.start"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TSS.different\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.end"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TTS.different\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.end" and single_isoform_values[1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "same.exon":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "diff.end":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.end" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif list(set(single_isoform_values)) == ["same.exon"]:
                final = "%s\t%s\tCheck.Manually\tShows.up.as.all.same.exons\n" % (p_key,gene)
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TTS\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon" and len(list(set(single_isoform_values))) > 1:
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (p_key,gene,str(descriptor_count))
                if p_key in pronephros_characterization_dict:
                    pronephros_characterization_dict[p_key].append(final)
                elif p_key not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({p_key:[final]})
        for final_p_iso in single_gene_pronephros:
            if final_p_iso not in comp_dict:
                final = "%s\t%s\tAlternative.Splicing\tNo.other.isoforms.with.same.number.of.exons\n" % (final_p_iso,gene)
                if final_p_iso in pronephros_characterization_dict:
                    pronephros_characterization_dict[final_p_iso].append(final)
                elif final_p_iso not in pronephros_characterization_dict:
                    pronephros_characterization_dict.update({final_p_iso:[final]})
    return pronephros_characterization_dict

#characterizing testis isoforms:
def characterize_testis():
    testis_isoforms = pull_testis_specific_isoforms()
    isoforms_to_compare = pull_testis_comparisons()
    exon_dict = read_isoseq_gtf()
    testis_characterization_dict = {}
    for gene in testis_isoforms:
        single_gene_testis = testis_isoforms[gene]
        single_gene_comparisons = isoforms_to_compare[gene]
        comp_exons = []
        #comp dict will have key == brain isoform and value == isoforms to be compared to brain with same number of exons
        comp_dict = {}
        for comp_iso in single_gene_comparisons:
            comp_exons.append([comp_iso,len(exon_dict[comp_iso])])
        for t_iso in single_gene_testis:
            single_t_iso_num_exons = len(exon_dict[t_iso])
            for val in comp_exons:
                if single_t_iso_num_exons == val[1]:
                    if t_iso in comp_dict:
                        comp_dict[t_iso].append(val[0])
                    elif t_iso not in comp_dict:
                        comp_dict.update({t_iso:[val[0]]})
        for t_key in comp_dict:
            testis_exons = exon_dict[t_key]
            single_key_isos = comp_dict[t_key]
            t_exons_all = []
            for t_exon in testis_exons:
                t_strand = testis_exons[0][0]
                if t_strand == "+":
                    t_exon_start = t_exon[1]
                    t_exon_end = t_exon[2]
                elif t_strand == "-":
                    t_exon_start = t_exon[2]
                    t_exon_end = t_exon[1]
                t_exons_all.append(t_exon_start)
                t_exons_all.append(t_exon_end)
            comp_exons = []
            for comp_key in single_key_isos:
                all_comp_exons = exon_dict[comp_key]
                comp_strand = all_comp_exons[0][0]
                for comp_exon in all_comp_exons:
                    if comp_strand == "+":
                        comp_exon_start = comp_exon[1]
                        comp_exon_end = comp_exon[2]
                    elif comp_strand == "-":
                        comp_exon_start = comp_exon[2]
                        comp_exon_end = comp_exon[1]
                    comp_exons.append(comp_exon_start)
                    comp_exons.append(comp_exon_end)
            single_t_iso_dict = {}
            x = 0
            y = 0
            while x < len(testis_exons)*2:
                single_comp_exon_start = comp_exons[x::len(testis_exons)*2]
                single_comp_exon_end = comp_exons[x+1::len(testis_exons)*2]
                single_testis_start = t_exons_all[x]
                single_testis_end = t_exons_all[x+1]
                if single_testis_start in single_comp_exon_start and single_testis_end in single_comp_exon_end:
                    comp_start_index = single_comp_exon_start.index(single_testis_start)
                    comp_end_index = single_comp_exon_end.index(single_testis_end)
                    if comp_start_index == comp_end_index or len(list(set(single_comp_exon_start))) == 1 or len(list(set(single_comp_exon_end))) == 1:
                        if y in single_t_iso_dict:
                            single_t_iso_dict[y].append("same.exon")
                        elif y not in single_t_iso_dict:
                            single_t_iso_dict.update({y:["same.exon"]})
                    elif comp_start_index != comp_end_index:
                        if y in single_t_iso_dict:
                            single_t_iso_dict[y].append("diff.exon")
                        elif y not in single_t_iso_dict:
                            single_t_iso_dict.update({y:["diff.exon"]})
                elif single_testis_start in single_comp_exon_start and single_testis_end not in single_comp_exon_end:
                    if y in single_t_iso_dict:
                        single_t_iso_dict[y].append("diff.end")
                    elif y not in single_t_iso_dict:
                        single_t_iso_dict.update({y:["diff.end"]})
                elif single_testis_start not in single_comp_exon_start and single_testis_end in single_comp_exon_end:
                    if y in single_t_iso_dict:
                        single_t_iso_dict[y].append("diff.start")
                    elif y not in single_t_iso_dict:
                        single_t_iso_dict.update({y:["diff.start"]})
                elif single_testis_start not in single_comp_exon_start and single_testis_end not in single_comp_exon_end:
                    if y in single_t_iso_dict:
                        single_t_iso_dict[y].append("diff.exon")
                    elif y not in single_t_iso_dict:
                        single_t_iso_dict.update({y:["diff.exon"]})
                x += 2
                y += 1
            single_isoform_values = []
            exon_numbers = list(range(0,len(testis_exons)-1))
            if len(exon_numbers) == 0:
                single_isoform_values = single_t_iso_dict[0]
            else:
                for val in exon_numbers:
                    single_isoform_values += single_t_iso_dict[val]
            if len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                final = "%s\t%s\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.different\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.start"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TSS.different\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.end"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TTS.different\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "same.exon":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "diff.end":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS,Alternate.TTS\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[0:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif list(set(single_isoform_values)) == ["same.exon"]:
                final = "%s\t%s\tCheck.Manually\tShows.up.as.all.same.exons\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon" and len(list(set(single_isoform_values))) > 1:
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (t_key,gene,str(descriptor_count))
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.end" and len(list(set(single_isoform_values))) > 1 and list(set(single_isoform_values[0:len(single_isoform_values)-1])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TTS\n" % (t_key,gene)
                if t_key in testis_characterization_dict:
                    testis_characterization_dict[t_key].append(final)
                elif t_key not in testis_characterization_dict:
                    testis_characterization_dict.update({t_key:[final]})
        for final_t_iso in single_gene_testis:
            if final_t_iso not in comp_dict:
                final = "%s\t%s\tAlternative.Splicing\tNo.other.isoforms.with.same.number.of.exons\n" % (final_t_iso,gene)
                if final_t_iso in testis_characterization_dict:
                    testis_characterization_dict[final_t_iso].append(final)
                elif final_t_iso not in testis_characterization_dict:
                    testis_characterization_dict.update({final_t_iso:[final]})
    return testis_characterization_dict

#characterizing pronephros isoforms:
def characterize_ovary():
    ovary_isoforms = pull_ovary_specific_isoforms()
    isoforms_to_compare = pull_ovary_comparisons()
    exon_dict = read_isoseq_gtf()
    ovary_characterization_dict = {}
    for gene in ovary_isoforms:
        single_gene_ovary = ovary_isoforms[gene]
        single_gene_comparisons = isoforms_to_compare[gene]
        comp_exons = []
        #comp dict will have key == brain isoform and value == isoforms to be compared to brain with same number of exons
        comp_dict = {}
        for comp_iso in single_gene_comparisons:
            comp_exons.append([comp_iso,len(exon_dict[comp_iso])])
        for o_iso in single_gene_ovary:
            single_o_iso_num_exons = len(exon_dict[o_iso])
            for val in comp_exons:
                if single_o_iso_num_exons == val[1]:
                    if o_iso in comp_dict:
                        comp_dict[o_iso].append(val[0])
                    elif o_iso not in comp_dict:
                        comp_dict.update({o_iso:[val[0]]})
        for o_key in comp_dict:
            ovary_exons = exon_dict[o_key]
            single_key_isos = comp_dict[o_key]
            o_exons_all = []
            for o_exon in ovary_exons:
                o_strand = ovary_exons[0][0]
                if o_strand == "+":
                    o_exon_start = o_exon[1]
                    o_exon_end = o_exon[2]
                elif o_strand == "-":
                    o_exon_start = o_exon[2]
                    o_exon_end = o_exon[1]
                o_exons_all.append(o_exon_start)
                o_exons_all.append(o_exon_end)
            comp_exons = []
            for comp_key in single_key_isos:
                all_comp_exons = exon_dict[comp_key]
                comp_strand = all_comp_exons[0][0]
                for comp_exon in all_comp_exons:
                    if comp_strand == "+":
                        comp_exon_start = comp_exon[1]
                        comp_exon_end = comp_exon[2]
                    elif comp_strand == "-":
                        comp_exon_start = comp_exon[2]
                        comp_exon_end = comp_exon[1]
                    comp_exons.append(comp_exon_start)
                    comp_exons.append(comp_exon_end)
            single_o_iso_dict = {}
            x = 0
            y = 0
            while x < len(ovary_exons)*2:
                single_comp_exon_start = comp_exons[x::len(ovary_exons)*2]
                single_comp_exon_end = comp_exons[x+1::len(ovary_exons)*2]
                single_ovary_start = o_exons_all[x]
                single_ovary_end = o_exons_all[x+1]
                if single_ovary_start in single_comp_exon_start and single_ovary_end in single_comp_exon_end:
                    comp_start_index = single_comp_exon_start.index(single_ovary_start)
                    comp_end_index = single_comp_exon_end.index(single_ovary_end)
                    if comp_start_index == comp_end_index or len(list(set(single_comp_exon_start))) == 1 or len(list(set(single_comp_exon_end))) == 1:
                        if y in single_o_iso_dict:
                            single_o_iso_dict[y].append("same.exon")
                        elif y not in single_o_iso_dict:
                            single_o_iso_dict.update({y:["same.exon"]})
                    elif comp_start_index != comp_end_index:
                        if y in single_o_iso_dict:
                            single_o_iso_dict[y].append("diff.exon")
                        elif y not in single_o_iso_dict:
                            single_o_iso_dict.update({y:["diff.exon"]})
                elif single_ovary_start in single_comp_exon_start and single_ovary_end not in single_comp_exon_end:
                    if y in single_o_iso_dict:
                        single_t_iso_dict[y].append("diff.end")
                    elif y not in single_o_iso_dict:
                        single_o_iso_dict.update({y:["diff.end"]})
                elif single_ovary_start not in single_comp_exon_start and single_ovary_end in single_comp_exon_end:
                    if y in single_o_iso_dict:
                        single_o_iso_dict[y].append("diff.start")
                    elif y not in single_o_iso_dict:
                        single_o_iso_dict.update({y:["diff.start"]})
                elif single_ovary_start not in single_comp_exon_start and single_ovary_end not in single_comp_exon_end:
                    if y in single_o_iso_dict:
                        single_o_iso_dict[y].append("diff.exon")
                    elif y not in single_o_iso_dict:
                        single_o_iso_dict.update({y:["diff.exon"]})
                x += 2
                y += 1
            single_isoform_values = []
            exon_numbers = list(range(0,len(ovary_exons)-1))
            if len(exon_numbers) == 0:
                single_isoform_values = single_o_iso_dict[0]
            else:
                for val in exon_numbers:
                    single_isoform_values += single_o_iso_dict[val]
            if len(single_isoform_values) == 1 and single_isoform_values == ["diff.exon"]:
                final = "%s\t%s\tAlternative.Splicing\tSingle.Exon.Both.TSS.TTS.different\n" % (o_key,gene)
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif len(single_isoform_values) == 1 and single_isoform_values == ["diff.start"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tSingle.Exon.TSS.different\n" % (o_key,gene)
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.start" and single_isoform_values[1] == "same.exon":
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (o_key,gene)
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "diff.start" and list(set(single_isoform_values[1:len(single_isoform_values)])) == ["same.exon"]:
                final = "%s\t%s\tTSS.TTS.UTR.diffs\tAlternate.TSS\n" % (o_key,gene)
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif len(single_isoform_values) == 2 and single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "diff.start" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "diff.start"and single_isoform_values[len(single_isoform_values)-1] == "diff.end":
                descriptor_count = 0
                for descriptor in single_isoform_values[1:len(single_isoform_values)-1]:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS,Alternate.TTS\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "same.exon":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "diff.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing,TSS.TTS.UTR.diffs\t%s.exons.alternatively.spliced,Alternate.TSS\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif single_isoform_values[0] == "same.exon" and single_isoform_values[len(single_isoform_values)-1] == "diff.start":
                descriptor_count = 0
                for descriptor in single_isoform_values:
                    if descriptor != "same.exon":
                        descriptor_count += 1
                final = "%s\t%s\tAlternative.Splicing\t%s.exons.alternatively.spliced\n" % (o_key,gene,str(descriptor_count))
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
            elif list(set(single_isoform_values)) == ["same.exon"]:
                final = "%s\t%s\tCheck.Manually\tShows.up.as.all.same.exons\n" % (o_key,gene)
                if o_key in ovary_characterization_dict:
                    ovary_characterization_dict[o_key].append(final)
                elif o_key not in ovary_characterization_dict:
                    ovary_characterization_dict.update({o_key:[final]})
        for final_o_iso in single_gene_ovary:
            if final_o_iso not in comp_dict:
                final = "%s\t%s\tAlternative.Splicing\tNo.other.isoforms.with.same.number.of.exons\n" % (final_o_iso,gene)
                if final_o_iso in ovary_characterization_dict:
                    ovary_characterization_dict[final_o_iso].append(final)
                elif final_o_iso not in ovary_characterization_dict:
                    ovary_characterization_dict.update({final_o_iso:[final]})
    return ovary_characterization_dict

#write output
def write_brain_isoforms():
    brain_isoforms = characterize_brain()
    output = sys.argv[12]
    with open(output, 'a') as out:
        for key in brain_isoforms:
            single_key = brain_isoforms[key]
            out.write(single_key[0])

def write_liver_isoforms():
    liver_isoforms = characterize_liver()
    output = sys.argv[13]
    with open(output, 'a') as out:
        for key in liver_isoforms:
            single_key = liver_isoforms[key]
            out.write(single_key[0])

def write_pronephros_isoforms():
    pronephros_isoforms = characterize_pronephros()
    output = sys.argv[14]
    with open(output, 'a') as out:
        for key in pronephros_isoforms:
            single_key = pronephros_isoforms[key]
            out.write(single_key[0])

def write_testis_isoforms():
    testis_isoforms = characterize_testis()
    output = sys.argv[15]
    with open(output, 'a') as out:
        for key in testis_isoforms:
            single_key = testis_isoforms[key]
            out.write(single_key[0])

def write_ovary_isoforms():
    ovary_isoforms = characterize_ovary()
    output = sys.argv[16]
    with open(output, 'a') as out:
        for key in ovary_isoforms:
            single_key = ovary_isoforms[key]
            out.write(single_key[0])


#call_functions
def call():
    brain = write_brain_isoforms()
    liver = write_liver_isoforms()
    pronephros = write_pronephros_isoforms()
    testis = write_testis_isoforms()
    ovary = write_ovary_isoforms()

call()

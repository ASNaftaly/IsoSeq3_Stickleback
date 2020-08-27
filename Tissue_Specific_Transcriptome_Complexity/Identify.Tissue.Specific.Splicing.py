#identifying isoforms that are tissue-specific splicing events vs tissue-specific due to differences in gene expression (present in only one tissue) vs isoforms found in two or more tissues
#tissue specific splicing events will also include TSS, TTS differences (will have a separate script to compare the specifics)
#combined sexes for somatic tissues but kept testis and ovary separate
#to run script: python3 Identify.Tissue.Specific.Splicing.py <brain specific isoforms file> <liver specific isoforms file> <pronephros specific isoforms file> <testis specific isoforms file> <ovary specific isoforms file> <output 1; brain specific AS isoforms> <output 2; brain specific expression isoforms> <output 3; liver specific AS isoforms> <output 4; liver specific expression isoforms> <output 5; pronephros specific AS isoforms> <output 6; pronephros specific expression isoforms> <output 7; testis specific AS isoforms> <output 8; testis specific expression isoforms> <output 9; ovary specific AS isoforms> <output 10; ovary specific expression isoforms>
#Author: Alice Naftaly, Aug 2020

import sys


#read in isoforms specific to brain
#returns dictionary with key == gene and value == isoform
def pull_brain_specific_isoforms():
    brain_file = sys.argv[1]
    brain_isoforms = {}
    with open(brain_file, 'r') as brain_only:
        for line in brain_only:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[1]
                if gene_id in brain_isoforms:
                    brain_isoforms[gene_id].append(isoform)
                elif gene_id not in brain_isoforms:
                    brain_isoforms.update({gene_id:[isoform]})
    return brain_isoforms

#read in isoforms specific to liver
#returns dictionary with key == gene and value == isoform
def pull_liver_specific_isoforms():
    liver_file = sys.argv[2]
    liver_isoforms = {}
    with open(liver_file, 'r') as liver_only:
        for line in liver_only:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[1]
                if gene_id in liver_isoforms:
                    liver_isoforms[gene_id].append(isoform)
                elif gene_id not in liver_isoforms:
                    liver_isoforms.update({gene_id:[isoform]})
    return liver_isoforms

#read in isoforms specific to pronephros
#returns dictionary with key == gene and value == isoform
def pull_pronephros_specific_isoforms():
    pronephros_file = sys.argv[3]
    pronephros_isoforms = {}
    with open(pronephros_file, 'r') as pronephros_only:
        for line in pronephros_only:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[1]
                if gene_id in pronephros_isoforms:
                    pronephros_isoforms[gene_id].append(isoform)
                elif gene_id not in pronephros_isoforms:
                    pronephros_isoforms.update({gene_id:[isoform]})
    return pronephros_isoforms

#read in isoforms specific to testis
#returns dictionary with key == gene and value == isoform
def pull_testis_specific_isoforms():
    testis_file = sys.argv[4]
    testis_isoforms = {}
    with open(testis_file, 'r') as testis_only:
        for line in testis_only:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[1]
                if gene_id in testis_isoforms:
                    testis_isoforms[gene_id].append(isoform)
                elif gene_id not in testis_isoforms:
                    testis_isoforms.update({gene_id:[isoform]})
    return testis_isoforms

#read in isoforms specific to ovary
#returns dictionary with key == gene and value == isoform
def pull_ovary_specific_isoforms():
    ovary_file = sys.argv[5]
    ovary_isoforms = {}
    with open(ovary_file, 'r') as ovary_only:
        for line in ovary_only:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene_id = new_line[1]
                if gene_id in ovary_isoforms:
                    ovary_isoforms[gene_id].append(isoform)
                elif gene_id not in ovary_isoforms:
                    ovary_isoforms.update({gene_id:[isoform]})
    return ovary_isoforms


#genes shared between tissues can be shared with as few as two tissues up to all tissues
def sort_isoforms():
    brain_isoforms = pull_brain_specific_isoforms()
    liver_isoforms = pull_liver_specific_isoforms()
    pronephros_isoforms = pull_pronephros_specific_isoforms()
    testis_isoforms = pull_testis_specific_isoforms()
    ovary_isoforms = pull_ovary_specific_isoforms()
    pronephros_shared_genes = {}
    pronephros_unique_genes = {}
    brain_shared_genes = {}
    brain_unique_genes = {}
    liver_shared_genes = {}
    liver_unique_genes = {}
    testis_shared_genes = {}
    testis_unique_genes = {}
    ovary_shared_genes = {}
    ovary_unique_genes = {}
    #there should be no shared isoforms as these isoforms are already specific to a tissue
    for p_gene in pronephros_isoforms:
        single_gene_pronephros = pronephros_isoforms[p_gene]
        if p_gene in brain_isoforms:
            single_gene_brain = brain_isoforms[p_gene]
            #shared gene and different isoform == alternative splicing candidate
            for p_iso in single_gene_pronephros:
                final = [p_iso, "shared.gene.with.brain"]
                if p_gene in pronephros_shared_genes:
                    pronephros_shared_genes[p_gene].append(final)
                elif p_gene not in pronephros_shared_genes:
                    pronephros_shared_genes.update({p_gene:[final]})
            for b_iso in single_gene_brain:
                final = [b_iso, "shared.gene.with.pronephros"]
                if p_gene in brain_shared_genes:
                    brain_shared_genes[p_gene].append(final)
                elif p_gene not in brain_shared_genes:
                    brain_shared_genes.update({p_gene:[final]})
        if p_gene in liver_isoforms:
            single_gene_liver = liver_isoforms[p_gene]
            for p_iso in single_gene_pronephros:
                final = [p_iso, "shared.gene.with.liver"]
                if p_gene in pronephros_shared_genes:
                    pronephros_shared_genes[p_gene].append(final)
                elif p_gene not in pronephros_shared_genes:
                    pronephros_shared_genes.update({p_gene:[final]})
            for l_iso in single_gene_liver:
                final = [l_iso, "shared.gene.with.pronephros"]
                if p_gene in liver_shared_genes:
                    liver_shared_genes[p_gene].append(final)
                elif p_gene not in liver_shared_genes:
                    liver_shared_genes.update({p_gene:[final]})
        if p_gene in testis_isoforms:
            single_gene_testis = testis_isoforms[p_gene]
            for p_iso in single_gene_pronephros:
                final = [p_iso, "shared.gene.with.testis"]
                if p_gene in pronephros_shared_genes:
                    pronephros_shared_genes[p_gene].append(final)
                elif p_gene not in pronephros_shared_genes:
                    pronephros_shared_genes.update({p_gene:[final]})
            for t_iso in single_gene_testis:
                final = [t_iso, "shared.gene.with.pronephros"]
                if p_gene in testis_shared_genes:
                    testis_shared_genes[p_gene].append(final)
                elif p_gene not in testis_shared_genes:
                    testis_shared_genes.update({p_gene:[final]})
        if p_gene in ovary_isoforms:
            single_gene_ovary = ovary_isoforms[p_gene]
            for p_iso in single_gene_pronephros:
                final = [p_iso, "shared.gene.with.ovary"]
                if p_gene in pronephros_shared_genes:
                    pronephros_shared_genes[p_gene].append(final)
                elif p_gene not in pronephros_shared_genes:
                    pronephros_shared_genes.update({p_gene:[final]})
            for o_iso in single_gene_ovary:
                final = [o_iso, "shared.gene.with.pronephros"]
                if p_gene in ovary_shared_genes:
                    ovary_shared_genes[p_gene].append(final)
                elif p_gene not in ovary_shared_genes:
                    ovary_shared_genes.update({p_gene:[final]})
    for final_p_gene in pronephros_isoforms:
        if final_p_gene not in pronephros_shared_genes:
            pronephros_unique_genes.update({final_p_gene:pronephros_isoforms[final_p_gene]})
    for b_gene in brain_isoforms:
        single_gene_brain = brain_isoforms[b_gene]
        if b_gene in liver_isoforms:
            single_gene_liver = liver_isoforms[b_gene]
            for b_iso in single_gene_brain:
                final = [b_iso, "shared.gene.with.liver"]
                if b_gene in brain_shared_genes:
                    brain_shared_genes[b_gene].append(final)
                elif b_gene not in brain_shared_genes:
                    brain_shared_genes.update({b_gene:[final]})
            for l_iso in single_gene_liver:
                final = [l_iso, "shared.gene.with.brain"]
                if b_gene in liver_shared_genes:
                    liver_shared_genes[b_gene].append(final)
                elif b_gene not in liver_shared_genes:
                    liver_shared_genes.update({b_gene:[final]})
        if b_gene in testis_isoforms:
            single_gene_testis = testis_isoforms[b_gene]
            for b_iso in single_gene_brain:
                final = [b_iso, "shared.gene.with.testis"]
                if b_gene in brain_shared_genes:
                    brain_shared_genes[b_gene].append(final)
                elif b_gene not in brain_shared_genes:
                    brain_shared_genes.update({b_gene:[final]})
            for t_iso in single_gene_testis:
                final = [t_iso, "shared.gene.with.brain"]
                if b_gene in testis_shared_genes:
                    testis_shared_genes[b_gene].append(final)
                elif b_gene not in testis_shared_genes:
                    testis_shared_genes.update({b_gene:[final]})
        if b_gene in ovary_isoforms:
            single_gene_ovary = ovary_isoforms[b_gene]
            for b_iso in single_gene_brain:
                final = [b_iso, "shared.gene.with.ovary"]
                if b_gene in brain_shared_genes:
                    brain_shared_genes[b_gene].append(final)
                elif b_gene not in brain_shared_genes:
                    brain_shared_genes.update({b_gene:[final]})
            for o_iso in single_gene_ovary:
                final = [o_iso, "shared.gene.with.brain"]
                if b_gene in ovary_shared_genes:
                    ovary_shared_genes[b_gene].append(final)
                elif b_gene not in ovary_shared_genes:
                    ovary_shared_genes.update({b_gene:[final]})
    for final_b_gene in brain_isoforms:
        if final_b_gene not in brain_shared_genes:
            brain_unique_genes.update({final_b_gene:brain_isoforms[final_b_gene]})
    for l_gene in liver_isoforms:
        single_gene_liver = liver_isoforms[l_gene]
        if l_gene in testis_isoforms:
            single_gene_testis = testis_isoforms[l_gene]
            for l_iso in single_gene_liver:
                final = [l_iso, "shared.gene.with.testis"]
                if l_gene in liver_shared_genes:
                    liver_shared_genes[l_gene].append(final)
                elif l_gene not in liver_shared_genes:
                    liver_shared_genes.update({l_gene:[final]})
            for t_iso in single_gene_testis:
                final = [t_iso, "shared.gene.with.liver"]
                if l_gene in testis_shared_genes:
                    testis_shared_genes[l_gene].append(final)
                elif l_gene not in testis_shared_genes:
                    testis_shared_genes.update({l_gene:[final]})
        if l_gene in ovary_isoforms:
            single_gene_ovary = ovary_isoforms[l_gene]
            for l_iso in single_gene_liver:
                final = [l_iso, "shared.gene.with.ovary"]
                if l_gene in liver_shared_genes:
                    liver_shared_genes[l_gene].append(final)
                elif l_gene not in liver_shared_genes:
                    liver_shared_genes.update({l_gene:[final]})
            for o_iso in single_gene_ovary:
                final = [o_iso, "shared.gene.with.liver"]
                if l_gene in ovary_shared_genes:
                    ovary_shared_genes[l_gene].append(final)
                elif l_gene not in ovary_shared_genes:
                    ovary_shared_genes.update({l_gene:[final]})
    for final_l_gene in liver_isoforms:
        if final_l_gene not in liver_shared_genes:
            liver_unique_genes.update({final_l_gene:liver_isoforms[final_l_gene]})
    for t_gene in testis_isoforms:
        single_gene_testis = testis_isoforms[t_gene]
        if t_gene in ovary_isoforms:
            single_gene_ovary = ovary_isoforms[t_gene]
            for t_iso in single_gene_testis:
                final = [t_iso, "shared.gene.with.ovary"]
                if t_gene in testis_shared_genes:
                    testis_shared_genes[t_gene].append(final)
                elif t_gene not in testis_shared_genes:
                    testis_shared_genes.update({t_gene:[final]})
            for o_iso in single_gene_ovary:
                final = [o_iso, "shared.gene.with.testis"]
                if t_gene in ovary_shared_genes:
                    ovary_shared_genes[t_gene].append(final)
                elif t_gene not in ovary_shared_genes:
                    ovary_shared_genes.update({t_gene:[final]})
    for final_t_gene in testis_isoforms:
        if final_t_gene not in testis_shared_genes:
            testis_unique_genes.update({final_t_gene:testis_isoforms[final_t_gene]})
    for o_gene in ovary_isoforms:
        if o_gene not in ovary_shared_genes:
            ovary_unique_genes.update({o_gene:ovary_isoforms[o_gene]})
    return brain_shared_genes, brain_unique_genes, liver_shared_genes, liver_unique_genes, pronephros_shared_genes, pronephros_unique_genes, testis_shared_genes, testis_unique_genes, ovary_shared_genes, ovary_unique_genes

#write isoforms to output file
#each output file will have 2 columns:
#isoform.id \t gene id \n
def write_output():
    brain_shared_genes, brain_unique_genes, liver_shared_genes, liver_unique_genes, pronephros_shared_genes, pronephros_unique_genes, testis_shared_genes, testis_unique_genes, ovary_shared_genes, ovary_unique_genes = sort_isoforms()
    brain_shared_output = sys.argv[6]
    brain_unique_output = sys.argv[7]
    liver_shared_output = sys.argv[8]
    liver_unique_output = sys.argv[9]
    pronephros_shared_output = sys.argv[10]
    pronephros_unique_output = sys.argv[11]
    testis_shared_output = sys.argv[12]
    testis_unique_output = sys.argv[13]
    ovary_shared_output = sys.argv[14]
    ovary_unique_output = sys.argv[15]
    with open(brain_shared_output, 'a') as out1, open(brain_unique_output, 'a') as out2, open(liver_shared_output, 'a') as out3, open(liver_unique_output, 'a') as out4, open(pronephros_shared_output, 'a') as out5, open(pronephros_unique_output, 'a') as out6, open(testis_shared_output, 'a') as out7, open(testis_unique_output, 'a') as out8, open(ovary_shared_output, 'a') as out9, open(ovary_unique_output, 'a') as out10:
        for key in brain_shared_genes:
            single_key = brain_shared_genes[key]
            for v in single_key:
                final = "%s\t%s\t%s\n" % (str(key), str(v[0]), str(v[1]))
                out1.write(final)
        for key2 in brain_unique_genes:
            single_key = brain_unique_genes[key2]
            for v in single_key:
                final = "%s\t%s\n" % (str(key2), str(v))
                out2.write(final)
        for key3 in liver_shared_genes:
            single_key = liver_shared_genes[key3]
            for v in single_key:
                final = "%s\t%s\t%s\n" % (str(key3), str(v[0]), str(v[1]))
                out3.write(final)
        for key4 in liver_unique_genes:
            single_key = liver_unique_genes[key4]
            for v in single_key:
                final = "%s\t%s\n" % (str(key4), str(v))
                out4.write(final)
        for key5 in pronephros_shared_genes:
            single_key = pronephros_shared_genes[key5]
            for v in single_key:
                final = "%s\t%s\t%s\n" % (str(key5), str(v[0]), str(v[1]))
                out5.write(final)
        for key6 in pronephros_unique_genes:
            single_key = pronephros_unique_genes[key6]
            for v in single_key:
                final = "%s\t%s\n" % (str(key6), str(v))
                out6.write(final)
        for key7 in testis_shared_genes:
            single_key = testis_shared_genes[key7]
            for v in single_key:
                final = "%s\t%s\t%s\n" % (str(key7),str(v[0]), str(v[1]))
                out7.write(final)
        for key8 in testis_unique_genes:
            single_key = testis_unique_genes[key8]
            for v in single_key:
                final = "%s\t%s\n" % (str(key8), str(v))
                out8.write(final)
        for key9 in ovary_shared_genes:
            single_key = ovary_shared_genes[key9]
            for v in single_key:
                final = "%s\t%s\t%s\n" % (str(key9), str(v[0]), str(v[1]))
                out9.write(final)
        for key10 in ovary_unique_genes:
            single_key = ovary_unique_genes[key10]
            for v in single_key:
                final = "%s\t%s\n" % (str(key10), str(v))
                out10.write(final)
write_output()

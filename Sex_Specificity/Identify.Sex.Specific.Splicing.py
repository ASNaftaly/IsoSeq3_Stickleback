#identifying isoforms that are sex-specific splicing events vs sex-specific due to differences in gene expression (present in only one sex) vs isoforms found in both sexes
#sex specific splicing events will also include TSS, TTS differences (will have a separate script to compare the specifics)
#to run script: python3 Identify.Sex.Specific.Splicing.py <female isoforms exon counts file; has isoform id in one column and exon count in the other> <male isoforms exon counts; same format as female exon counts file> <all isoform ids file; file has two columns with isoform id and gene id> <output 1; shared isoforms> <output 2; male specific isoforms by expression> <output 3; female specific isoforms by expression> <output 4; male specific isoforms by splicing> <output 5; female specific isoforms by splicing>
#Author: Alice Naftaly, Aug 2020

import sys

#read in isoforms from female exon counts file:
#returns list of isoform ids
def pull_female_specific_isoforms():
    female_file = sys.argv[1]
    female_isoforms = []
    with open(female_file, 'r') as female_only:
        for line in female_only:
            if line.startswith("PB"):
                new_line = line.split()
                female_isoforms.append(new_line[0])
    return female_isoforms


#read in isoforms from male exon counts file:
#returns list of isoform ids
def pull_male_specific_isoforms():
    male_file = sys.argv[2]
    male_isoforms = []
    with open(male_file, 'r') as male_only:
        for line in male_only:
            if line.startswith("PB"):
                new_line = line.split()
                male_isoforms.append(new_line[0])
    return male_isoforms


#read all isoforms
#returns dictionary with key == gene and value == isoform ids that match that gene
def pull_all_isoforms():
    all_file = sys.argv[3]
    all_isoforms = {}
    with open(all_file, 'r') as all:
        for line in all:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                gene = new_line[1]
                if gene in all_isoforms:
                    all_isoforms[gene].append(isoform)
                elif gene not in all_isoforms:
                    all_isoforms.update({gene:[isoform]})
    return all_isoforms

#convert female and male isoform lists to dictionary that matches all isoforms dictionary
def convert_female_list():
    all_isoforms = pull_all_isoforms()
    female_isoforms = pull_female_specific_isoforms()
    converted_female_dict = {}
    for gene in all_isoforms:
        single_gene = all_isoforms[gene]
        for iso in single_gene:
            if iso in female_isoforms:
                if gene in converted_female_dict:
                    converted_female_dict[gene].append(iso)
                elif gene not in converted_female_dict:
                    converted_female_dict.update({gene:[iso]})
    return converted_female_dict

def convert_male_list():
    all_isoforms = pull_all_isoforms()
    male_isoforms = pull_male_specific_isoforms()
    converted_male_dict = {}
    for gene in all_isoforms:
        single_gene = all_isoforms[gene]
        for iso in single_gene:
            if iso in male_isoforms:
                if gene in converted_male_dict:
                    converted_male_dict[gene].append(iso)
                elif gene not in converted_male_dict:
                    converted_male_dict.update({gene:[iso]})
    return converted_male_dict


#sort isoforms as shared between sexes, male specific or female specific due to expression, male or female specific due to alternative splicing
def sort_isoforms():
    converted_female_dict = convert_female_list()
    converted_male_dict = convert_male_list()
    shared_between_sexes = {}
    male_specific_expression = {}
    female_specific_expression = {}
    male_specific_splicing = {}
    female_specific_splicing = {}
    for gene in converted_female_dict:
        #examining isoforms that match shared genes between the sexes
        if gene in converted_male_dict:
            single_female_gene = converted_female_dict[gene]
            single_male_gene = converted_male_dict[gene]
            for iso in single_female_gene:
                if iso in single_male_gene:
                    #if isoform is present in both sexes = add to shared_between_sexes dictionary
                    if gene in shared_between_sexes:
                        shared_between_sexes[gene].append(iso)
                    elif gene not in shared_between_sexes:
                        shared_between_sexes.update({gene:[iso]})
                elif iso not in single_male_gene:
                    #this means the isoform is female specific due to splicing
                    if gene in female_specific_splicing:
                        female_specific_splicing[gene].append(iso)
                    elif gene not in female_specific_splicing:
                        female_specific_splicing.update({gene:[iso]})
            for i in single_male_gene:
                if i not in single_female_gene:
                    #this means the isoform is male specific due to splicing
                    if gene in male_specific_splicing:
                        male_specific_splicing[gene].append(i)
                    elif gene not in male_specific_splicing:
                        male_specific_splicing.update({gene:[i]})
        elif gene not in converted_male_dict:
            single_female_gene = converted_female_dict[gene]
            for iso2 in single_female_gene:
                if gene in female_specific_expression:
                    female_specific_expression[gene].append(iso2)
                elif gene not in female_specific_expression:
                    female_specific_expression.update({gene:[iso2]})
    #to identify isoforms that are sex specific likely due to expression
    for g in converted_male_dict:
        if g not in converted_female_dict:
            single_male_gene = converted_male_dict[g]
            for isoform in single_male_gene:
                if g in male_specific_expression:
                    male_specific_expression[g].append(isoform)
                elif g not in male_specific_expression:
                    male_specific_expression.update({g:[isoform]})
    return shared_between_sexes, male_specific_expression, female_specific_expression, male_specific_splicing, female_specific_splicing

#write isoforms to output file
#each output file will have 2 columns:
#isoform.id \t gene id \n
def write_output():
    shared_between_sexes, male_specific_expression, female_specific_expression, male_specific_splicing, female_specific_splicing = sort_isoforms()
    output_shared = sys.argv[4]
    output_male_expression = sys.argv[5]
    output_female_expression = sys.argv[6]
    output_male_splicing = sys.argv[7]
    output_female_splicing = sys.argv[8]
    with open(output_shared, 'a') as out1, open(output_male_expression, 'a') as out2, open(output_female_expression, 'a') as out3, open(output_male_splicing, 'a') as out4, open(output_female_splicing, 'a') as out5:
        for key in shared_between_sexes:
            single_key = shared_between_sexes[key]
            for v in single_key:
                final = "%s\t%s\n" % (str(key), str(v))
                out1.write(final)
        for key2 in male_specific_expression:
            single_key = male_specific_expression[key2]
            for v in single_key:
                final = "%s\t%s\n" % (str(key2), str(v))
                out2.write(final)
        for key3 in female_specific_expression:
            single_key = female_specific_expression[key3]
            for v in single_key:
                final = "%s\t%s\n" % (str(key3), str(v))
                out3.write(final)
        for key4 in male_specific_splicing:
            single_key = male_specific_splicing[key4]
            for v in single_key:
                final = "%s\t%s\n" % (str(key4), str(v))
                out4.write(final)
        for key5 in female_specific_splicing:
            single_key = female_specific_splicing[key5]
            for v in single_key:
                final = "%s\t%s\n" % (str(key5), str(v))
                out5.write(final)

write_output()

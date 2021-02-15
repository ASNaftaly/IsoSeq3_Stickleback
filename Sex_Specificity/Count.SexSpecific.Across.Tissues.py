#counting how many sex-specific spliced isoforms are present in multiple tissues
#need to read in isoform ids from sex-specific spliced isoforms and each individual tissue for that sex
#to run script: python3 Count.SexSpecific.Across.Tissues.py <all sex specific spliced isoforms (male or female)> <brain isoforms (male or female)> <liver isoforms (male or female)> <pronephros isoforms (male or female)> <gonad isoforms (male or female)> <output file>
#all input files have one column with one isoform id per line (may or may not include a header line)
#output format = 2 columns (Isoform ID \t which tissues the isoform is in)
#author: Alice Naftaly, Feb 2021

import sys

#read in all sex-specific spliced isoforms for one sex
def read_alltissues():
    iso_file = sys.argv[1]
    all_sexspecific_isoforms = []
    with open(iso_file, 'r') as all_isos:
        for line in all_isos:
            if line.startswith("PB"):
                all_sexspecific_isoforms.append(line.strip("\n"))
    return all_sexspecific_isoforms

#read brain isoforms:
def read_brain():
    brain_file = sys.argv[2]
    brain_isoforms = []
    with open(brain_file, 'r') as brain_isos:
        for line in brain_isos:
            if line.startswith("PB"):
                brain_isoforms.append(line.strip("\n"))
    return brain_isoforms

#read liver isoforms:
def read_liver():
    liver_file = sys.argv[3]
    liver_isoforms = []
    with open(liver_file, 'r') as liver_isos:
        for line in liver_isos:
            if line.startswith("PB"):
                liver_isoforms.append(line.strip("\n"))
    return liver_isoforms

#read pronephros isoforms:
def read_pronephros():
    pronephros_file = sys.argv[4]
    pronephros_isoforms = []
    with open(pronephros_file, 'r') as pronephros_isos:
        for line in pronephros_isos:
            if line.startswith("PB"):
                pronephros_isoforms.append(line.strip("\n"))
    return pronephros_isoforms

#read gonad isoforms:
def read_gonad():
    gonad_file = sys.argv[5]
    gonad_isoforms = []
    with open(gonad_file, 'r') as gonad_isos:
        for line in gonad_isos:
            if line.startswith("PB"):
                gonad_isoforms.append(line.strip("\n"))
    return gonad_isoforms

#comparing isoforms
def compare():
    all_sexspecific_isoforms = read_alltissues()
    brain_isoforms = read_brain()
    liver_isoforms = read_liver()
    pronephros_isoforms = read_pronephros()
    gonad_isoforms = read_gonad()
    sex_specific_dict = {}
    for iso in all_sexspecific_isoforms:
        if iso in brain_isoforms:
            if iso in sex_specific_dict:
                sex_specific_dict[iso].append("brain")
            elif iso not in sex_specific_dict:
                sex_specific_dict.update({iso:["brain"]})
        if iso in liver_isoforms:
            if iso in sex_specific_dict:
                sex_specific_dict[iso].append("liver")
            elif iso not in sex_specific_dict:
                sex_specific_dict.update({iso:["liver"]})
        if iso in pronephros_isoforms:
            if iso in sex_specific_dict:
                sex_specific_dict[iso].append("pronephros")
            elif iso not in sex_specific_dict:
                sex_specific_dict.update({iso:["pronephros"]})
        if iso in gonad_isoforms:
            if iso in sex_specific_dict:
                sex_specific_dict[iso].append("gonad")
            elif iso not in sex_specific_dict:
                sex_specific_dict.update({iso:["gonad"]})
    return sex_specific_dict

#writing output
def write():
    sex_specific_summary = compare()
    output = sys.argv[6]
    with open(output, 'a') as out:
        for iso in sex_specific_summary:
            single_iso = sex_specific_summary[iso]
            final_iso_format = ",".join(single_iso)
            final = "%s\t%s\n" % (str(iso), str(final_iso_format))
            out.write(final)


write()

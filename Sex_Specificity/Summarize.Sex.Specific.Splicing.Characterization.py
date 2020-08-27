#script to summarize output from Categorize.Sex.Specific.Splice.Variants.py
#will count up number of alternative splicing events and differences in TSS,TTS,UTRs
#will need to read in output from Categorize.Sex.Specific.Splice.Variants.py
#to run script: python3 Summarize.Sex.Specific.Splicing.Characterization.py <output from Categorize.Sex.Specific.Splicing.Variants.py>

import sys

#read in categorized output
#returns 2 dictionaries with counts of category type and category specifics
def read_categorized_output():
    input_file = sys.argv[1]
    diff_category_dict = {}
    category_specifics_dict = {}
    with open(input_file, 'r') as categories:
        for line in categories:
            new_line = line.split("\t")
            isoform = new_line[0]
            gene = new_line[1]
            diff_category_full = new_line[3]
            category_specifics = new_line[4].strip("\n")
            if diff_category_full in diff_category_dict:
                diff_category_dict[diff_category_full].append(isoform)
            elif diff_category_full not in diff_category_dict:
                diff_category_dict.update({diff_category_full:[isoform]})
            if category_specifics in category_specifics_dict:
                category_specifics_dict[category_specifics].append(isoform)
            elif category_specifics not in category_specifics_dict:
                category_specifics_dict.update({category_specifics:[isoform]})
    return diff_category_dict, category_specifics_dict

#count differences in Alternative Splicing vs TSS.TTS.UTR.diffs
def count_splicing_vs_alternateTS():
    diff_category_dict, category_specifics_dict = read_categorized_output()
    for key in diff_category_dict:
        print(key)
        print(len(diff_category_dict[key]))

#count differences in category specifics
def count_category_specifics():
    diff_category_dict, category_specifics_dict = read_categorized_output()
    final_dict = {}
    for key in category_specifics_dict:
        split_key_by_period = key.split(".")
        single_key = category_specifics_dict[key]
        if split_key_by_period[1].startswith("exon"):
            split_key_by_type = key.split(",")
            if len(split_key_by_type) == 1:
                if split_key_by_period[0] != '0':
                    if "alternative.splicing.only" in final_dict:
                        final_dict["alternative.splicing.only"].append(len(single_key))
                    elif "alternative.splicing.only" not in final_dict:
                        final_dict.update({"alternative.splicing.only":[len(single_key)]})
                elif split_key_by_period[0] == '0':
                    print("isoform shows up as no differences in exons or TSS/TTS")
                    print(single_key)
            elif len(split_key_by_type) == 2:
                if split_key_by_period[0] == '0':
                    if split_key_by_type[1] == "Alternate.TSS":
                        if "alternate.tss.only" in final_dict:
                            final_dict["alternate.tss.only"].append(len(single_key))
                        elif "alternate.tss.only" not in final_dict:
                            final_dict.update({"alternate.tss.only":[len(single_key)]})
                    elif split_key_by_type[1] == "Alternate.TTS":
                        if "alternate.tts.only" in final_dict:
                            final_dict["alternate.tts.only"].append(len(single_key))
                        elif "alternate.tts.only" not in final_dict:
                            final_dict.update({"alternate.tts.only":[len(single_key)]})
                elif split_key_by_period[0] != "0":
                    if split_key_by_type[1] == "Alternate.TSS":
                        if "alternative.splicing.and.alternate.tss" in final_dict:
                            final_dict["alternative.splicing.and.alternate.tss"].append(len(single_key))
                        elif "alternative.splicing.and.alternate.tss" not in final_dict:
                            final_dict.update({"alternative.splicing.and.alternate.tss":[len(single_key)]})
                    elif split_key_by_type[1] == "Alternate.TTS":
                        if "alternative.splicing.and.alternate.tts" in final_dict:
                            final_dict["alternative.splicing.and.alternate.tts"].append(len(single_key))
                        elif "alternative.splicing.and.alternate.tts" not in final_dict:
                            final_dict.update({"alternative.splicing.and.alternate.tts":[len(single_key)]})
            elif len(split_key_by_type) == 3:
                if split_key_by_period[0] == '0':
                    if "alternate.tss.and.alternate.tts" in final_dict:
                        final_dict["alternate.tss.and.alternate.tts"].append(len(single_key))
                    elif "alternate.tss.and.alternate.tts" not in final_dict:
                        final_dict.update({"alternate.tss.and.alternate.tts":[len(single_key)]})
                elif split_key_by_period[0] != "0":
                    if "alternate.splicing.and.alternate.tss.and.tts" in final_dict:
                        final_dict["alternate.splicing.and.alternate.tss.and.tts"].append(len(single_key))
                    elif "alternate.splicing.and.alternate.tss.and.tts" not in final_dict:
                        final_dict.update({"alternate.splicing.and.alternate.tss.and.tts":[len(single_key)]})
        elif split_key_by_period[0].startswith("Both"):
            split_key_by_type = key.split(",")
            if len(split_key_by_type) == 2:
                if split_key_by_type[1] == "Alternate.TSS":
                    if "alternative.splicing.and.alternate.tss" in final_dict:
                        final_dict["alternative.splicing.and.alternate.tss"].append(len(single_key))
                    elif "alternative.splicing.and.alternate.tss" not in final_dict:
                        final_dict.update({"alternative.splicing.and.alternate.tss":[len(single_key)]})
                elif split_key_by_type[1] == "Alternate.TTS":
                    if "alternative.splicing.and.alternate.tts" in final_dict:
                        final_dict["alternative.splicing.and.alternate.tts"].append(len(single_key))
                    elif "alternative.splicing.and.alternate.tts" not in final_dict:
                        final_dict.update({"alternative.splicing.and.alternate.tts":[len(single_key)]})
            elif len(split_key_by_type) == 3:
                if "alternate.splicing.and.alternate.tss.and.tts" in final_dict:
                    final_dict["alternate.splicing.and.alternate.tss.and.tts"].append(len(single_key))
                elif "alternate.splicing.and.alternate.tss.and.tts" not in final_dict:
                    final_dict.update({"alternate.splicing.and.alternate.tss.and.tts":[len(single_key)]})
        elif split_key_by_period[0].startswith("Single"):
            if key == "Single.Exon.Both.TSS.TTS.different":
                if "alternate.splicing.and.alternate.tss.and.tts" in final_dict:
                    final_dict["alternate.splicing.and.alternate.tss.and.tts"].append(len(single_key))
                elif "alternate.splicing.and.alternate.tss.and.tts" not in final_dict:
                    final_dict.update({"alternate.splicing.and.alternate.tss.and.tts":[len(single_key)]})
            elif key == "Single.Exon.Different.TSS":
                if "alternate.tss.only" in final_dict:
                    final_dict["alternate.tss.only"].append(len(single_key))
                elif "alternate.tss.only" not in final_dict:
                    final_dict.update({"alternate.tss.only":[len(single_key)]})
            elif key == "Single.Exon.Different.TTS":
                if "alternate.tts.only" in final_dict:
                    final_dict["alternate.tts.only"].append(len(single_key))
                elif "alternate.tts.only" not in final_dict:
                    final_dict.update({"alternate.tts.only":[len(single_key)]})
        elif key.startswith("Alternate"):
            split_key_by_type = key.split(",")
            if len(split_key_by_type) == 1:
                if split_key_by_type == ["Alternate.TSS"]:
                    if "alternate.tss.only" in final_dict:
                        final_dict["alternate.tss.only"].append(len(single_key))
                    elif "alternate.tss.only" not in final_dict:
                        final_dict.update({"alternate.tss.only":[len(single_key)]})
                elif split_key_by_type == ["Alternate.TTS"]:
                    if "alternate.tts.only" in final_dict:
                        final_dict["alternate.tts.only"].append(len(single_key))
                    elif "alternate.tts.only" not in final_dict:
                        final_dict.update({"alternate.tts.only":[len(single_key)]})
            elif len(split_key_by_type) == 2:
                if "alternate.tss.and.alternate.tts" in final_dict:
                    final_dict["alternate.tss.and.alternate.tts"].append(len(single_key))
                elif "alternate.tss.and.alternate.tts" not in final_dict:
                    final_dict.update({"alternate.tss.and.alternate.tts":[len(single_key)]})
        elif key.startswith("No"):
            if "alternative.splicing.only" in final_dict:
                final_dict["alternative.splicing.only"].append(len(single_key))
            elif "alternative.splicing.only" not in final_dict:
                final_dict.update({"alternative.splicing.only":[len(single_key)]})
        else:
            print(key)
    for type in final_dict:
        single_type = final_dict[type]
        total = sum(single_type)
        print(type)
        print(total)


#call functions()
def call():
    categories_count = count_splicing_vs_alternateTS()
    categories_specifics = count_category_specifics()

call()

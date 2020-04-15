#classifing where ncRNAs are located (while also taking into account annotated vs novel transcripts)
#input will be classification file and novel/annotated transcripts (isoform ids)
#to run script: python3 Locating_ncRNAs.py <noncoding isoforms classification file> <annotated transcripts file from Annotated_vs_Novel_ncRNAS.py> <novel transcripts file from Annotated_vs_Novel_ncRNAS.py> <output breakdown of SQANTI splice type for annotated and novel transcripts>
#Author: Alice Naftaly, April 2020

import sys

#read classification file
#returns dictionary with key = isoform and value = splice type
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                splice_type = new_line[5]
                class_dict.update({isoform:splice_type})
    return class_dict


#read in annotated transcripts file
#returns list of isoform ids for all annotated transcripts
def read_annotated_transcripts():
    input_file = sys.argv[2]
    annotated_isoforms = []
    with open(input_file, 'r') as transcripts:
        for line in transcripts:
            final = line.strip("\n")
            annotated_isoforms.append(final)
    return annotated_isoforms


#read in novel transcripts file
#returns list of isoform ids for all novel transcripts
def read_novel_transcripts():
    input_file = sys.argv[3]
    novel_isoforms = []
    with open(input_file, 'r') as transcripts:
        for line in transcripts:
            final = line.strip("\n")
            novel_isoforms.append(final)
    return novel_isoforms


#classifying annotated transcripts by location
#returns dictionary with key == splice type and value == list of isoform ids for each splice type
def classify_annotated_transcripts():
    annotated_transcripts = read_annotated_transcripts()
    class_dict = read_class()
    annotated_transcripts_splice_dict = {}
    for isoform in annotated_transcripts:
        splice_type = class_dict[isoform]
        if splice_type in annotated_transcripts_splice_dict:
            annotated_transcripts_splice_dict[splice_type].append(isoform)
        elif splice_type not in annotated_transcripts_splice_dict:
            annotated_transcripts_splice_dict.update({splice_type:[isoform]})
    return annotated_transcripts_splice_dict


def classify_novel_transcripts():
    novel_transcripts = read_novel_transcripts()
    class_dict = read_class()
    novel_transcripts_splice_dict = {}
    for isoform in novel_transcripts:
        splice_type = class_dict[isoform]
        if splice_type in novel_transcripts_splice_dict:
            novel_transcripts_splice_dict[splice_type].append(isoform)
        elif splice_type not in novel_transcripts_splice_dict:
            novel_transcripts_splice_dict.update({splice_type:[isoform]})
    return novel_transcripts_splice_dict


#print summary counts for annotated and novel transcripts:
def summary_counts():
    classified_annotated_transcripts = classify_annotated_transcripts()
    classified_novel_transcripts = classify_novel_transcripts()
    print("Annotated Transcripts\n")
    for splice_type in classified_annotated_transcripts:
        print(splice_type)
        print(len(classified_annotated_transcripts[splice_type]))
    print("Novel Transcripts\n")
    for splice_type_2 in classified_novel_transcripts:
        print(splice_type_2)
        print(len(classified_novel_transcripts[splice_type_2]))

#write annotated and novel transcripts to output file (1 file)
#header = Isoform.ID \t Annotated/Novel \t Splice.Type \n
def write_classification():
    classified_annotated_transcripts = classify_annotated_transcripts()
    classified_novel_transcripts = classify_novel_transcripts()
    output = sys.argv[4]
    with open(output, 'a') as out:
        header = "Isoform.ID\tAnnotated.Or.Novel.Transcript\tSplice.Type\n"
        out.write(header)
        for splice in classified_annotated_transcripts:
            single_splice = classified_annotated_transcripts[splice]
            for isoform in single_splice:
                final = "%s\t%s\t%s\n" % (str(isoform), "Annotated",str(splice))
                out.write(final)
        for splice_2 in classified_novel_transcripts:
            single_splice = classified_novel_transcripts[splice_2]
            for isoform_2 in single_splice:
                final = "%s\t%s\t%s\n" % (str(isoform_2), "Annotated",str(splice_2))
                out.write(final)

#call all functions
def call():
    summary = summary_counts()
    write = write_classification()

call()

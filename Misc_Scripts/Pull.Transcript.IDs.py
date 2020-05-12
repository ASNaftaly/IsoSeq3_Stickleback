#Need transcript ids and isoform ids for Isoforms
#if Ensembl transcript is annotated, need that ENSGACT number otherwise can use the isoform id
#to run script: python3 Pull.Transcript.IDs.py <classification file> <output file, with 1 value per line either ENSGACT or isoform id number>

import sys


#read in classification file
#returns list of transcripts
def read_class():
    class_file = sys.argv[1]
    class_dict = {}
    final_transcripts = []
    with open(class_file, 'r') as class_info:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split()
                isoform = new_line[0]
                transcript_id = new_line[7]
                if transcript_id.startswith("ENSGACT"):
                    final_transcripts.append(transcript_id)
                elif transcript_id == "novel":
                    final_transcripts.append(isoform)
    return final_transcripts


#write output file
def write():
    transcripts = read_class()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for transcript in transcripts:
            final = "%s\n" % str(transcript)
            out.write(final)

write()

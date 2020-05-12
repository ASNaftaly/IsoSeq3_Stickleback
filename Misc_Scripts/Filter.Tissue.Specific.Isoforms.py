#pull transcript/isoform ids from tissue specific isoform file
#tissue specific isoform file has the following format: Isoform.ID \t Gene.ID \t Transcript.ID \t Exon.Number \t Splice.Type \t Sub.Type
#need to pull transcript ID if ENSGACT or isoform ID if transcript is Novel
#to run script: python3 Filter.Tissue.Specific.Isoforms.py <tissue specific isoforms file> <output file with transcript or isoform id>
#Author: Alice Naftaly, May 2020

import sys


#read in tissue specific file
#returns list of transcript and isoform ids
def read_tissue_specific_isoforms():
    tissue_file = sys.argv[1]
    final_transcripts = []
    with open(tissue_file, 'r') as tissue_isoforms:
        for line in tissue_isoforms:
            if line.startswith("PB"):
                new_line = line.split()
                isoform_id = new_line[0]
                transcript_id = new_line[2]
                if transcript_id.startswith("ENSGACT"):
                    final_transcripts.append(transcript_id)
                elif transcript_id == "novel":
                    final_transcripts.append(isoform_id)
    return final_transcripts


#write transcripts to file
def write():
    final_transcripts = read_tissue_specific_isoforms()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for transcript in final_transcripts:
            final = "%s\n" % str(transcript)
            out.write(final)

write()

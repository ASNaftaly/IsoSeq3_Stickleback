#checking all start and stop sequences from 17226 isoforms that had start and stop codons
#to run script: python3 Check.start.stop.seqs.py <seqs file as output from PullSeq_Check.py <output1> <output2> <output3> <output4> <output5> <output6> <output7> <output8> <output9>
#author: Alice Naftaly, Aug 2020

import sys

#reads file with isoform ids and start and stop seqs
#returns list where order is isoform, start codon seq, stop codon seq
def read_seqs():
    seqs_file = sys.argv[1]
    isoform_dict = {}
    with open(seqs_file, 'r') as seqs:
        for line in seqs:
            new_line = line.strip()
            if new_line.startswith("PB"):
                isoform = new_line
            else:
                if isoform in isoform_dict:
                    isoform_dict[isoform].append(new_line)
                elif isoform not in isoform_dict:
                    isoform_dict.update({isoform:[new_line]})
    return isoform_dict

#checking sequences
#returns several lists of isoforms
#lists include:
#isoforms that no sequnces were pulled *should recheck
#isoforms with the correct sequences (ATG and TAA,TAG, or TGA or CAT and TTA,CTA, or TCA)
#isoforms with ATG start, but with stop TTA,CTA, or TCA
#isoforms with CAT start, but with stop TAA, TAG, TGA
#isoforms with the correct start seq, but wrong stop sequence (not TTA, CTA, TCA, TAA, TAG, or TGA)
#isoforms with no start sequences
#isoforms with completely wrong start and stop sequences *should recheck these
def check_seqs():
    isoform_seqs = read_seqs()
    isoforms_with_no_seqs = []
    isoforms_with_correct_start_stop_seqs = []
    isoforms_with_ATG_start_minusstrand_stop = []
    isoforms_with_CAT_start_plusstrand_stop = []
    isoforms_with_plusstrand_stop_wrong_start = []
    isoforms_with_minusstrand_stop_wrong_start = []
    isoforms_with_correct_start_wrong_stop_seq = []
    no_start_seq = []
    wrong_seqs = []
    for iso in isoform_seqs:
        single_isoform = isoform_seqs[iso]
        if len(single_isoform) == 2:
            isoforms_with_no_seqs.append(iso)
        elif len(single_isoform) == 4:
            start_seq = single_isoform[1]
            start_seq_strip = start_seq.strip("[]")
            start_seq_split = start_seq_strip.split(",")
            new_start_seq = []
            for val in start_seq_split:
                val_strip = val.strip("'")
                val_strip_2 = val_strip.strip(" ")
                final_val_strip = val_strip_2.strip("'")
                new_start_seq.append(final_val_strip)
            stop_seq = single_isoform[3]
            stop_seq_strip = stop_seq.strip("[]")
            stop_seq_split = stop_seq_strip.split(",")
            new_stop_seq = []
            for v in stop_seq_split:
                v_strip = v.strip("'")
                v_strip_2 = v_strip.strip(" ")
                final_v_strip = v_strip_2.strip("'")
                new_stop_seq.append(final_v_strip)
            joined_start = "".join(new_start_seq)
            joined_stop = "".join(new_stop_seq)
            if joined_start == "ATG":
                if joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_correct_start_stop_seqs.append(iso)
                elif joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_ATG_start_minusstrand_stop.append(iso)
                else:
                    isoforms_with_correct_start_wrong_stop_seq.append(iso)
            elif joined_start == "CAT":
                if joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_correct_start_stop_seqs.append(iso)
                elif joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_CAT_start_plusstrand_stop.append(iso)
                else:
                    isoforms_with_correct_start_wrong_stop_seq.append(iso)
            else:
                if joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_plusstrand_stop_wrong_start.append(iso)
                elif joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_minusstrand_stop_wrong_start.append(iso)
                else:
                    wrong_seqs.append(iso)
        elif len(single_isoform) == 6:
            start_seq = single_isoform[3]
            start_seq_strip = start_seq.strip("[]")
            start_seq_split = start_seq_strip.split(",")
            new_start_seq = []
            for val in start_seq_split:
                val_strip = val.strip("'")
                val_strip_2 = val_strip.strip(" ")
                final_val_strip = val_strip_2.strip("'")
                new_start_seq.append(final_val_strip)
            stop_seq = single_isoform[5]
            stop_seq_strip = stop_seq.strip("[]")
            stop_seq_split = stop_seq_strip.split(",")
            new_stop_seq = []
            for v in stop_seq_split:
                v_strip = v.strip("'")
                v_strip_2 = v_strip.strip(" ")
                final_v_strip = v_strip_2.strip("'")
                new_stop_seq.append(final_v_strip)
            joined_start = "".join(new_start_seq)
            joined_stop = "".join(new_stop_seq)
            if joined_start == "ATG":
                if joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_correct_start_stop_seqs.append(iso)
                elif joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_ATG_start_minusstrand_stop.append(iso)
                else:
                    isoforms_with_correct_start_wrong_stop_seq.append(iso)
            elif joined_start == "CAT":
                if joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_correct_start_stop_seqs.append(iso)
                elif joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_CAT_start_plusstrand_stop.append(iso)
                else:
                    isoforms_with_correct_start_wrong_stop_seq.append(iso)
            else:
                if joined_stop == "TAA" or joined_stop == "TAG" or joined_stop == "TGA":
                    isoforms_with_plusstrand_stop_wrong_start.append(iso)
                elif joined_stop == "TTA" or joined_stop == "CTA" or joined_stop == "TCA":
                    isoforms_with_minusstrand_stop_wrong_start.append(iso)
                else:
                    wrong_seqs.append(iso)
        elif len(single_isoform) == 3:
            no_start_seq.append(iso)
        elif len(single_isoform) == 7:
            isoforms_with_minusstrand_stop_wrong_start.append(iso)
    return isoforms_with_no_seqs, isoforms_with_correct_start_stop_seqs, isoforms_with_ATG_start_minusstrand_stop, isoforms_with_CAT_start_plusstrand_stop, isoforms_with_plusstrand_stop_wrong_start, isoforms_with_minusstrand_stop_wrong_start, isoforms_with_correct_start_wrong_stop_seq, no_start_seq, wrong_seqs


#write output
#will be 9 output files
def write():
    isoforms_with_no_seqs, isoforms_with_correct_start_stop_seqs, isoforms_with_ATG_start_minusstrand_stop, isoforms_with_CAT_start_plusstrand_stop, isoforms_with_plusstrand_stop_wrong_start, isoforms_with_minusstrand_stop_wrong_start, isoforms_with_correct_start_wrong_stop_seq, no_start_seq, wrong_seqs = check_seqs()
    output1 = sys.argv[2]
    output2 = sys.argv[3]
    output3 = sys.argv[4]
    output4 = sys.argv[5]
    output5 = sys.argv[6]
    output6 = sys.argv[7]
    output7 = sys.argv[8]
    output8 = sys.argv[9]
    output9 = sys.argv[10]
    with open(output1, 'a') as out1, open(output2, 'a') as out2, open(output3, 'a') as out3, open(output4, 'a') as out4, open(output5, 'a') as out5, open(output6, 'a') as out6, open(output7, 'a') as out7, open(output8, 'a') as out8, open(output9, 'a') as out9:
        for iso1 in isoforms_with_no_seqs:
            final = "%s\n" % str(iso1)
            out1.write(final)
        for iso2 in isoforms_with_correct_start_stop_seqs:
            final = "%s\n" % str(iso2)
            out2.write(final)
        for iso3 in isoforms_with_ATG_start_minusstrand_stop:
            final = "%s\n" % str(iso3)
            out3.write(final)
        for iso4 in isoforms_with_CAT_start_plusstrand_stop:
            final = "%s\n" % str(iso4)
            out4.write(final)
        for iso5 in isoforms_with_plusstrand_stop_wrong_start:
            final = "%s\n" % str(iso5)
            out5.write(final)
        for iso6 in isoforms_with_minusstrand_stop_wrong_start:
            final = "%s\n" % str(iso6)
            out6.write(final)
        for iso7 in isoforms_with_correct_start_wrong_stop_seq:
            final = "%s\n" % str(iso7)
            out7.write(final)
        for iso8 in no_start_seq:
            final = "%s\n" % str(iso8)
            out8.write(final)
        for iso9 in wrong_seqs:
            final = "%s\n" % str(iso9)
            out9.write(final)


write()

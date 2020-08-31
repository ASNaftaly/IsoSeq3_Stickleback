#Count coverage of splice junctions
#need to read in both exon positions and splice junctions coverage for each tissue
#will combine all splice junctions coverage to confirm full transcriptome
#to run script: python3 Count.SJ.Coverage.All.Isoforms.py <exon positions gtf> <splice junctions combined file> <output summary file>
#Author: Alice Naftaly, June 2020

import sys

#read in gtf file with exon positions
#returns dictionary with key == isoform id and value == [chr num, strand, exon start, exon end]
def read_isoseq_gtf():
    gtf_file = sys.argv[1]
    exon_dict = {}
    final_exon_dict = {}
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            new_line = line.split()
            chr_num = new_line[0]
            exon_start = new_line[3]
            exon_end = new_line[4]
            strand = new_line[6]
            isoform_id_full = new_line[9].strip(";")
            isoform_id = isoform_id_full.strip("\"")
            dict_value = [chr_num, strand, exon_start, exon_end]
            if isoform_id in exon_dict:
                exon_dict[isoform_id].append(dict_value)
            elif isoform_id not in exon_dict:
                exon_dict.update({isoform_id:[dict_value]})
    for isoform in exon_dict:
        single_isoform = exon_dict[isoform]
        if len(single_isoform) == 1:
            single = single_isoform[0]
            final_exon_dict.update({isoform:[single]})
        elif len(single_isoform) > 1:
            strand = single_isoform[0][1]
            if strand == "+":
                final_exon_dict.update({isoform:single_isoform})
            elif strand == "-":
                single_isoform.reverse()
                final_exon_dict.update({isoform:single_isoform})
    return final_exon_dict


#get intron positions for all isoforms with more than 1 exon
#returns dictionary with key == isoform and value == [chr num, strand, intron pos (farthest upstream; smaller number), intron pos (farthest downstream; larger number)]
#intron start and end will include 1bp of exon
def get_intron_positions():
    exon_positions = read_isoseq_gtf()
    introns_dict = {}
    for isoform in exon_positions:
        single_isoform = exon_positions[isoform]
        if len(single_isoform) == 1:
            continue
        elif len(single_isoform) > 1:
            chr_num = single_isoform[0][0]
            strand = single_isoform[0][1]
            exon_count = 0
            while exon_count < len(single_isoform)-1:
                exon_1 = single_isoform[exon_count]
                exon_2 = single_isoform[exon_count+1]
                if strand == "+":
                    intron_start = int(exon_1[3])+1
                    intron_end = int(exon_2[2])-1
                    new_dict_value = [chr_num, strand, intron_start, intron_end]
                elif strand == "-":
                    intron_start = int(exon_1[2])-1
                    intron_end = int(exon_2[3])+1
                    new_dict_value = [chr_num, strand, intron_end, intron_start]
                if isoform in introns_dict:
                    introns_dict[isoform].append(new_dict_value)
                elif isoform not in introns_dict:
                    introns_dict.update({isoform:[new_dict_value]})
                exon_count += 1
    return introns_dict

#read in splice junction Coverage
#returns dictionary with key == chr number and value == [intron start, intron end, number of uniquely mapped reads, number of multimapped reads] *everything is increasing in position number
def read_sj_coverage():
    sj_file = sys.argv[2]
    sj_coverage = {}
    with open(sj_file, 'r') as splice_junctions:
        for line in splice_junctions:
            new_line = line.split()
            chr_num = new_line[0]
            intron_start = int(new_line[1])
            intron_end = int(new_line[2])
            num_unique_mapped_reads = int(new_line[6])
            num_multi_mapped_reads = int(new_line[7])
            dict_value = [intron_start, intron_end, num_unique_mapped_reads, num_multi_mapped_reads]
            if chr_num in sj_coverage:
                sj_coverage[chr_num].append(dict_value)
            elif chr_num not in sj_coverage:
                sj_coverage.update({chr_num:[dict_value]})
    return sj_coverage


#get read coverage for all junctions
def combine_sj_coverage_allsj():
    all_junctions = get_intron_positions()
    sj_coverage = read_sj_coverage()
    junction_coverage = {}
    count = 0
    for chr in sj_coverage:
        single_chr = sj_coverage[chr]
        for sj in single_chr:
            sj_start = sj[0]
            sj_end = sj[1]
            unique_reads = sj[2]
            multi_reads = sj[3]
            for isoform in all_junctions:
                single_isoform = all_junctions[isoform]
                for intron in single_isoform:
                    intron_start = intron[2]
                    intron_end = intron[3]
                    if intron_start == sj_start and intron_end == sj_end:
                        dict_value = [intron_start, intron_end, unique_reads, multi_reads]
                        if isoform in junction_coverage:
                            junction_coverage[isoform].append(dict_value)
                        elif isoform not in junction_coverage:
                            junction_coverage.update({isoform:[dict_value]})
    final_junction_coverage = {}
    for iso in junction_coverage:
        single_iso = junction_coverage[iso]
        single_iso_dict = {}
        for val in single_iso:
            iso_val = [val[2], val[3]]
            if str(val[0]) in single_iso_dict:
                single_iso_dict[str(val[0])].append(iso_val)
            elif str(val[0]) not in single_iso_dict:
                single_iso_dict.update({str(val[0]):[iso_val]})
        for read_count in single_iso_dict:
            unique_mapped = 0
            multi_mapped = 0
            single_read_count = single_iso_dict[read_count]
            for v in single_read_count:
                unique_mapped += v[0]
                multi_mapped += v[1]
            final_dict_value = [unique_mapped, multi_mapped]
            if iso in final_junction_coverage:
                final_junction_coverage[iso].append(final_dict_value)
            elif iso not in final_junction_coverage:
                final_junction_coverage.update({iso:[final_dict_value]})
    return final_junction_coverage

#create summary file
#returns dictionary with key == isoform id and value == [isoform, number of exons, number of introns, number of introns with coverage, number of introns without coverage, list of introns without coverage, average unique reads coverage of introns, average multi mapped reads coverage of introns, minimum unique reads coverage, minimum multi mapped reads coverage, maximum unique reads coverage, maximum multimapped reads coverage]
def create_summary():
    all_junctions = read_isoseq_gtf()
    intron_coverage = combine_sj_coverage_allsj()
    summary_dict = {}
    for isoform in all_junctions:
        single_isoform = all_junctions[isoform]
        if isoform in intron_coverage:
            single_isoform_coverage = intron_coverage[isoform]
            number_exons = len(single_isoform)
            number_introns = len(single_isoform)-1
            number_introns_covered = len(single_isoform_coverage)
            introns_with_coverage = 0
            introns_without_coverage = 0
            total_unique_coverage = []
            total_multi_coverage = []
            covered_introns = []
            if number_introns == number_introns_covered:
                for intron in single_isoform_coverage:
                    intron_unique_reads = intron[0]
                    intron_multi_reads = intron[0]
                    total_unique_coverage.append(intron_unique_reads)
                    total_multi_coverage.append(intron_multi_reads)
                    if intron_unique_reads == 0 and intron_multi_reads == 0:
                        introns_without_coverage += 1
                    else:
                        introns_with_coverage += 1
                average_unique_reads = round(sum(total_unique_coverage)/len(total_unique_coverage),2)
                average_mutli_reads = round(sum(total_multi_coverage)/len(total_multi_coverage),2)
                min_unique_reads = min(total_unique_coverage)
                min_multi_reads = min(total_multi_coverage)
                max_unique_reads = max(total_unique_coverage)
                max_multi_reads = max(total_multi_coverage)
                dict_value = [isoform, str(number_exons), str(number_introns), str(introns_without_coverage),  str(average_unique_reads), str(average_mutli_reads), str(min_unique_reads), str(min_multi_reads), str(max_unique_reads), str(max_multi_reads)]
                summary_dict.update({isoform:dict_value})
                #want to know which introns are not covered
            elif number_introns != number_introns_covered:
                num_introns_missing = number_introns - number_introns_covered
                introns_without_coverage += num_introns_missing
                for intron in single_isoform_coverage:
                    covered_introns.append(intron[0])
                    intron_unique_reads = intron[0]
                    intron_multi_reads = intron[1]
                    total_unique_coverage.append(intron_unique_reads)
                    total_multi_coverage.append(intron_multi_reads)
                    if intron_unique_reads == 0 and intron_multi_reads == 0:
                        introns_without_coverage += 1
                    else:
                        introns_with_coverage += 1
                average_unique_reads = round(sum(total_unique_coverage)/len(total_unique_coverage),2)
                average_mutli_reads = round(sum(total_multi_coverage)/len(total_multi_coverage),2)
                min_unique_reads = min(total_unique_coverage)
                min_multi_reads = min(total_multi_coverage)
                max_unique_reads = max(total_unique_coverage)
                max_multi_reads = max(total_multi_coverage)
                dict_value = [isoform, str(number_exons), str(number_introns), str(introns_without_coverage), str(average_unique_reads), str(average_mutli_reads), str(min_unique_reads), str(min_multi_reads), str(max_unique_reads), str(max_multi_reads)]
                summary_dict.update({isoform:dict_value})
        elif isoform not in intron_coverage:
            number_exons = len(single_isoform)
            number_introns = number_exons - 1
            dict_value = [isoform, str(number_exons), str(number_introns),".", ".", ".", ".", ".", ".", "."]
            summary_dict.update({isoform:dict_value})
    return summary_dict

#write output file
def write():
    summary_dict = create_summary()
    output = sys.argv[3]
    with open(output, 'a') as out:
        header = "Isoform.ID\tNo.of.Exons\tNo.of.Introns\tNo.Introns.without.Cov\tAve.Unique.ReadCov\tAve.MultiMapped.Read.Cov\tMin.Unique.Read.Cov\tMin.MultiMapped.Read.Cov\tMax.Unique.Read.Cov\tMax.MultiMapped.Read.Cov\n"
        out.write(header)
        for isoform in summary_dict:
            single_isoform = summary_dict[isoform]
            final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(single_isoform[0]),str(single_isoform[1]),str(single_isoform[2]),str(single_isoform[3]), str(single_isoform[4]), str(single_isoform[5]), str(single_isoform[6]), str(single_isoform[7]),str(single_isoform[8]), str(single_isoform[9]))
            out.write(final)

write()

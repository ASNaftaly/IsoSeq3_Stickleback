#pulling high confidence junctions from all junctions
#basically if the second column of the confidence file is Y and the third column is N; that junction has high confidence in the PSI levels
#to run script: python3 Pull.High.Confidence.Junctions.py <confidence file> <all junctions PSI values file> <output file>
#Author: Alice Naftaly, June 2020

import sys

#reads in confidence file with format: Junction id \t Sufficient Read Coverage (Y or N) \t Imbalance in read counts (Y or N)
#returns dictionary with key == junction id and value == [ sufficient read coverage, imbalance in read counts]
def read_confidence_file():
    confidence_file = sys.argv[1]
    confidence_dict = {}
    with open(confidence_file, 'r') as confidence_levels:
        for line in confidence_levels:
            if line.startswith("Junction"):
                continue
            else:
                new_line = line.split("\t")
                junction_id = new_line[0]
                sufficient_read_cov = new_line[1]
                imbalance = new_line[2].strip("\n")
                dict_value = [sufficient_read_cov, imbalance]
                confidence_dict.update({junction_id:dict_value})
    return confidence_dict

#read all Junctions PSI values with format: junction id \t PSI value
#returns PSI dictionary with key == junction id and value == PSI value
def read_PSI_values():
    psi_file = sys.argv[2]
    psi_dict = {}
    with open(psi_file, 'r') as psi_calcs:
        for line in psi_calcs:
            new_line = line.split("\t")
            junction_id = new_line[0]
            psi_value = new_line[1].strip("\n")
            psi_dict.update({junction_id:psi_value})
    return psi_dict


#pull high confidence junctions:
def pull_high_confidence_junctions():
    psi_values = read_PSI_values()
    confidence_levels = read_confidence_file()
    confident_psi_values = {}
    for junction in confidence_levels:
        if junction in psi_values:
            single_confidence = confidence_levels[junction]
            read_cov = single_confidence[0]
            imbalance = single_confidence[1]
            single_psi_value = psi_values[junction]
            if read_cov == "Y" and imbalance == "N":
                confident_psi_values.update({junction:single_psi_value})
    return confident_psi_values


#write confident psi values to new file in same format as original psi values file
def write():
    confident_psi = pull_high_confidence_junctions()
    output = sys.argv[3]
    with open(output, 'a') as out:
        for junction in confident_psi:
            final = "%s\t%s\n" % (str(junction), str(confident_psi[junction]))
            out.write(final)


write()

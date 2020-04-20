#Pooling counts from InterProScan databases
#basically, for each identifier in a database, I want to count up the number of occurrences and write this as a new output file to be read by R for plotting
#uses 2nd output file from Summarize_All_InterProScan_Results.py
#to run script: python3 Interproscan_database_counts.py <for R file from Summarize_All_InterProScan_Results.py> <output file>
#Author: Alice Naftaly, April 2020

import sys

#read forR input file
#returns dictionary with key == database name and value == accession number/identifier (will be different for each program)
def read_input():
    input_file = sys.argv[1]
    interproscan_dict = {}
    with open(input_file, 'r') as info:
        for line in info:
            new_line = line.split("\t")
            database = new_line[0]
            identifier = new_line[2]
            if database in interproscan_dict:
                interproscan_dict[database].append(identifier)
            elif database not in interproscan_dict:
                interproscan_dict.update({database:[identifier]})
    return interproscan_dict


#counting occurrences:
def write_database_counts():
    interproscan_dict = read_input()
    output = sys.argv[2]
    with open(output, 'a') as out:
        for database in interproscan_dict:
            single_database = interproscan_dict[database]
            counts_dict = {}
            for entry in single_database:
                if entry in counts_dict:
                    counts_dict[entry].append("1")
                elif entry not in counts_dict:
                    counts_dict.update({entry:["1"]})
            print("Number of unique identifiers for database")
            print(database, len(counts_dict))
            for key in counts_dict:
                key_counts = len(counts_dict[key])
                final = "%s\t%s\t%s\n" % (str(database), str(key), str(key_counts))
                out.write(final)


write_database_counts()

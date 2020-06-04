#converting classification file into an R readable format with the following columns
#isoform id, structural_category, associated gene, associated transcript, isoform expression, gene expression, ratio expression, coding potential
#to run script: python3 Convert.Class.to.R.py <classification file with short read expression> <output>

import sys

#read and convert classification file
def convert_class():
    class_file = sys.argv[1]
    output = sys.argv[2]
    class_dict = {}
    with open(class_file, 'r') as class_info, open(output, 'a') as out:
        for line in class_info:
            if line.startswith("PB"):
                new_line = line.split("\t")
                isoform = new_line[0]
                structural_category = new_line[5]
                gene = new_line[6]
                transcript = new_line[7]
                isoform_expression = new_line[23]
                gene_expression = new_line[24]
                ratio_isoform_to_gene_expression = new_line[25]
                coding_potential = new_line[27]
                final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(isoform), str(structural_category), str(gene), str(transcript), str(isoform_expression), str(gene_expression), str(ratio_isoform_to_gene_expression), str(coding_potential))
                out.write(final)

convert_class()

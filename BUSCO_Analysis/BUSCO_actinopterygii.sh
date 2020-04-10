#program: BUSCO 3.0.2

#command:
#need to remove <> prior to running
#downloaded database to local machine from http://busco.ezlab.org/ 

python3 run_BUSCO.py -i <input faa file> -o <output file> -l <actinopterygii_odb9> -m proteins -c <number of cores to use> -sp zebrafish

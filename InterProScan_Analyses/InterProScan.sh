#program: InterProScan v5.32

#command:
#need to remove <> prior to running script
#protein sequences were used as input

sh interproscan.sh -b <output file> -cpu <number of cpus to use> -i <input faa file> -t p -T <temporary directory for files>

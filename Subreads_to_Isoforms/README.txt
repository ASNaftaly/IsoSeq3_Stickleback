#scripts to go from raw subreads to predicted transcripts using Isoseq3 v3.1 pipeline
#completed through a combination of shell/bash commands and python scripts

#general pipeline order:
#1. build CCS reads from subreads files
#1a. may need to create xml file if multiple subreads files before running 1
#2. Classify full length reads
#2a. create xml from multiple CCS bam files (only needed if more than one CCS file is used)
#2b. remove primers using Lima
#2c. classify full length reads using refine from Isoseq3
#3. Cluster full length reads using Cluster from Isoseq3
#4. Full length reads were polished using polish from isoseq3
#5. Map full length polished reads to genome using minimap2
#5a. sort sam to bam using Samtools
#6. Collapse identical isoforms using cDNA CupCake
#7. get count information using cDNA CupCake - needed for SQANTI
#8. Run SQANTI
#8a. run SQANTI classification = classifies isoforms
#8b. run SQANTI qc = filters out artifacts through machine learning
#8c. remove artifacts in 8a files using FilterSQANTI_input.py
#8d. rerun SQANTI classification = final classification files


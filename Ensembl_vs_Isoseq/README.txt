#README for Ensembl vs Isoseq
#scripts to compare various aspects of Ensembl annotations and Isoseq annotations

EnsemblvsIsoseq_summary.txt
  this script counts the number of isoforms per gene for both ensembl and isoseq annotations
    splits these into 1 isoform, 2-3 isoforms, 4-5 isoforms, and more than 6 isoforms
  prints out the number of coding vs noncoding isoforms and structural categories

  input files needed: classification files for ensembl and isoseq
  output files: prints to standard output


FSM_TSS_TSS_differences.py
  this script pulls the full splice matches from the isoseq annotations and compares the TSSs and TTSs to ensembl annotations to determine how different the positions are

  input files needed: isoseq classification file, exon isoseq gtf file
  output files: output files with TSS differences and TTS differences


NovelIsoform_Comparison_to_Ensembl.py
  not completed and not used in paper

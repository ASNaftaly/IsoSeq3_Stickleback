#program: picard 2.16

#command:
#will need to remove <> prior to running script

java -jar picard.jar DownsampleSam I=<input bam file> O=<output bam file> P=<% of reads desired> VALIDATION_STRINGENCY=SILENT

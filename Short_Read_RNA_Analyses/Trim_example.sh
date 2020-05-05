#program: Trimmomatic/0.36


for file in *_R1_*
  do
    export output=${file%%_R*}
    java -jar trimmomatic-0.36.jar PE -threads 16 -phred33 -trimlog *.log -basein ${file} -baseout ${output}_trim.fastq.gz TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:50
done

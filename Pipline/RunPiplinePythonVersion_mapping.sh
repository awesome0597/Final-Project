#!/usr/bin/env tcsh

set path1="Insert pathe to raw data"
set path2="Insert path to script location"
set type="Insert type of data"
set snrna="Insert fata file or path to fasta file"
set file1="Insert pathe to genome file or genome file"

####################################################################################

set i = 1
foreach fastq_file ($path1/*PRS*_R1.fastq.gz)
  echo $fastq_file
  set base=`basename $fastq_file | cut -f1,2 -d\_`
  echo $base
  set fastq_file2=`echo $fastq_file | sed -e 's/_R1/_R2/'`
  echo $fastq_file2

  smalt map $snrna $fastq_file $fastq_file2 > $base\_vs_$type$i.sam

  samtools view -f 0x02 -bS $base\_vs_$type$i.sam > $base\_vs_$type$i\_good_pairs.bam
  rm $base\_vs_$type$i.sam

  set i = `expr $i + 1`
end
###################################################################################



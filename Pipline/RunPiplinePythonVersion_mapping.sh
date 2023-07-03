#!/usr/bin/env tcsh

set path1="/home/user/brucea/PRS_Data/"
set path2="/home/user/brucea/Scripts/"
set type="smallRNA"
set snrna="/home/user/brucea/PRS_Data/TB_small_RNAs_DB_w_praveen.fa"
set file1="/home/user/tirza/DB/TB_small_RNAs_DB_w_praveen.genome"

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



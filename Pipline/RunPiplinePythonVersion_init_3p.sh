#!/usr/bin/env tcsh

set path1="Insert pathe to raw data"
set path2="Insert path to script location"
set type="Insert type of data"
set snrna="Insert fata file or path to fasta file"
set file1="Insert pathe to genome file or genome file"

###################################################################################
foreach bam_file (*.bam)
    echo $bam_file
    set prefix=`basename $bam_file|cut -f1,2 -d\_`
    ##BAM_TO_BED
    bedtools bamtobed -bedpe -i  $prefix\_vs_$type\_good_pairs.bam  | awk '{print $1 "\t" $2  "\t" $6 "\t" $7 "\t" $8 "\t" $9}' | sort -k1,1 -k2,2n > $prefix\_vs_$type\_good_pairs.sorted.bed
    ##GenomeCoverage
    bedtools genomecov -d -g $file1 -i $prefix\_vs_$type\_good_pairs.sorted.bed  >  $prefix\_vs_$type\_good_pairs.sorted.genomecov

    ##init
    python3 $path2/CountInitiating_prototype.py $file1 $prefix\_vs_smallRNA_good_pairs.sorted.bed > $prefix\_vs_smallRNA_good_pairs_prototype.sorted.init

    python3 $path2/count3p_prototype.py $file1 $prefix\_vs_smallRNA_good_pairs.sorted.bed > $prefix\_vs_smallRNA_good_pairs_prototype.sorted.3p
end


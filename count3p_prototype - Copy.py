import sys

with open(sys.argv[1], 'r') as file1:
    genome_hash = {}
    for line in file1:
        chrom, length = line.strip().split('\t')
        genome_hash[chrom] = [0] * length
file1.close()
with open(sys.argv[2], 'r') as file2:
    for line in file2:
        chrom, start, end, location, quality, strand = line.strip().split('\t')
        genome_hash[chrom][end] += 1
file2.close()
for keys, value in genome_hash:
    print(keys, value, sep='\t')
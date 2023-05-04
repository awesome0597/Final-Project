import sys

genome_hash = {}

# actual genome file
with open(sys.argv[1], 'r') as file1:
    for line in file1:
        chrom, length = line.strip().split('\t')
        genome_hash[chrom] = [0] * int(length)
file1.close()
with open(sys.argv[2], 'r') as file2:
    for line in file2:
        chrom, start, end, location, quality, strand = line.strip().split('\t')
        genome_hash[chrom][int(start)] += 1
file2.close()
for keys, array in genome_hash.items():
    count = 0
    for index, value in enumerate(array):
        print(keys, count, value, sep='\t')
        count += 1

import statistics
import sys
import numpy as np
import pandas as pd
#from Bio import SeqIO


def createFile(gene_list, score_a, score_b, score_c, my_init, my_3p, my_cov):
    # , fasta, sno_data
    # Need to generalize for multiple files, currently works only for one file, add fasta file when accessible
    # sno_rna = sno_data['snoRNA'].tolist()
    # modification = sno_data['modified'].tolist()
    number_list = []
    for i in range(1, len(score_b) + 1):
        number_list.append(i)
    # fasta = fasta[:len(my_init)]
    genes_list = gene_list[:len(my_init)]
    df = pd.DataFrame({'Gene': genes_list, ' ': number_list, '5p': my_init, '3p': my_3p,
                       'cov': my_cov, 'Sa': score_a, 'Sb': score_b, 'Sc': score_c})
    # 'bp': fasta, 'modification': modification[:len(my_init)], 'snoRNA': sno_rna[:len(my_init)]}
    df.to_excel("output.xlsx", index=False)


def stats(my_cov, start, end):
    mean = statistics.fmean(my_cov[int(start):int(end)])
    std = statistics.stdev(my_cov[int(start):int(end)])
    return mean, std


def calculateScores(my_cov, my_length):
    # we can consider using a numpy array instead of a list
    score_a = [0] * my_length
    score_b = [0] * my_length
    score_c = [0] * my_length

    for i in range(win_size + 1, my_length - win_size):
        # A Score
        m_l, s_l = stats(my_cov, i - win_size / 2, i)
        m_r, s_r = stats(my_cov, i, i + win_size / 2)

        score_a[i] = max(0, 1 - (2 * my_cov[i] + 1) / (0.5 * abs(m_l - s_l) + my_cov[i] + 0.5 * abs(m_r - s_r) + 1))

        # B + C Score
        s1 = 0
        for j in range(i - win_size, i):
            s1 += (1 - 0.1 * (j - 1)) * my_cov[i - j]

        s1 = s1 / W
        s2 = 0
        for j in range(1, win_size):
            s2 += (1 - 0.1 * (j - 1)) * my_cov[i + j]
        s2 = s2 / W
        try:
            score_c[i] = max(0, 1 - (2 * my_cov[i]) / (s1 + s2))
        except ZeroDivisionError:
            score_c[i] = np.nan
        try:
            score_b[i] = abs((my_cov[i] - 0.5 * (s1 + s2)) / my_cov[i] + 1)
        except ZeroDivisionError:
            score_b[i] = np.nan
    return score_a, score_b, score_c


def covAndLen(init_library, three_p_library):
    pre_my_init = pd.read_table(init_library, header=None, usecols=[2])
    pre_my3p = pd.read_table(three_p_library, header=None, usecols=[2])
    # shifting reads by 1 bp
    my_new_init = pre_my_init[2].tolist()
    for i in range(0, len(my_new_init) - 1):
        my_new_init[i] = my_new_init[i + 1]
    my_new_3p = pre_my3p[2].tolist()
    cov = []
    for i in range(0, len(my_new_init)):
        cov.append(my_new_init[i] + my_new_3p[i])

    return my_new_init, my_new_3p, cov, len(my_new_init)


if __name__ == '__main__':
    # The program takes a genome file, init file, 3p file, fasta file path and known snoRNA info as arguments
    win_size = 6
    W = ((1 + (1 - 0.1 * win_size)) * win_size) / 2

    # process genomes to work by size
    gene_list_per_base_pair = []
    with open(sys.argv[1], 'r') as file1:
        for line in file1:
            chrom, rna_length = line.strip().split('    ')
            genes_to_add = [chrom] * (int(rna_length) + 1)
            gene_list_per_base_pair.extend(genes_to_add)
    file1.close()
    init_file, trep_file = sys.argv[2], sys.argv[3]
    # handle fasta file
    # fasta_file_path = sys.argv[4]
    myfasta = []
    # with open(fasta_file_path) as handle:
    #    for record in SeqIO.parse(handle, "fasta"):
    #        myfasta.append(str(record.seq))
    # handle.close()
    # fasta_as_list = []
    # fasta_as_list[:0] = myfasta[0]
    # sno_df = pd.read_table(sys.argv[5], delimiter="\t")
    # new_sno_df = sno_df[['modified', 'snoRNA']]

    # would consider using the code from count init\3p that found the different libraries and
    # the length of the RNA, and then save them as a tuple in a list

    # also need to see how we generalize the code for the fasta file, maybe consider putting it
    # as part of the array? done in theory

    myinit, my3p, mycov, mylength = covAndLen(init_file, trep_file)
    Sa, Sb, Sc = calculateScores(mycov, mylength)
    createFile(gene_list_per_base_pair, Sa, Sb, Sc, myinit, my3p, mycov)
    # , fasta_as_list
# , new_sno_df

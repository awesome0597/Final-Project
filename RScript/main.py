import statistics
import sys
import re
import os
import pandas as pd
from Bio import SeqIO


def createFile(library, Sa, Sb, Sc, myinit, my3p, mycov):
    ## Need to generalize for multiple files, currently works only for one file

    df = pd.DataFrame({'5p': myinit, '3p': my3p,
                       'cov': mycov, 'Sa': Sa, 'Sb': Sb, 'Sc': Sc})
    df.to_csv(f'{library}.csv', index=False)


def stats(mycov, start, end):
    mean = statistics.fmean(mycov[start:end])
    std = statistics.stdev(mycov[start:end])
    return mean, std


def calculateScores(mycov, mylength, win_size, W):
    # we can consider using a numpy array instead of a list
    Sa = [0] * mylength
    Sb = [0] * mylength
    Sc = [0] * mylength

    for i in range(win_size + 1, mylength - win_size):
        # A Score
        M_l, S_l = stats(mycov, i - win_size / 2, i)
        M_r, S_r = stats(mycov, i, i + win_size / 2)

        Sa[i] = max(0, 1 - (2 * mycov[i] + 1) / (0.5 * abs(M_l - S_l) + mycov[i] + 0.5 * abs(M_r - S_r) + 1))

        # B + C Score
        S1 = 0
        for j in range(i - win_size, i):
            S1 += (1 - 0.1 * (j - 1)) * mycov[i - j]

        S1 = S1 / W
        S2 = 0
        for j in range(1, win_size):
            S2 += (1 - 0.1 * (j - 1)) * mycov[i + j]
        S2 = S2 / W
        Sc[i] = max(0, 1 - (2 * mycov[i]) / (S1 + S2))
        Sb[i] = abs((mycov[i] - 0.5 * (S1 + S2)) / mycov[i] + 1)
    return Sa, Sb, Sc


def libraries(directory_path):
    file_pattern = re.compile(r'sorted\.init$')
    files = [f for f in os.listdir(directory_path) if file_pattern.match(f)]

    return [re.sub(r'\.sorted\.init$', '', f) for f in files]


def covAndLen(library, rna_length):
    pre_myinit = pd.read_table(library + ".sorted.init", header=None, usecols=[2])
    pre_my3p = pd.read_table(library + ".sorted.3p", header=None, usecols=[2])
    # shifting reads by 1 bp
    myinit = [0] + pre_myinit + ([0] * (rna_length - len(pre_myinit) + 1))
    my3p = pre_my3p[2:] + ([0] * (rna_length - len(pre_my3p) + 1))

    return myinit, my3p, myinit + my3p, len(myinit)


if __name__ == '__main__':
    ## The program takes a window size and a directory path as arguments
    win_size = int(sys.argv[1])
    W = ((1 + (1 - 0.1 * win_size)) * win_size) / 2

    directory_path = sys.argv[2]

    # would consider using the code from count init\3p that found the different libraries and
    # the length of the RNA, and then save them as a tuple in a list
    # also need to see how we generalize the code for the fasta file, maybe consider putting it
    # as part of the array?

    mylibs = libraries(directory_path)
    rna_length = int(sys.argv[2])
    file_path = "C:/Users/User/OneDrive - Bar Ilan University/Bioinfo Computer/Shula/HydraPsi/Files/TB_rRNA_chr2.fa"
    myfasta = []
    with open(file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            myfasta.append(str(record.seq))


    for library in mylibs:
        myinit, my3p, mycov, mylength = covAndLen(library, rna_length)
        Sa, Sb, Sc = calculateScores(mycov, mylength, win_size, W)
        createFile(library, Sa, Sb, Sc, myinit, my3p, mycov, mylength)

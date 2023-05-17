import statistics
import sys
import numpy as np
import pandas as pd
# from Bio import SeqIO
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QFormLayout, QGroupBox, QLineEdit, QDialogButtonBox, \
    QComboBox


def createFile(gene_list, score_a, score_b, score_c, my_init, my_3p, my_cov, number_list_appended):
    # , fasta, sno_data
    # Need to generalize for multiple files, currently works only for one file, add fasta file when accessible
    # sno_rna = sno_data['snoRNA'].tolist()
    # modification = sno_data['modified'].tolist()
    # fasta = fasta[:len(my_init)]
    genes_list = gene_list[:len(my_init)]
    df = pd.DataFrame({'Gene': genes_list, ' ': number_list_appended, '5p': my_init, '3p': my_3p,
                       'cov': my_cov, 'Sa': score_a, 'Sb': score_b, 'Sc': score_c})
    # 'bp': fasta, 'modification': modification[:len(my_init)], 'snoRNA': sno_rna[:len(my_init)]}
    df.to_excel("Ribosomal_test.xlsx", index=False)


def stats(my_cov, start, end):
    mean = statistics.fmean(my_cov[int(start):int(end)])
    std = statistics.stdev(my_cov[int(start):int(end)])
    return mean, std


def calculateScores(my_number_list, my_cov, my_length, win_size, W):
    # we can consider using a numpy array instead of a list
    score_a = [0] * my_length
    score_b = [0] * my_length
    score_c = [0] * my_length

    for i in range(0, my_length - win_size):
        if my_number_list[i] < 5 or my_number_list[i + 5] < 5:
            score_a[i] = "NA"
            score_b[i] = "NA"
            score_c[i] = "NA"
        else:
            # A Score
            m_l, s_l = stats(my_cov, i - win_size / 2, i)
            m_r, s_r = stats(my_cov, i + 1, (i + 1) + win_size / 2)
            if i == 3436:
                print("hi")

            score_a[i] = max(0, 1 - (2 * my_cov[i] + 1) / (0.5 * abs(m_l - s_l) + my_cov[i] + 0.5 * abs(m_r - s_r) + 1))

            # B + C Score
            s1 = 0
            for j in range(1, win_size + 1):
                s1 += (1 - 0.1 * (j - 1)) * my_cov[i - j]

            s1 = s1 / W
            s2 = 0
            for j in range(1, win_size + 1):
                s2 += (1 - 0.1 * (j - 1)) * my_cov[i + j]
            s2 = s2 / W
            try:
                score_c[i] = max(0, 1 - (2 * my_cov[i]) / (s1 + s2))
            except ZeroDivisionError:
                score_c[i] = "NA"
            try:
                score_b[i] = abs((my_cov[i] - 0.5 * (s1 + s2)) / (my_cov[i] + 1))
            except ZeroDivisionError:
                score_b[i] = "NA"
    return score_a, score_b, score_c


def covAndLen(init_library, three_p_library):
    pre_my_init = pd.read_table(init_library, header=None, usecols=[2])
    pre_my3p = pd.read_table(three_p_library, header=None, usecols=[2])
    init_to_move = pre_my_init[2].tolist()
    my_new_init = pre_my_init[2].tolist()
    my_new_3p = pre_my3p[2].tolist()
    my_new_3p.append(0)
    init_to_move.append(0)
    my_new_init.append(0)
    for i in range(0, len(my_new_init) - 2):
        my_new_init[i + 2] = init_to_move[i]

    cov = []
    for i in range(0, len(my_new_init)):
        cov.append(my_new_init[i] + my_new_3p[i])

    return my_new_init, my_new_3p, cov, len(my_new_init)


def runScript(sequencing_type, window_size, genome_file_path, init_file_path, three_p_file_path, fasta_file_path,
              output_file_name):
    # The program takes a genome file, init file, 3p file, fasta file path and known snoRNA info as arguments
    W = (1 + (1 - 0.1 * window_size)) * window_size / 2

    # process genomes to work by size
    gene_list_per_base_pair = []
    number_list = []
    with open(genome_file_path, 'r') as file1:
        for line in file1:
            chrom, rna_length = line.strip().split('\t')
            genes_to_add = [chrom] * (int(rna_length))
            gene_list_per_base_pair.extend(genes_to_add)
            for p in range(0, int(rna_length)):
                number_list.append(p)
    file1.close()
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

    myinit, my3p, mycov, mylength = covAndLen(init_file_path, three_p_file_path)
    Sa, Sb, Sc = calculateScores(number_list, mycov, mylength, window_size,W)
    createFile(gene_list_per_base_pair, Sa, Sb, Sc, myinit, my3p, mycov, number_list)
    # , fasta_as_list
    # , new_sno_df


class MyWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        # Create the main layout
        layout = QVBoxLayout()

        # Create the group box
        group_box = QGroupBox("Group Box")
        form_layout = QFormLayout()

        # Create the QLineEdit and QComboBox widgets
        combo_box = QComboBox()
        combo_box.addItem("PRS")
        combo_box.addItem("Total RNA")
        labels_and_inputs = [("Sequencing Type:", combo_box), ("Window Size:", QLineEdit()),
                             ("Genome File Path:", QLineEdit()),
                             ("init File Path:", QLineEdit()), ("3p File Path:", QLineEdit()),
                             ("Fasta File Path:", QLineEdit()), ("Output File Name:", QLineEdit())]

        for i in range(len(labels_and_inputs)):
            label, input_var = labels_and_inputs[i]

            form_layout.addRow(label, input_var)

        group_box.setLayout(form_layout)

        # Add the group box to the main layout
        layout.addWidget(group_box)

        # Create the QDialogButtonBox
        button_box = QDialogButtonBox()
        button_box.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        # Connect the Ok button's clicked signal to the on_ok_button_clicked function
        button_box.accepted.connect(self.on_ok_button_clicked)

        # Add the button box to the main layout
        layout.addWidget(button_box)

        self.setLayout(layout)

    def on_ok_button_clicked(self):
        # Retrieve the inputs from the QLineEdit widgets and the selected option from the QComboBox
        inputs = [self.findChild(QComboBox).currentText()]

        line_edits = self.findChildren(QLineEdit)
        for line_edit in line_edits:
            inputs.append(line_edit.text())  # Get the text from each QLineEdit

        # Pass the inputs to another function for further processing
        process_inputs(inputs)

        # Close the dialog or perform other actions as needed
        self.close()


def process_inputs(inputs):
    # Perform your logic on the inputs
    sequencing_type = inputs[0]
    window_size = int(inputs[1])
    genome_file_path = inputs[2]
    init_file_path = inputs[3]
    three_p_file_path = inputs[4]
    fasta_file_path = inputs[5]
    output_file_name = inputs[6]

    runScript(sequencing_type, window_size, genome_file_path, init_file_path, three_p_file_path, fasta_file_path,
              output_file_name)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    widget = MyWidget()
    widget.show()
    sys.exit(app.exec_())

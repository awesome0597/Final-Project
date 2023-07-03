import statistics
import sys
import pandas as pd
from Bio import SeqIO
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QFormLayout, QGroupBox, QLineEdit, QDialogButtonBox, \
    QComboBox


def createFile(gene_list, score_a, score_b, score_c, my_init, my_3p,
               my_cov, number_list_appended, output_file_name, fasta):
    #Create a dataframe with the sequencing for all the different genes
    df = pd.DataFrame({'Gene': gene_list, 'Position': number_list_appended, 'bp': fasta, '5p': my_init, '3p': my_3p,
                       'cov': my_cov, 'Sa': score_a, 'Sb': score_b, 'Sc': score_c})
    # Save file to excel
    df.to_excel(output_file_name, index=False)


def stats(my_cov, start, end):
    #calculate the mean and standard deviation of the coverage
    mean = statistics.fmean(my_cov[int(start):int(end)])
    std = statistics.stdev(my_cov[int(start):int(end)])
    return mean, std


def calculateScores(my_number_list, my_cov, my_length, win_size, w):
    #calculate the scores for each position based on the article
    score_a = [0] * my_length
    score_b = [0] * my_length
    score_c = [0] * my_length

    for i in range(0, my_length - win_size):
        # Set NA for the first and last 5 positions of the gene since they cannot be properly calculated
        if my_number_list[i] < 6 or my_number_list[i + 5] < 7:
            score_a[i] = "NA"
            score_b[i] = "NA"
            score_c[i] = "NA"
        else:
            # A Score
            m_l, s_l = stats(my_cov, i - win_size / 2, i)
            m_r, s_r = stats(my_cov, i + 1, (i + 1) + win_size / 2)

            score_a[i] = max(0, 1 - (2 * my_cov[i] + 1) / (0.5 * abs(m_l - s_l) + my_cov[i] + 0.5 * abs(m_r - s_r) + 1))

            # B + C Score
            s1 = 0
            for j in range(1, win_size + 1):
                s1 += (1 - 0.1 * (j - 1)) * my_cov[i - j]

            s1 = s1 / w
            s2 = 0
            for j in range(1, win_size + 1):
                s2 += (1 - 0.1 * (j - 1)) * my_cov[i + j]
            s2 = s2 / w
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
    my_new_init = pre_my_init[2].tolist()
    my_holder_init = pre_my_init[2].tolist()
    my_holder_3p = pre_my3p[2].tolist()
    my_new_3p = pre_my3p[2].tolist()
    for i in range(0, len(my_new_init) - 1):
        my_new_init[i + 1] = my_holder_init[i]
    my_new_init[0] = 0
    for i in range(len(my_new_3p) - 1, 1, -1):
        my_new_3p[i-1] = my_holder_3p[i]
    my_new_3p[len(my_new_3p) - 1] = 0

    cov = []
    for i in range(0, len(my_new_init)):
        cov.append(my_new_init[i] + my_new_3p[i])

    return my_new_init, my_new_3p, cov, len(my_new_init)

#TODO: add the difference between the two sequencing types (PRS and Total RNA)
def runScript(sequencing_type, window_size, genome_file_path, init_file_path, three_p_file_path, fasta_file_path,
              output_file_name):
    # The program takes a genome file, init file, 3p file, fasta file path and known snoRNA info as arguments

    #  set the value of W as per the article
    w = (1 + (1 - 0.1 * window_size)) * window_size / 2

    # process genomes to work by size
    #list that holds all the different genes and their lengths
    gene_list_per_base_pair = []
    #list that holds the index for different positions in each gene
    number_list = []
    with open(genome_file_path, 'r') as file1:
        for line in file1:
            chrom, rna_length = line.strip().split()
            genes_to_add = [chrom] * (int(rna_length) + 1)
            gene_list_per_base_pair.extend(genes_to_add)
            for p in range(0, int(rna_length) + 1):
                number_list.append(p)
    file1.close()
    # handle fasta file
    myfasta = []
    fasta_as_list = []
    with open(fasta_file_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            myfasta.append(str(record.seq))
    handle.close()
    for fasta_string in myfasta:
        fasta_string = "<" + fasta_string
        fasta_as_list.extend(list(fasta_string))

    myinit, my3p, mycov, mylength = covAndLen(init_file_path, three_p_file_path)
    sa, sb, sc = calculateScores(number_list, mycov, mylength, window_size, w)
    createFile(gene_list_per_base_pair, sa, sb, sc, myinit, my3p, mycov, number_list, output_file_name, fasta_as_list)


class MyWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        # Create the main layout
        layout = QVBoxLayout()
        # Name for window
        self.setWindowTitle("Coverage Score Calculations")
        self.resize(400, 300)

        # Create the group box
        group_box = QGroupBox("Please Enter the Following Parameters:")
        form_layout = QFormLayout()

        # Create the QLineEdit and QComboBox widgets
        combo_box = QComboBox()
        combo_box.addItem("PRS")
        combo_box.addItem("Total RNA")
        # set default value for window size
        window_size = QLineEdit()
        window_size.setText('6')
        # set default value for output file name
        output_file_name = QLineEdit()
        output_file_name.setText('temp.xlsx')
        # set default genome file
        genome_file_path = QLineEdit()
        genome_file_path.setText('PRS_genome.txt')
        # set fasta file
        fasta_file_path = QLineEdit()
        fasta_file_path.setText('TB_small_RNAs_DB_w_praveen.fa')
        labels_and_inputs = [("Sequencing Type:", combo_box), ("Window Size:", window_size),
                             ("Genome File Path:", genome_file_path),
                             ("init File Path:", QLineEdit()), ("3p File Path:", QLineEdit()),
                             ("Fasta File Path:", fasta_file_path), ("Output File Name:", output_file_name)]

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
        button_box.rejected.connect(self.on_cancel_button_clicked)

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

    def on_cancel_button_clicked(self):
        self.close()


def process_inputs(inputs):
    # receive the inputs from the user and run the script
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

# 942PRS16012023_S7_vs_smallRNA_good_pairs_prototype.sorted.init
# TB_small_RNAs_DB_w_praveen.fa
# actual_genome.txt

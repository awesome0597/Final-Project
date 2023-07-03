# PRS Pipeline and Scoring System
## Pipeline:
This program is a NGS pipeline for processing RNAseq files.
### Sequencing Types:
* Total RNA
* PRS<br>
### Program Instructions:
This script is designed to process FASTQ files and perform various operations on them. Follow the steps below to run the script:

1. Make sure you have the required dependencies installed:
   - `smalt` for read mapping
   - `samtools` for manipulating SAM/BAM files
   - `bedtools` for converting and processing BED files
   - `python3` with the necessary Python packages for executing the Python scripts

2. Set the necessary paths and variables:
   - Edit the script file and modify the following variables at the beginning:
     - `path1`: Path to the directory containing the FASTQ files
     - `path2`: Path to the directory containing additional scripts (in this case the init and 3p scripts)
     - `type`: Specify the type of data being processed (e.g., "smallRNA")
     - `snrna`: Path to the small RNA file for read mapping
     - `file1`: Path to the genome file
3. Run the script:
   - Open a terminal and navigate to the directory containing the script file.
   - Execute the script using the following command:
     ```
     ./RunPiplinePythonVersion_mapping.sh
     ```
   - The script will process the FASTQ files, perform read mapping, filter and convert the resulting SAM files, and generate intermediate and final output files.
   - Then execute the second part of the script using the following command:
     ```
     ./RunPiplinePythonVersion_init_3p.sh
     ```
   - This step will process the BAM files generated in the previous step and give us the necesary 3p and init files.

## Additional Information

- Make sure the paths to the input files and directories are correct.
- Ensure that you have the necessary permissions to read and write files in the specified directories.
- This program is intended to run in a linux enviroment.<br>

#### output files<br>
* bam file for every fastq file the program runs on
* bed file extracted from the recived bam from the mapping part of the pipline
* genomecov file computed from the recived bed file
* init file presenting read ends per nucliotide of every mapped file.
* 3p file presenting read beginings per nucliotide of every mapped file.
____________________________________________________________________________________________________________________________________________________________________
## Scoring system:
The Coverage Score Calculations script is a Python program that calculates coverage scores for different genes based on provided sequencing data. It takes into account various parameters such as the type of sequencing, window size, and input files, and generates an output file with the calculated scores.

### Required Packages:
Python 3.x
pandas
BioPython
PyQt5
### Installation
Clone the repository or download the script file.
Install the required dependencies using pip:
```
pip install pandas biopython pyqt5
```

This program returns scoring based on the methods mentioned in the file: anie_201408362_sm_miscellaneous_information(4).pdf attached.<br>
### Sequencing Types:
* Total RNA
* PRS<br>
### Program Instructions:
Run the script using the following command:
```
python3 Scores_script.py
```
The script will launch a graphical user interface (GUI) window.

Enter the required parameters in the provided fields:
- Sequencing Type: Select the type of sequencing from the dropdown menu (PRS or Total RNA).
- Window Size: Enter the window size for score calculations (default is 6).
- Genome File Path: Enter the path to the genome file containing chromosome and RNA length information.
- init File Path: Enter the path to the init file containing initial library data.
- 3p File Path: Enter the path to the 3p file containing 3' library data.
- Fasta File Path: Enter the path to the FASTA file containing genomic sequences.
- Output File Name: Enter the desired name for the output Excel file (default is temp.xlsx).

The script will process the inputs and generate the output file with coverage scores.

Once the calculations are complete, the GUI window will close automatically.

## Additional Information
*It is recommended to run the program in a windows enviroment

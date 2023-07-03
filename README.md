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

## Scoring system:
This program returns scoring based on the methods mentioned in the article: [insert article].<br>
### Sequencing Types:
* Total RNA
* PRS<br>
### Required Packages:
1. Total RNA
2. PRS<br>
### Program Instructions:
#### file orginization<br>
#### output files<br>
#### run instructions<br>

*It is recommended to run the program in a windows enviroment

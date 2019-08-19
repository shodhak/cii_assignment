## CII Programming Assignment
### Submitted by : Shreyas Joshi

When working with any new dataset, the most important thing is to understand the data structure. Here there are two files from 454 sequencing output. One is a fasta file (test.fna), and a qual file (test.qual). Whenever provided with a file containing data, the first step is to parse the file. Based on the questions that were asked, my initial approach was to make a dataframe containing columns with sequence identifiers, lengths of sequences, and the actual sequences themselves. Since most of my experience is in R, I would have done this by scraping the file first, and then separate the individual records using ">" as delimiter. This would be followed by appending information from individual reads into the dataframe. Similar approach would be used for qual file followed by combining sequence and qual information based on identifiers. I have implemented a version of this approach via pandas for one part of this assignment.

For this assignment, I have used Biopython to answer the questions based on the provided dataset. Biopython contains great resources for using python for bioinformatics and especially for input/output, handling and analysis of NGS data. There are 3 python scripts to generate the required output and complete majority of tasks, and one bash script for standalone blast. I have tried to explain programming logic and code output for all of them.

### *Script 1: cii_assignment.py*
This is the main script that provides output from read and quality realted parameters from original dataset, and also generates the main output file. This script uses **SeqIO module from Bio package** to read, parse, and write files. In the first step, the program asks the user to provide path to the directory containing the required files. This will be set as the working directory. It also asks the user for the names of fasta and qual files, which are then stored in their own variables.

```python
#Ask the user for path to the working directory, and names of fasta and qual files
wrk_dir = input("Enter path to working directory that contains fasta and qual files: ")
fasta_file = input("Enter name of the fasta file containing the reads: ")
qual_file = input("Enter name of qual file containing the quality scores for reads: ")

#Assign user provided path as working directory
os.chdir(wrk_dir)
#Print current working directory and the names of fasta and qual files
print("The current working directory is %s" % os.getcwd())
print("The fasta file name is %s" % fasta_file)
print("The qual file name is %s" % qual_file)

#Store original dataset filename in a variable
original_dataset = fasta_file
```

### 1) Total number of reads in the original dataset

Create a list object of SeqIO parser output of the fasta file, and use len function to compute the number of reads. There are **43090** reads in the original dataset file.

```python
from Bio import SeqIO
#Find total number of reads in the file
reads = list(SeqIO.parse(original_dataset,"fasta"))
print("Total number of reads in the original dataset : %s" % len(reads))
```
For the remaining questions, the identifier, length and sequences were parsed and saved into lists. This is same programming logic as mentioned earlier, just that here there are lists instead of a data frame.

```python
#Read quality scores for sequences from the qual file
#Initalize variables for sequence id and base quality data
seqid_qual = []
qual_data = []
average_qual = []
for ind_qual in SeqIO.parse(qual_file, "qual"):
    seqid_qual.append(ind_qual.id)
    qual_data.append(ind_qual.letter_annotations["phred_quality"])
    average_qual.append(mean(ind_qual.letter_annotations["phred_quality"]))
```

### 2) Total number of reads greater than 100 bp in the original dataset

Filter reads based on length greater than 100. There are **38638** reads with sequence lengths greater than 100.
```python
#Compute total number of reads where bp length is greater than 100
seq_len_100 = sum(length > 100 for length in lengths)
print("Total number of reads greater than 100 bp in the original dataset: %s" % seq_len_100)
```

### 3) Total number of reads with average quality scores greater than 20 in the original dataset

SeqIO reads qual files as well and this time, lists of identifiers and quality scores were made, and average quaity scores were calculated using mean function. These were then filtered for values above 20. There are **43056** reads with average quality scores greater than 20.
```python
#Read quality scores for sequences from the qual file
#Initalize variables for sequence id and base quality data
seqid_qual = []
qual_data = []
average_qual = []
for ind_qual in SeqIO.parse(qual_file, "qual"):
    seqid_qual.append(ind_qual.id)
    qual_data.append(ind_qual.letter_annotations["phred_quality"])
    average_qual.append(mean(ind_qual.letter_annotations["phred_quality"]))

#Compute total number of reads where average score quality is more than 20
len_quality_20 = sum(average_quality > 20 for average_quality in average_qual )
print("Total number of reads with average quality scores greater than 20 in the original dataset: %s" % len_quality_20)
```

### 4) Total number of reads with primer sequences

Check which reads have primer sequence at the start. There are **34059** reads with primers at the begining of sequences.
```python
#Identify the reads with primer sequences
#Assign variable for primer sequence
primer = "CGCCGTTTCCCAGTAGGTCTC"
reads_with_primer = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if sequence.seq.startswith(primer))
#reads_with_primer is a generator expression. Convert it first to list to get the number of sequences containing primer.
print("Total number of reads with primer sequences : %s" % len(list(reads_with_primer)))
```

### 5) Total number of reads with adaptor sequences

When the adaptor sequence provided ("ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG") was searched for matches in the reads, there were no hits. The blast output gave matches for 20 nucleotides between positions 14 and 33 on the adaptor sequence. This string with 20 nucleotides ("CAAGGCACACAGGGGATAGG") was used for adaptor trimming.There are **3440** reads containing the adaptor sequence.

```python
#Identiy reads with adaptor sequence
#Assign variable for adaptor sequence
adaptor = "ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG"
reads_with_adaptor = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if sequence.seq.find(adaptor[13:33]) != -1)
#reads_with_adaptor is a generator expression. Convert it first to list to get the number of sequences containing adaptor.
print("Total number of reads with adaptor sequence : %s" % len(list(reads_with_adaptor)))
```

### 6) Total number of reads with both primer and adaptor sequences

Conditions for both primer and adaptor searching from previous two steps were combined to identify reads containing both these sequences. There are **3339** reads that contain both primer and adaptor sequences.

```python
#Identify reads with both primer and adaptor sequence
reads_with_both = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if (sequence.seq.startswith(primer) and (sequence.seq.find(adaptor[13:33]) != -1)))
#reads_with_both is a generator expression. Convert it first to list to get the number of sequences containing both primer and adaptor.
print("Total number of reads with both primer and adaptor sequences : %s" % len(list(reads_with_both)))
```
### Generating main output
The remaining part of the script generates Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed. This is achieved in a stepwise manner:
1. *Filter reads for sequence length greater than 100*
2. *Filter reads for quality scores greater than 20*
3. *Trim primer*
4. *Trim adaptor*

The functions for trimming primer and adaptor sequences used in the script are from the [Biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html). The final output file (main_output.fasta) contains **38053** reads.

### *Script 2: comparison_plots.py*
Plots for length and quality data in original and output datasets provide a good visual representation of how those parameters have changed, and is also a good sanity check to see if the reads have been filtered correctly. **Matplotlib** package was used to make these plots. 

### Read lengths
The frequency of reads was plotted by incrementally binning them in read length sizes of 20 nucleotides. Plots for both original and output datasets were made using this approach. The range of lengths in main output fasta dataset begins at 100 which means that our code worked correctly for filtering reads based on lengths.

#### Read lengths in original dataset
![Lengths Original](https://drive.google.com/uc?export=view&id=1aocL51WY7oSVyyENtF8h4YievPce9bpC)

#### Read lengths in main output
![Lengths Output](https://drive.google.com/uc?export=view&id=1UIp7i5yFBIcI6YL2FxdY8mF9N2sCv3DH)

#### Average quality scored in original dataset
![Quality Original](https://drive.google.com/uc?export=view&id=1lYtzw9XMqJOl5d_1oV65Szlz6iI4olap)

#### Average quality scored in main output
![Quality Output](https://drive.google.com/uc?export=view&id=18BwO4LCDLjRUNqqmd_dimUwkLgLerer-)

### Start and end positions of primer and adaptor
The goal here was to make a dataframe containing sequence identifiers for all reads and check whether primer and/or adaptor were present or not. In case they were present, their start and end positions were stored. **Pandas** package was used to make this dataframe. As mentioned earlier, majority of questions pertaining to the lengths and presence and absence of primer and adaptor sequences can be answered through such a data structure. 

### *Script 3: df_make.py*
The output of this script provides start and end positions for primer and adaptor sequences if they are present in the read sequence. If they are absent, then "NA" is assigned. The sequence identifier is the index column. The initial dataframe also contains sequence lengths and actual sequences and the results from *cii_assignment.py* were double checked. The values were same as what we found earlier:
* Number of reads = 43090
* Number of reads with lengths greater than 100 = 38638
* Number of reads containing the primer = 34059
* Number of reads containing adapter = 3440
* Number of reads containing both = 3339

The length and sequence column are dropped from the final dataframe which is then saved as a csv file *"start_end.csv"*.

### BLAST
Standalone blast was performed for primer and adapter sequences against the database built from original fasta file *"test.fna"*. The output was stored in tabular format using *-outfmt 6* option in *blastn*.

### *Script 4: blast_run.sh*
Blast database was built for *"test.fna"*, and primer and adaptor sequences were stored in a separate fasta file *"prim_adapt.fasta"*. Blast output in tabular format was stored in *"test_blast.tab"*.

```shell
#!/bin/bash

blast_path="/usr/local/ncbi/blast/bin"
#Make blast database with primer and adaptor sequences
$blast_path/makeblastdb -in test.fna -title testdb -dbtype nucl -parse_seqids -out testdb

#Run nucleotide blast and produce output in tabular format
#m8 is legacy option in blast for tabular output. In BLAST+, -outfmt 6 does the same job
$blast_path/blastn -query prim_adapt.fasta -db testdb -task "blastn-short" -outfmt 6 -out test_blast.tab
```


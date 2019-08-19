import os
from Bio import SeqIO
from statistics import mean
import numpy as np
import matplotlib.pyplot as plt

#Set current working directory
os.chdir("/Users/sj577/Documents/cii_assignment/Files_for_test")
#Check current working directory
print(os.getcwd())

#LENGTH PLOT FOR ORIGINAL DATASET
#Store original dataset and main output filenames in a variables
original_dataset = "test.fna"
main_output = "main_output.fasta"

#Initiate variables that will store the sequence id, length, and sequence in the fasta file
seqid_fna = []
lengths = []

#Extract information for sequence id, length, and sequence for individual sequences in the fasta file
for ind_seq in SeqIO.parse(original_dataset,"fasta"):
     seqid_fna.append(ind_seq.id)
     lengths.append(len(ind_seq))

#Plot distribution of basepair lengths of sequences in original dataset
plt.figure()
#Assign ranges for x and y axis
plt.axis([0, 600, 0, 5000])
#Define size of bins
# The maximum value is 6000
bin_size = np.linspace(0,600,21)
plt.hist(lengths, bins=bin_size)
#Define plot parameters
plt.title("Original dataset with %i reads\nRange of lengths: %i to %i" % (len(lengths),min(lengths),max(lengths)))
plt.xlabel("Sequence Length")
plt.ylabel('Frequency')
plt.savefig("len_original_dataset.png")

#LENGTH PLOT FOR MAIN OUTPUT

#Initiate variables that will store the sequence id, length, and sequence in the fasta file
seqid_fna = []
lengths = []

#Extract information for sequence id, length, and sequence for individual sequences in the fasta file
for ind_seq in SeqIO.parse(main_output,"fasta"):
     seqid_fna.append(ind_seq.id)
     lengths.append(len(ind_seq))

#Plot distribution of basepair lengths of sequences in original dataset
plt.figure()
#Assign ranges for x and y axis
plt.axis([0, 600, 0, 5000])
#Define size of bins
# The maximum value is 6000
bin_size = np.linspace(0,600,21)
plt.hist(lengths, bins=bin_size)
#Define plot parameters
plt.title("Main output with %i reads\nRange of lengths: %i to %i" % (len(lengths),min(lengths),max(lengths)))
plt.xlabel("Sequence Length")
plt.ylabel('Frequency')
plt.savefig("len_main_output.png")

#QUALITY PLOT FOR ORIGINAL DATASET
#Read quality scores for sequences from the qual file
#Initalize variables for average base quality data
average_qual = []
for ind_qual in SeqIO.parse("test.qual", "qual"):
    average_qual.append(mean(ind_qual.letter_annotations["phred_quality"]))

#Plot average quality scores in the original dataset
plt.figure()
plt.plot(average_qual)
plt.ylim(0,40)
#Define plot parameters
plt.title("Original dataset quality plot for %i reads\nRange of scores: %i to %i" \
          % (len(average_qual),min(average_qual),max(average_qual)))
plt.xlabel("Reads")
plt.ylabel('Quality scores')
plt.savefig("qual_original_dataset.png")

#QUALITY PLOT FOR MAIN OUTPUT
#Read quality scores for sequences from the qual file
#Initalize variable for average base quality data
average_qual = []
for ind_qual in SeqIO.parse("filtered.qual", "qual"):
    average_qual.append(mean(ind_qual.letter_annotations["phred_quality"]))

#Plot average quality scores in main output
plt.figure()
plt.plot(average_qual)
plt.ylim(0,40)
#Define plot parameters
plt.title("Main Output quality plot for %i reads\nRange of scores: %i to %i" \
          % (len(average_qual),min(average_qual),max(average_qual)))
plt.xlabel("Reads")
plt.ylabel('Quality scores')
plt.savefig("qual_main_output.png")
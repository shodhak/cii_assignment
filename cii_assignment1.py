"""Set current working directory"""
import os
from Bio import SeqIO
from statistics import mean
import numpy as np
import matplotlib.pyplot as plt

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

#Find total number of reads in the file
reads = list(SeqIO.parse(original_dataset,"fasta"))
print("Total number of reads in the original dataset : %s" % len(reads))

#Initiate variables that will store the sequence id, length, and sequence in the fasta file
seqid_fna = []
lengths = []
sequences = []

#Extract information for sequence id, length, and sequence for individual sequences in the fasta file
for ind_seq in SeqIO.parse(original_dataset,"fasta"):
     seqid_fna.append(ind_seq.id)
     lengths.append(len(ind_seq))
     sequences.append(repr(ind_seq.seq))

#Compute total number of reads where bp length is greater than 100
seq_len_100 = sum(length > 100 for length in lengths)
print("Total number of reads greater than 100 bp in the original dataset: %s" % seq_len_100)

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

#Identify the reads with primer sequences
#Assign variable for primer sequence
primer = "CGCCGTTTCCCAGTAGGTCTC"
reads_with_primer = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if sequence.seq.startswith(primer))
#reads_with_primer is a generator expression. Convert it first to list to get the number of sequences containing primer.
print("Total number of reads with primer sequence : %s" % len(list(reads_with_primer)))

#Identiy reads with adaptor sequence
#Assign variable for adaptor sequence
adaptor = "ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG"
reads_with_adaptor = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if sequence.seq.find(adaptor[13:33]) != -1)
#reads_with_adaptor is a generator expression. Convert it first to list to get the number of sequences containing adaptor.
print("Total number of reads with adaptor sequence : %s" % len(list(reads_with_adaptor)))

#Identify reads with both primer and adaptor sequence
reads_with_both = (sequence for sequence in SeqIO.parse(original_dataset,"fasta") \
    if (sequence.seq.startswith(primer) and (sequence.seq.find(adaptor[13:33]) != -1)))
#reads_with_both is a generator expression. Convert it first to list to get the number of sequences containing both primer and adaptor.
print("Total number of reads with both primer and adaptor sequences : %s" % len(list(reads_with_both)))

#GENERATE MAIN OUTPUT FILE
#Filter reads with bp lengths greater than 100
#Identify indexes where length (number of bp) greater than 100
ind_len_100 = [len_index for len_index, value in enumerate(lengths) if value > 100]

#Extract identifiers for reads for which number of bp greater than 100
seqid_gt_100 = [seqid_fna[i] for i in ind_len_100]

#Generate a fasta file with reads greater than 100bp
len100 = (read for read in SeqIO.parse(original_dataset, "fasta") if read.id in seqid_gt_100)
count = SeqIO.write(len100, "main_len100.fasta", "fasta")
print("Number of reads in length filtered file where bp > 100: %s" % count)

#Filter reads with average quality scores higher than 20 in previous output which is the file generated from filtering by length
#Identify indexes where average quality scores are greater than 20
ind_qual_20 = [qual_index for qual_index, value in enumerate(average_qual) if value > 20]

#Extract identifiers for reads for which average quality scores greater than 20
seqqual_gt_20 = [seqid_fna[i] for i in ind_qual_20]

#From the previously generated main output, shortlist reads with average quality scores greater than 20bp
qual20 = (read for read in SeqIO.parse("main_len100.fasta", "fasta") if read.id in seqqual_gt_20)
count = SeqIO.write(qual20, "main_qual20.fasta", "fasta")
print("Number of reads in quality filtered file where average quality scores > 20: %s" % count)

#TRIM PRIMERS
#Function for primer trimming (from Biopython tutorial)
def trim_primer(record, primer):
    if record.seq.startswith(primer):
        return record[len(primer):]
    else:
        return record

primer_trimmed = (trim_primer(read, primer) for read in SeqIO.parse("main_qual20.fasta", "fasta"))
count = SeqIO.write(primer_trimmed, "main_primer_trimmed.fasta", "fasta")
print("Saved %i primer trimmed reads" % count)

#TRIM ADAPTORS
#function for adaptor trimming (from Biopython tutorial)
def trim_adaptors(records, adaptor, min_len):
    len_adaptor = len(adaptor) #cache this for later
    for record in records:
    #Check for sequence length
        len_record = len(record) #cache this for later
        if len(record) < min_len:
            #Too short to keep
            continue
        index = record.seq.find(adaptor)
        if index == -1:
            #adaptor not found, so won’t trim
            yield record
        elif len_record - index - len_adaptor >= min_len:
            #after trimming this will still be long enough
            yield record[index+len_adaptor:]

filtered_reads = SeqIO.parse("main_primer_trimmed.fasta", "fasta")
adaptor_trimmed = trim_adaptors(filtered_reads, adaptor, 100)
count = SeqIO.write(adaptor_trimmed, "main_output.fasta", "fasta")
print("Saved %i adaptor trimmed reads in main output" % count)

#Extract sequence identifiers from main output so as to select those reads from qual file.
#This file would be required for generating quality score plots in comparison_plots.py
mo_seqid = []
for fil_seq in SeqIO.parse("main_output.fasta","fasta"):
     mo_seqid.append(fil_seq.id)

mo_qual20 = (read for read in SeqIO.parse("test.qual", "qual") if read.id in mo_seqid)
count = SeqIO.write(mo_qual20, "filtered.qual", "qual")
print("Saved %i filtered quality records" % count)



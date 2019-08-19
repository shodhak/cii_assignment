import os
import pandas as pd
import numpy as np
from Bio import SeqIO
os.chdir("/Users/sj577/Documents/cii_assignment/Files_for_test")
"""Check current working directory"""
print(os.getcwd())

#The goal is to make a tab-delimited file containing sequence identifiers and starting and end positions for primer and adaptor

#Store original dataset filename in a variable
original_dataset = "test.fna"
primer = "CGCCGTTTCCCAGTAGGTCTC"
adaptor = "ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG"
print("Length of primer : %i" % len(primer))
print("Length of adaptor : %i" % len(adaptor[13:33]))

#Initiate variables that will store the sequence id, length, and sequence in the fasta file
seqid_fna = []
lengths = []
sequences = []

#Extract information for sequence id, length, and sequence for individual sequences in the fasta file
for ind_seq in SeqIO.parse(original_dataset,"fasta"):
     seqid_fna.append(ind_seq.id)
     lengths.append(len(ind_seq))
     sequences.append(str(ind_seq.seq))

#index = sequences.str.find(adaptor)
#Make a zipped list of sequence info
original_zl = list(zip(seqid_fna, lengths, sequences))

#Make data frame from zipped list
original_df = pd.DataFrame(original_zl, columns=['Identifier','Length','Sequence'])
#Set sequence identifier as index column
original_df.set_index('Identifier')
original_df['primer_start'] = original_df['Sequence'].str.find(primer)
original_df['primer_start'] = pd.to_numeric(original_df['primer_start'])
original_df['primer_start'] = original_df['primer_start'] + 1
original_df['primer_end'] = np.where(original_df['primer_start'] == 0, 0, original_df['primer_start'] + (len(primer)-1))
original_df['adaptor_start'] = original_df['Sequence'].str.find(adaptor[13:33])
original_df['adaptor_start'] = pd.to_numeric(original_df['adaptor_start'])
original_df['adaptor_start'] = original_df['adaptor_start'] + 1
original_df['adaptor_end'] = 0
original_df['adaptor_end'] = np.where(original_df['adaptor_start'] == 0, 0, original_df['adaptor_start'] + (len(adaptor[13:33])-1))

#Check the number of reads
print("Total number of reads : %i" % len(original_df.index))
#Check number of reads with lengths greater than 100
len_gt_100 = original_df.apply(lambda x: True if x['Length'] > 100 else False, axis=1)
nlen_gt_100 = len(len_gt_100[len_gt_100 == True].index)
print("Number of reads with lengths greater than 100 : %i" % nlen_gt_100)
#The length and sequence columns won't be part of final output, so they can be dropped
original_df = original_df.drop(columns=['Length','Sequence'])
#Double checking results from SeqIO
primer_present = original_df.apply(lambda x: True if x['primer_start'] == 1 else False, axis=1)
nprimer = len(primer_present[primer_present == True].index)
print("Number of reads containing the primer : %i" % nprimer)

adaptor_present = original_df.apply(lambda x: True if x['adaptor_start'] > 0 else False, axis=1)
nadaptor = len(adaptor_present[adaptor_present == True].index)
print("Number of reads containing the adaptor : %i" % nadaptor)

both_present = original_df.apply(lambda x: True if ((x['primer_start'] == 1) and (x['adaptor_start'] > 0)) else False, axis=1)
nboth = len(both_present[both_present == True].index)
print("Number of reads containing both primer and adaptor : %i" % nboth)

#Replace 0 with NA
original_df['primer_start'] = original_df['primer_start'].replace(0, 'NA')
original_df['primer_end'] = original_df['primer_end'].replace(0, 'NA')
original_df['adaptor_start'] = original_df['adaptor_start'].replace(0, 'NA')
original_df['adaptor_end'] = original_df['adaptor_end'].replace(0, 'NA')

#Save dataframe to csv file
original_df.to_csv(r'start_end.csv')


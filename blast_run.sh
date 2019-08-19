#!/bin/bash

blast_path="/usr/local/ncbi/blast/bin"
#Make blast database with primer and adaptor sequences
$blast_path/makeblastdb -in test.fna -title testdb -dbtype nucl -parse_seqids -out testdb

#Run nucleotide blast and produce output in tabular format
#m8 is legacy option in blast for tabular output. In BLAST+, -outfmt 6 does the same job
$blast_path/blastn -query prim_adapt.fasta -db testdb -task "blastn-short" -outfmt 6 -out test_blast.tab


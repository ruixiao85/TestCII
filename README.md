# TestCII

You need to write a small program using a programming language of your choice that will take:

a) Fasta and quality files containing the reads and corresponding quality scores. You can download the files from the link below. Please let me know if there are problems in downloading files.

https://drive.google.com/file/d/0BzTpYi_vYHmldXZfbGRYNHkzaVE/edit?usp=sharing

b) Adaptor and primer sequences are:

Primer Sequence:CGCCGTTTCCCAGTAGGTCTC
Adaptor Sequence:ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG

The program should generate the following output:

1) Total number of reads in the dataset
2) Total number of reads greater than 100 bp
3) Total number of reads with average quality scores greater than 20
4) Total number of reads with primer sequences
5) Total number of reads with adaptor sequences
6) Total number of reads with both primer and adaptor sequences

In addition, your program needs to generate the following files:

1) Blast output file in m8 format.
2) Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
3) Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences.

Please tar and gzip your program along with a README file and email it. Please feel free to contact me if you have any questions.

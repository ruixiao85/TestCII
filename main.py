import os
import timeit
from statistics import mean

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

if __name__=='__main__':
	fna,qual,fastq="test.fna","test.qual","test.fastq"
	unit='bp'
	qc2=100
	qc3=20
	primer=Seq("CGCCGTTTCCCAGTAGGTCTC")
	adaptor=Seq("ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG")

	if not os.path.exists(fastq): # pair sequence and quality files into one fastq file
		with open(fna) as f_handle,open(qual) as q_handle:
			records=PairedFastaQualIterator(f_handle,q_handle)
			count=SeqIO.write(records,fastq,"fastq")

	fq=SeqIO.parse(fastq,"fastq")
	c1,c2,c3,c4,c5,c6,c7=0,0,0,0,0,0,0
	now=timeit.timeit()
	for rec in fq:
		c1+=1 # increment every entry # print(rec.seq)
		if len(rec)>qc2: c2+=1 # increment if length met predefined criteria
		scores=rec.letter_annotations["phred_quality"] # print(scores)
		ave_score=mean(scores) # print(ave_score)
		if ave_score>qc3: c3+=1 # increment if average quality score met predefined criteria
		primer_in=primer in rec.seq
		adaptor_in=adaptor in rec.seq
		if primer_in: c4+=1
		if adaptor_in: c5+=1
		if primer_in and adaptor_in: c6+=1

	print(f'Processed in {timeit.timeit()-now:.5f} seconds.')
	print(f'O1. Total number of matching reads: {c1}.')
	print(f'O2. Number of reads greater than {qc2} {unit}: {c2}')
	print(f'O3. Number of reads with average quality score greater than {qc3}: {c3}')
	print(f'O4. Number of reads with primer sequences: {c4}')
	print(f'O5. Number of reads with adapter sequences: {c5}')
	print(f'O6. Number of reads with both primer and adapter sequences: {c6}')


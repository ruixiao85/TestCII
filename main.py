import os
import timeit
from statistics import mean

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

if __name__=='__main__':
	primer=Seq("CGCCGTTTCCCAGTAGGTCTC")
	adaptor=Seq("ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG")
	fna,qual,fastq="test.fna","test.qual","test.fastq"
	unit='bp'
	qo2=100
	qo3=20

	if not os.path.exists(fastq):
		with open(fna) as f_handle,open(qual) as q_handle:
			records=PairedFastaQualIterator(f_handle,q_handle)
			count=SeqIO.write(records,fastq,"fastq")

	now=timeit.timeit()
	fq=SeqIO.parse(fastq,"fastq")
	o1,o2,o3,o4,o5,o6,o7=0,0,0,0,0,0,0
	for rec in fq:
		o1+=1 # print(rec.seq)
		if len(rec)>qo2: o2+=1
		scores=rec.letter_annotations["phred_quality"] # print(scores)
		mean_score=mean(scores) # print(mean_score)
		if mean_score>qo3: o3+=1
		primer_in=primer in rec.seq
		adaptor_in=adaptor in rec.seq
		if primer_in: o4+=1
		if adaptor_in: o5+=1
		if primer_in and adaptor_in: o6+=1

	print(f'time lapsed {timeit.timeit()-now} seconds.')
	print(f'O1. Total number of matching reads: {o1}.')
	print(f'O2. Number of reads greater than {qo2} {unit}: {o2}')
	print(f'O3. Number of reads with average quality score greater than {qo3}: {o3}')
	print(f'O4. Number of reads with primer sequences: {o4}')
	print(f'O5. Number of reads with adapter sequences: {o5}')
	print(f'O6. Number of reads with both primer and adapter sequences: {o6}')


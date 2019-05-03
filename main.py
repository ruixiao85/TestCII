import os

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

	fq=SeqIO.parse(fastq,"fastq")

	sizes=[len(rec) for rec in fq]
	print(f'Out 1. Total number of matching reads: {len(sizes)} ({min(sizes)} {unit} ~ {max(sizes)} {unit}).')
	print(f'Out 2. Number of reads greater than {qo2} {unit}: {sum(1 for r in sizes if r>qo2)}')

	scores=[s for r in fq for s in r.letter_annotations["phred_quality"]]
	print(scores)

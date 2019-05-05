import os
import timeit
from statistics import mean

from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

if __name__=='__main__':
	fna,qual,fastq="test.fna","test.qual","test.fastq"
	qc2,qc3,unit=100,20,'bp'
	primer=Seq("CGCCGTTTCCCAGTAGGTCTC")
	adaptor=Seq("ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG")
	blast_file_folder,blast_file_prefix="blast","test"
	if not os.path.exists(blast_file_folder): os.mkdir(blast_file_folder)
	blast_result_file="1.blast_m8.txt"
	filter_trim_file="2.filter_trim.fna"
	primer_adaptor_file="3.primer_adaptor_loc.txt"

	if not os.path.exists(fastq): # pair sequence and quality files into one fastq file
		with open(fna) as f_handle,open(qual) as q_handle:
			records=PairedFastaQualIterator(f_handle,q_handle)
			count=SeqIO.write(records,fastq,"fastq")

	fq=SeqIO.parse(fastq,"fastq")
	c1,c2,c3,c4,c5,c6,c7=0,0,0,0,0,0,0
	now=timeit.timeit()
	fb,cb=open(blast_result_file,'w'),0
	fb.write('\t'.join(["query","subject","%id","alignment_length","mismatches","gap_openings",
		"query_start","query_end","subject_start","subject_end","E_value","bit_score"])+'\n')
	for rec in fq:
		c1+=1 # increment every entry # print(rec.seq)
		if len(rec)>qc2: c2+=1 # increment if length met predefined criteria
		scores=rec.letter_annotations["phred_quality"] # print(scores)
		if mean(scores)>qc3: c3+=1 # increment if average quality score met predefined criteria
		primer_in=primer in rec.seq
		adaptor_in=adaptor in rec.seq
		if primer_in: c4+=1
		if adaptor_in: c5+=1
		if primer_in and adaptor_in:
			c6+=1

		blast_file=os.path.join(blast_file_folder,f'{blast_file_prefix}_{rec.name}.xml')
		if not os.path.exists(blast_file): # qblast ncbi for xml format if not found
			qblast=NCBIWWW.qblast("blastn","nt",rec.seq)
			with open(blast_file,"w") as file:
				file.write(qblast.read())
			print(f'blast completed for {rec.name}')
			qblast.close()
		if os.path.exists(blast_file): # if qblast and download successfully
			qblast=open(blast_file)
			for r in NCBIXML.parse(qblast):
				for alignment in r.alignments:
					for hsp in alignment.hsps:
						_al,_id,_gap=hsp.align_length,hsp.identities,hsp.gaps
						fields=[rec.name,alignment.hit_id,f'{_id/len(rec.seq):.1%}',f'{_al}',f'{_al-_id-_gap}',f'{_gap}',
							f'{hsp.query_start}',f'{hsp.query_end}',f'{hsp.sbjct_start}',f'{hsp.sbjct_end}',
							f'{hsp.expect:.1e}',f'{hsp.bits:.1e}']
						fb.write("\t".join(fields)+"\n")
				cb+=1 # count of sequences blasted and written down
		# if c1>3: break # debug only, quickly check the first few entries
	fb.close()
	print(f'Processed in {timeit.timeit()-now:.5f} seconds.')

	print('The program should generate the following output:')
	print(f'1.Total number of matching reads: {c1}.')
	print(f'2.Number of reads greater than {qc2} {unit}: {c2}')
	print(f'3.Number of reads with average quality score greater than {qc3}: {c3}')
	print(f'4.Number of reads with primer sequences: {c4}')
	print(f'5.Number of reads with adapter sequences: {c5}')
	print(f'6.Number of reads with both primer and adapter sequences: {c6}')

	print('In addition, your program needs to generate the following files:')
	print(f'1.Total number of reads blasted and written into the m8 format file [{blast_result_file}]: {cb} ({cb/c1:.1%}).')
	print(f'2.Fasta file containing reads greater than {qc2} {unit}, average read quality scores greater than {qc3}, primers and adaptors trimmed (assuming "and" condition) [{filter_trim_file}].')
	print(f'3.Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences [{primer_adaptor_file}].')



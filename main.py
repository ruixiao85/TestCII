import argparse
import os
import timeit
import datetime
from statistics import mean

from Bio import SeqIO,motifs
from Bio.Blast import NCBIWWW,NCBIXML
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

def find_print_trim(_rec,_seq,_seq_name,_writer):
	_len=len(_seq)
	_motifs=motifs.create([_seq]) # potentially include ,_seq.complement(),_seq.reverse_complement()
	_list=[]
	for pos,seq in _motifs.instances.search(_rec.seq):
		_writer.write('\t'.join([_rec.name,_seq_name,f'{pos}',f'{pos+_len}'])+'\n')
		_list.append(pos)
	if _list:
		_list.reverse() # reverse in place
		for pos in _list: # position move to the left or become smaller as we delete, so delete from the right side
			_rec=_rec[:pos]+_rec[pos+_len:] # delete in place # TODO consider updating description length=?

		return True,_rec # target sequence found
	return False,rec # target sequence not found

if __name__=='__main__':
	parser=argparse.ArgumentParser()
	parser.add_argument('-f','--file',dest='file',action='store',default='test',help='~.fna, ~.qual, ~.fastq')
	parser.add_argument('-l','--length',dest='length',action='store',default=100,help='keep reads with length greater than this value')
	parser.add_argument('-q','--quality',dest='quality',action='store',default=20,help='keep reads with average quality score greater than this value')
	parser.add_argument('-u','--unit',dest='unit',action='store',default='bp',help='basic unit')
	parser.add_argument('-p','--primer',dest='primer',action='store',default='CGCCGTTTCCCAGTAGGTCTC',help='primer sequence')
	parser.add_argument('-a','--adaptor',dest='adaptor',action='store',default='ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG',help='adaptor sequence')
	args=parser.parse_args()

	fna,qual,fastq=args.file+".fna",args.file+".qual",args.file+".fastq"
	qc2,qc3,unit=args.length,args.quality,args.unit
	primer=Seq(args.primer)
	adaptor=Seq(args.adaptor)
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
	cor_len,cor_qual,primer_in,adaptor_in,pos_list=None,None,None,None,None # avoid temporary variables in the loop
	now=timeit.timeit()
	fb,cb=open(blast_result_file,'w'),0
	fb.write('\t'.join(["query","subject","%id","alignment_length","mismatches","gap_openings",
		"query_start","query_end","subject_start","subject_end","E_value","bit_score"])+'\n')
	fft=open(filter_trim_file,'w') # fasta no header needed
	fpa=open(primer_adaptor_file,'w')
	fpa.write('\t'.join(["identifier","match","start","end"])+'\n')
	for rec in fq:
		c1+=1 # increment every entry # print(rec.seq)
		cor_len=len(rec)>qc2
		if cor_len: c2+=1 # increment if length met predefined criteria
		scores=rec.letter_annotations["phred_quality"] # print(scores)
		cor_qual=mean(scores)>qc3
		if cor_qual: c3+=1 # increment if average quality score met predefined criteria
		primer_in,rec=find_print_trim(rec,primer,"primer",fpa) # return whether primer found and trimmed sequence
		adaptor_in,rec=find_print_trim(rec,adaptor,"adaptor",fpa) # return whether adaptor found and trimmed sequence
		if primer_in:
			c4+=1
		if adaptor_in: c5+=1
		if primer_in and adaptor_in: c6+=1

		if cor_len and cor_qual:
			SeqIO.write(rec,fft,"fasta")
		blast_file=os.path.join(blast_file_folder,f'{blast_file_prefix}_{rec.name}.xml')
		if not os.path.exists(blast_file): # qblast ncbi for xml format if not found
			qblast=NCBIWWW.qblast("blastn","nt",rec.seq)
			with open(blast_file,"w") as file:
				file.write(qblast.read())
			print(f'blast completed for {rec.name} @ {datetime.datetime.now()}')
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
		# if c1>4000: break
	fb.close()
	fft.close()
	fpa.close()
	print(f'Processed in {timeit.timeit()-now:.5f} seconds.')
	print()
	print('The program should generate the following output:')
	print(f'1. Total number of matching reads: {c1}.')
	print(f'2. Number of reads greater than {qc2} {unit}: {c2}')
	print(f'3. Number of reads with average quality score greater than {qc3}: {c3}')
	print(f'4. Number of reads with primer sequences: {c4}')
	print(f'5. Number of reads with adapter sequences: {c5}')
	print(f'6. Number of reads with both primer and adapter sequences: {c6}')
	print()
	print('In addition, your program needs to generate the following files:')
	print(f'1. Total number of reads blasted and written into the m8 format file "{blast_result_file}"": {cb} ({cb/c1:.1%}).')
	print(f'2. Fasta file containing reads greater than {qc2} {unit}, average read quality scores greater than {qc3}, primers and adaptors trimmed (assuming and/&& condition) "{filter_trim_file}".')
	print(f'3. Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences "{primer_adaptor_file}".')



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
   _motifs=motifs.create([_seq])  # potentially include _seq.complement(),_seq.reverse_complement()
   _list=[]
   for pos,seq in _motifs.instances.search(_rec.seq):
      _writer.write('\t'.join([_rec.name,_seq_name,f'{pos}',f'{pos+_len}'])+'\n')
      _list.append(pos)
   if _list:
      _list.reverse()  # reverse in place
      for pos in _list:  # position move to the left or become smaller as we delete, so delete from the right side
         _rec=_rec[:pos]+_rec[pos+_len:]  # delete in place # TODO consider updating description length=?
      return True,_rec  # target sequence found
   return False,_rec  # target sequence not found


if __name__=='__main__':
   parser=argparse.ArgumentParser()
   parser.add_argument('-f','--filename',dest='filename',action='store',default='test',
                       help='file name without extension, ~.fna, ~.qual, ~.fastq')
   parser.add_argument('-l','--length',dest='length',action='store',default=100,
                       help='keep reads with length greater than this value')
   parser.add_argument('-q','--quality',dest='quality',action='store',default=20,
                       help='keep reads with average quality score greater than this value')
   parser.add_argument('-u','--unit',dest='unit',action='store',default='bp',help='unit, default: bp')
   parser.add_argument('-p','--primer',dest='primer',action='store',default='CGCCGTTTCCCAGTAGGTCTC')
   parser.add_argument('-a','--adaptor',dest='adaptor',action='store',default='ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG')
   parser.add_argument('-w','--web',dest='web',action='store',default=False,
                       help='whether to query ncbi database for blast results')
   args=parser.parse_args()

   # initialize variables and create folders if needed
   basefile=args.filename
   fna,qual,fastq=basefile+".fna",basefile+".qual",basefile+".fastq"
   qc2,qc3,unit=args.length,args.quality,args.unit
   primer=Seq(args.primer)
   adaptor=Seq(args.adaptor)
   web_access=args.web # set to false to only process existing blast results
   blast_folder="blast"
   if not os.path.exists(blast_folder): os.mkdir(blast_folder)
   blast_result_file="1.blast_m8.txt"
   filter_trim_file="2.filter_trim.fna"
   primer_adaptor_file="3.primer_adaptor_loc.txt"

   # merge fna and qual files and write into one fastq file is not found, otherwise directly parse fastq
   if not os.path.exists(fastq):  # pair sequence and quality files into one fastq file, skip if available
      with open(fna) as f_handle,open(qual) as q_handle:
         records=PairedFastaQualIterator(f_handle,q_handle)
         count=SeqIO.write(records,fastq,"fastq")
      print(f'{count} entries was written to {fastq}.')
   fq=SeqIO.parse(fastq,"fastq") # once the fastq is generated, this step directly parse the fastq file

   # set counters to zeros, initialize filewriters
   c1,c2,c3,c4,c5,c6,c7=0,0,0,0,0,0,0
   fb,cb=open(blast_result_file,'w'),0
   fb.write('\t'.join(["query","subject","%id","alignment_length","mismatches","gap_openings","query_start",
                       "query_end","subject_start","subject_end","E_value","bit_score"])+'\n')
   fft=open(filter_trim_file,'w')  # fasta no header needed
   fpa=open(primer_adaptor_file,'w')
   fpa.write('\t'.join(["identifier","match","start","end"])+'\n')

   now=timeit.default_timer() # start timer
   for rec in fq: # loop only once for all console and file outputs to be efficient
      c1+=1  # increment every entry # print(rec.seq)
      cor_len=len(rec)>qc2 # whether has the correct length
      if cor_len: c2+=1  # increment if length met predefined criteria
      scores=rec.letter_annotations["phred_quality"]  # print(scores)
      cor_qual=mean(scores)>qc3 # whether has the correct average quality score
      if cor_qual: c3+=1  # increment if average quality score met predefined criteria
      has_primer,rec=find_print_trim(rec,primer,"primer",fpa)  # return whether has primer and trimmed sequence
      has_adaptor,rec=find_print_trim(rec,adaptor,"adaptor",fpa)  # return whether has adaptor and trimmed sequence
      if has_primer: c4+=1 # increment if has primer
      if has_adaptor: c5+=1 # increment if has adaptor
      if has_primer and has_adaptor: c6+=1  # increment if has both

      if cor_len and cor_qual: # if length and quality meet the criteria
         SeqIO.write(rec,fft,"fasta") # write trimmed sequence to fasta file
      blast_file=os.path.join(blast_folder,f'{basefile}_{rec.name}.xml') # path to local blast file (xml format)
      if not os.path.exists(blast_file) and web_access: # if not found and use internet
         qblast=NCBIWWW.qblast("blastn","nt",rec.seq) # perform ncbi qblast
         with open(blast_file,"w") as file:  # write blast result to local computer
            file.write(qblast.read())
         print(f'blast completed for {rec.name} @ {datetime.datetime.now()}.') # update console of the progress
         qblast.close()
      if os.path.exists(blast_file):  # if blast result found
         qblast=open(blast_file)
         for r in NCBIXML.parse(qblast): # iterate records in xml file
            for alignment in r.alignments: # iterate alignments in a record
               for hsp in alignment.hsps: # iterate high score pairs in a alignment
                  _al,_id,_gap=hsp.align_length,hsp.identities,hsp.gaps
                  fields=[rec.name,alignment.hit_id,f'{_id/len(rec.seq):.1%}',f'{_al}',f'{_al-_id-_gap}',f'{_gap}',
                          f'{hsp.query_start}',f'{hsp.query_end}',f'{hsp.sbjct_start}',f'{hsp.sbjct_end}',
                          f'{hsp.expect:.1e}',f'{hsp.bits:.1e}']
                  fb.write("\t".join(fields)+"\n")
         cb+=1  # count of blast file processed
   print(f'Processed in {timeit.default_timer()-now:.2f} seconds.') # time the loop

   fb.close()
   fft.close()
   fpa.close()
   print()
   print('The program should generate the following output:')
   print(f'1. Total number of matching reads: {c1}.')
   print(f'2. Number of reads greater than {qc2} {unit}: {c2}.')
   print(f'3. Number of reads with average quality score greater than {qc3}: {c3}.')
   print(f'4. Number of reads with primer sequences: {c4}.')
   print(f'5. Number of reads with adapter sequences: {c5}.')
   print(f'6. Number of reads with both primer and adapter sequences: {c6}.')
   print()
   print('In addition, your program needs to generate the following files:')
   print(f'1. Total number of reads blasted and written into the m8 format file":\n'
         f'   [{blast_result_file}] completed {cb} / {c1} ({cb/c1:.1%}) via NCBI qblast.')
   print(f'2. Fasta file containing reads greater than {qc2} {unit}, average read quality scores greater than {qc3},\n'
         f'   primers and adaptors trimmed.\n'
         f'   [{filter_trim_file}] (assuming and/&& condition, and trim all instances of exact match).')
   print(f'3. Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences:\n'
         f'   [{primer_adaptor_file}].')

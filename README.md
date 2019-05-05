# Test CII
---
## Question
You need to write a small program using a programming language of your choice that will take:

a) Fasta and quality files containing the reads and corresponding quality scores. You can download the files from the link below.
Please let me know if there are problems in downloading files.  
https://drive.google.com/file/d/0BzTpYi_vYHmldXZfbGRYNHkzaVE/edit?usp=sharing

b) Adaptor and primer sequences are
+ Primer Sequence: CGCCGTTTCCCAGTAGGTCTC  
+ Adaptor Sequence: ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG

The program should generate the following output:  
1. Total number of reads in the dataset
2. Total number of reads greater than 100 bp
3. Total number of reads with average quality scores greater than 20
4. Total number of reads with primer sequences
5. Total number of reads with adaptor sequences
6. Total number of reads with both primer and adaptor sequences

In addition, your program needs to generate the following files:
1. Blast output file in m8 format.
2. Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
3. Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences.

Please tar and gzip your program along with a README file and email it. Please feel free to contact me if you have any questions.

---
## Answer
Since python was recommended to solve the problem, I searched and found *biopython* to be most suitable for this task.

Noticeably, the primer and adaptor sequences and criteria for filtering reads can be made into arguments and become reusable and adjustable in future tasks.
>     parser=argparse.ArgumentParser()  
>     parser.add_argument('-l','--length',dest='length',action='store',default=100)  
>     parser.add_argument('-q','--quality',dest='quality',action='store',default=20)  
>     parser.add_argument('-p','--primer',dest='primer',action='store',
>                         default='CGCCGTTTCCCAGTAGGTCTC')  
>     parser.add_argument('-a','--adaptor',dest='adaptor',action='store',
>                         default='ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG')  
>     ...  
>     args=parser.parse_args()

Because combining fasta and quality files are computationally intensive, one-time process by pairing and writting to fastq file would be a good approach for future uses.
PairedFastaQualIterator in biopython can combine and write to a fastq file with the following code:
>     with open(fna) as f_handle, open(qual) as q_handle:  
>        records=PairedFastaQualIterator(f_handle,q_handle)  
>        count=SeqIO.write(records,fastq,"fastq")
   
Since all the console and file outputs can be accomplished in one pass of all sequence records, I will iterate through all records only once for optimal performance.
For each record, I will increment a counter (c1~c6) if the criteria are met.
>     c1+=1  # increment every entry  
>     cor_len=len(rec)>qc2 # whether has the correct length  
>     if cor_len: c2+=1  # increment if length met predefined criteria  
>     scores=rec.letter_annotations["phred_quality"]  
>     cor_qual=mean(scores)>qc3 # whether has the correct average quality score  
>     if cor_qual: c3+=1  # increment if average quality score met predefined criteria  
>     has_primer,rec=find_print_trim(rec,primer,"primer",fpa)  # return whether has primer and trimmed record 
>     has_adaptor,rec=find_print_trim(rec,adaptor,"adaptor",fpa)  # return whether has adaptor and trimmed record  
>     if has_primer: c4+=1 # increment if has primer  
>     if has_adaptor: c5+=1 # increment if has adaptor  
>     if has_primer and has_adaptor: c6+=1  # increment if has both  

Here the find_print_trim() function takes in SeqRecord, Seq(primer/adaptor), SeqName("primer"/"adaptor"), and FileWriter object.
It returns whether the SeqOfInterest is found in the SeqRecord and the SeqRecord after trimming the SeqOfInterest.
Here I made use of motif searching in the biopython package and reverse list of positions for sequence deletion.
>     _motifs=motifs.create([_seq])  # potentially include _seq.complement(),_seq.reverse_complement()  
>     _list=[]  
>     for pos,seq in _motifs.instances.search(_rec.seq):  
>        _writer.write('\t'.join([_rec.name,_seq_name,f'{pos}',f'{pos+_len}'])+'\n')  
>        _list.append(pos)  
>     if _list:  
>        _list.reverse()  # reverse in place  
>        for pos in _list:  # positions move to the left/smaller as we delete, so delete from the right side  
>           _rec=_rec[:pos]+_rec[pos+_len:]  # delete in place # TODO consider updating description length=?  

With the code block above, it also continuously writes to the *3rd* required file output.

For the *1st* required file, blast result can be retrieved either by query the ncbi database or perform blast on local computer.
I chose the former solution for this task because it's quicker to set up and the database are curated
without having to download and maintain the database as I do not run these tasks on a daily basis yet.

Retrieval of ncbi blast result in default xml format was done with the following code:
>     qblast=NCBIWWW.qblast("blastn","nt",rec.seq) # perform ncbi qblast  
>     with open(blast_file,"w") as file:  # write blast result to local computer  
>        file.write(qblast.read())

After the default xml files were retrieved, the files are parsed and saved into m8 format and the output can be highly customizable.
>     for r in NCBIXML.parse(qblast): # iterate records in xml file
>         for alignment in r.alignments: # iterate alignments in a record
>             for hsp in alignment.hsps: # iterate high score pairs in a alignment
>                 _al,_id,_gap=hsp.align_length,hsp.identities,hsp.gaps
>                 fields=[rec.name,alignment.hit_id,f'{_id/len(rec.seq):.1%}',f'{_al}',f'{_al-_id-_gap}',f'{_gap}',
>                         f'{hsp.query_start}',f'{hsp.query_end}',f'{hsp.sbjct_start}',f'{hsp.sbjct_end}',
>                         f'{hsp.expect:.1e}',f'{hsp.bits:.1e}']
>                 fb.write("\t".join(fields)+"\n")
 
For the *2nd* required file, since we already trimmed primers and adaptors in the find_print_trim() function, we can directly write SeqRecord to fasta file. 
>     if cor_len and cor_qual: # if length and quality meet the criteria
>        SeqIO.write(rec,fft,"fasta") # write trimmed sequence to fasta file
         
The final console output is as follows:
>     print('The program should generate the following output:')
>     print(f'1. Total number of matching reads: {c1}.')
>     print(f'2. Number of reads greater than {qc2} {unit}: {c2}.')
>     print(f'3. Number of reads with average quality score greater than {qc3}: {c3}.')
>     print(f'4. Number of reads with primer sequences: {c4}.')
>     print(f'5. Number of reads with adapter sequences: {c5}.')
>     print(f'6. Number of reads with both primer and adapter sequences: {c6}.')
>     print()
>     print('In addition, your program needs to generate the following files:')
>     print(f'1. Total number of reads blasted and written into the m8 format file":\n'
>           f'   [{blast_result_file}] completed {cb} / {c1} ({cb/c1:.1%}) via NCBI qblast.')
>     print(f'2. Fasta file containing reads greater than {qc2} {unit}, average read quality scores greater than {qc3},\n'
>           f'   primers and adaptors trimmed.\n'
>           f'   [{filter_trim_file}] (assuming and/&& condition, and trim all instances of exact match).')
>     print(f'3. Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences:\n'
>           f'   [{primer_adaptor_file}].')

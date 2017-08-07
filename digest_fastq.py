# Michael Olvera 2017
# This code was written to take in FASTQ reads, perform in-silico
# digest, and output the fragmented FASTQ files, to be used as part of
# a TLA analysis. 

"""
As per a TLA pipeline, initally unmapped reads can be remapped following an 
in-silico digest. This script takes in a .fastq file (exported from unmapped 
.sam reads) and outputs a .fastq file with reads digested by an input RE site.

Example:
	python3 digest_fastq.py tla_unmapped.fastq ATGC 4 tla_unmapped.digested -e -v

Inputs:
	input_fastq : A fastq file to be digested.

	digest_seq : The recognition site of where to make the cut. 

	cut_index : Where (relative to the first base of the digest site) 
	to make the cut. Index 0 is right befor ethe first base of the cut site,
	while index 1 is in between the first and second base of the site. Cuts
	can be make after the digest_seq. 
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import argparse


def highlight_cut(digest_seq, cut_index):
	"""
	Highlights were the cut is going to be made, and prints it out to the screen.
	TODO: Make the output optional.
	"""
	if cut_index >= len(digest_seq):
		digest_seq = digest_seq + 'N'*(cut_index - len(digest_seq) + 1)
	print 'Cutting here:', digest_seq[:cut_index],'*',digest_seq[cut_index:]

def fastq_const(seq, id, name, desc, phred):
	"""
	Constructor for new FASTQ objects (techincally SeqRecord objects in SeqIO).
	"""
	record = SeqRecord(Seq(seq, 'single_letter_alphabet'), id=id , name=name, description=desc)
	dict.__setitem__(record._per_letter_annotations, "phred_quality", phred)
	return record

def cut_read(fastq_read, digest_seq, cut_index, min_length):
	"""
	Cuts the fastq_reads at a specific site.
	"""
	fastq_seq = fastq_read.seq._data
	fastq_phred = fastq_read._per_letter_annotations['phred_quality']
	cut_indexes = [0] + [cut.start()+cut_index for cut in re.finditer(digest_seq, fastq_seq)]
	if len(cut_indexes) == 0:
		return [fastq_read] # No cutting
	else:
		# Split seq and phred
		digest_products = [fastq_seq[i:j] for i,j in zip(cut_indexes, cut_indexes[1:]+[None])]
		split_phred = [fastq_phred[i:j] for i,j in zip(cut_indexes, cut_indexes[1:]+[None])]
		# List new FASTQ fragments
		fragments = [fastq_const(digest_products[frag_i], fastq_read.id +'_%i'%(frag_i), 
			fastq_read.name +'_%i'%(frag_i), fastq_read.description, 
			split_phred[frag_i]) for frag_i in range(len(digest_products)) if len(digest_products[frag_i]) >= min_length]

		return fragments

def digest_fastq(fastq, digest_seq, cut_index, exclude_non_digested, min_length):
	"""
	Parse through the FASTQ file, and digest reads. 
	"""
	fq_list = list(SeqIO.parse(fastq, "fastq"))
	digested_list = list()
	for read in fq_list:
		cut_fragments = cut_read(read, digest_seq, cut_index, min_length)
		if exclude_non_digested and len(cut_fragments) != 1:
			digested_list += cut_fragments 

	if args.verbose:
		print "%i total reads processed, %i reads returned." %(len(fq_list) ,len(digested_list))
	return digested_list


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Digest a fastq file with a specific enzyme.')
	parser.add_argument('input_fastq', type=str, help='Fastq file of expected library reads.')
	parser.add_argument('digest_seq', type=str, help='Sequence of the digest site')
	parser.add_argument('cut_index', type=int, help='Where in the sequence to cut. For example, \
		0 = cut right before the sequence, 1 = cut between the firs tand second nucleotide of the recognition sequence.')
	parser.add_argument('out', type=str, help='Name of output fastq, minus \'.fastq\'.')
	parser.add_argument('--min_length', type=int, help='Minimum length of returned FASTQ reads', default=70)

	parser.add_argument('-e','--exclude_non_digested', help="Exclude reads that don't digest", action="store_true")
	parser.add_argument('-v', '--verbose', help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	if args.verbose:
		highlight_cut(args.digest_seq, args.cut_index)
	SeqIO.write( digest_fastq(args.input_fastq, args.digest_seq, args.cut_index, args.exclude_non_digested, args.min_length), args.out+".fastq", "fastq")


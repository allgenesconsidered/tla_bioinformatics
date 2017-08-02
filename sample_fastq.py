from Bio import SeqIO
from numpy.random import choice
import argparse


def sample_fastq(fastq, n, out):

	fq_list = list(SeqIO.parse(fastq, "fastq"))
	subset_index = choice(fq_list,n,replace=False)
	SeqIO.write(subset_index, out+".fastq", "fastq")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Randomly subset n number of reads from a fastq file.')
	parser.add_argument(
		'input_fastq', type=str, help='Fastq file of expected library reads.')
	parser.add_argument(
		'n', type=int, help='Number of reads to take.')
	parser.add_argument(
		'out', type=str, help='name of output fastq, minus \'.fastq\'.')
	args = parser.parse_args()

	sample_fastq(args.input_fastq, args.n, args.out)
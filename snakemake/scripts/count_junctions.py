import argparse, HTSeq
from collections import defaultdict

def main(args):

	gene_strand = {}
	genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)
	for line in open(args.genes):
		cur = line.rstrip().split('\t')
		gene_id = cur[3]
		if cur[5] == '+':
			orientation = '-'
		else:
			orientation = '+'
		iv = HTSeq.GenomicInterval(cur[0], int(cur[1])-1, int(cur[2]), orientation)
		genes[iv] += gene_id
		gene_strand[gene_id] = cur[5]


	five_prime = {}
	three_prime = {} 

	for line in open(args.exons):
		cur = line.rstrip().split('\t')
		chrom = cur[0]
		left = int(cur[1])
		right = int(cur[2])
		strand = cur[5]
		gene_id = cur[3].split(';')[0]
		exon_rank = int(cur[3].split(';')[1][1:])

		if strand == '+':
			five_prime[(chrom, left)] = (gene_id, exon_rank)
			three_prime[(chrom, right)] = (gene_id, exon_rank)
		else:
			five_prime[(chrom, right)] = (gene_id, exon_rank)
			three_prime[(chrom, left)] = (gene_id, exon_rank)

	junctions = defaultdict(int)

	for alignment in HTSeq.BAM_Reader(args.alignment):
		candidate_junctions = set()
		if alignment !=None and alignment.pe_which == 'first' and alignment.aligned and alignment.aQual > 5:
			for i, cigop in enumerate(alignment.cigar):
				if cigop.type=='N' and alignment.cigar[i-1].type=='M' and alignment.cigar[i-1].size>3 and alignment.cigar[i+1].type =='M' and alignment.cigar[i+1].size>3:
					candidate_junctions.add((cigop.ref_iv.chrom, cigop.ref_iv.start-1, cigop.ref_iv.end+1, cigop.ref_iv.strand))
		for candidate in candidate_junctions:
			overlaps = set()
			for iv, val in genes[HTSeq.GenomicInterval(candidate[0], candidate[1], candidate[2], candidate[3])].steps():
				overlaps |= val
			if len(overlaps) == 1:
				junctions[(candidate[0], candidate[1], candidate[2], list(overlaps)[0])] += 1
	output = []

	'''
	Iterates through junctions and classifies them as either consitutive, cryptic, exon skipping, or unannotated
	'''
	for junction in junctions:
		chrom, left, right, gene = junction
		if gene_strand[gene] == '+':
			if (chrom, left + 1) in three_prime and (chrom, right) in five_prime:
				low = three_prime[(chrom, left + 1)][1]
				high = five_prime[(chrom, right)][1]
				if low+1 == high:
					output.append('%s;%d;%d;%s;%d;%d;Constituitive\t%d' % (chrom, left, right, gene, low, high, junctions[junction]))
				else:
					output.append('%s;%d;%d;%s;%d;%d;Skipping\t%d' % (chrom, left, right, gene, low, high, junctions[junction]))
			elif (chrom, left + 1) in three_prime:
				output.append('%s;%d;%d;%s;%d;NA;Alt_3\t%d' % (chrom, left, right, gene, three_prime[(chrom, left + 1)][1], junctions[junction]))
			elif (chrom, right) in five_prime:
				output.append('%s;%d;%d;%s;NA;%d;Alt_5\t%d' % (chrom, left, right, gene, five_prime[(chrom, right)][1], junctions[junction]))
			else:
				output.append('%s;%d;%d;%s;NA;NA;Unannotated\t%d' % (chrom, left, right, gene, junctions[junction]))
		else:
			if (chrom, left + 1) in five_prime and (chrom, right) in three_prime:
				high = five_prime[(chrom, left + 1)][1]
				low = three_prime[(chrom, right)][1]
				if low+1 == high:
					output.append('%s;%d;%d;%s;%d;%d;Constituitive\t%d' % (chrom, left, right, gene, low, high, junctions[junction]))
				else:
					output.append('%s;%d;%d;%s;%d;%d;Skipping\t%d' % (chrom, left, right, gene, low, high, junctions[junction]))
			elif (chrom, left + 1) in five_prime:
				output.append('%s;%d;%d;%s;%d;NA;Alt_5\t%d' % (chrom, left, right, gene, five_prime[(chrom, left + 1)][1], junctions[junction]))
			elif (chrom, right) in three_prime:
				output.append('%s;%d;%d;%s;NA;%d;Alt_3\t%d' % (chrom, left, right, gene, three_prime[(chrom, right)][1], junctions[junction]))
			else:
				output.append('%s;%d;%d;%s;NA;NA;Unannotated\t%d' % (chrom, left, right, gene, junctions[junction]))

	with open(args.output, 'w') as out:
		out.write('\n'.join(output))

def parseArguments():
	parser = argparse.ArgumentParser(prog="count_junctions_reannotated.py", description='', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-a', '--alignment', required=True, help=' Alignment file (.bam)', metavar='', dest='alignment')
	required.add_argument('-g', '--genes', required=True, help=' File containing gene ranges (.bed)', metavar='', dest='genes')
	required.add_argument('-e', '--exons', required=True, help=' File containing exon ranges (.bed)', metavar='', dest='exons')
	required.add_argument('-o', '--output', required=True, help=' Output file containing counts', metavar='', dest='output')
	return parser.parse_args()

args = parseArguments()
main(args)
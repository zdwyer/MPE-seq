import argparse, HTSeq
from collections import defaultdict

def main():
	
	intron_set = set()
	targets = HTSeq.GenomicArrayOfSets("auto", stranded=True)
	for line in open(snakemake.input[2]):
		fields = line.rstrip().split('\t')
		if fields[5] == '+':
			orientation = '-'
		else:
			orientation = '+'
		iv = HTSeq.GenomicInterval(fields[0], int(fields[1])-1, int(fields[2]), orientation)
		targets[iv] += fields[3]
		intron_set.add(fields[3])

	premature = defaultdict(int)
	mature = defaultdict(int)

	for alignment in HTSeq.BAM_Reader(snakemake.input[0]):

		if alignment.pe_which == 'first' and alignment.aligned and alignment.aQual > 5:

			introns = set()
			junctions = set()

			for i, cigop in enumerate(alignment.cigar):
				if cigop.type == 'M':
					for iv, val in targets[cigop.ref_iv].steps():
						introns |= val
				elif cigop.type == 'N' and alignment.cigar[i-1].type =='M' and alignment.cigar[i-1].size > 3 and alignment.cigar[i+1].type =='M' and alignment.cigar[i+1].size > 3:
					for iv, val in targets[cigop.ref_iv].steps():
						junctions |= val

			candidate_genes = set()
			for i in introns:
				candidate_genes.add(i)
			for i in junctions:
				candidate_genes.add(i)

			if len(candidate_genes) == 1:
				for i in junctions:
					mature[i] += 1
				for i in introns:
					premature[i] += 1

	with open(snakemake.output[0], 'w') as out:
		out.write('Intron\tMature\tPremature\n')
		for intron in sorted(intron_set):
			out.write('%s\t%d\t%d\n' % (intron, mature[intron], premature[intron]))

main()
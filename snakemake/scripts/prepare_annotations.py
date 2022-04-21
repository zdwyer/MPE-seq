import sys, argparse
from collections import defaultdict

def main(args):

	targets = read_targets(args.targets)
	genes, alias, exons, introns, utr_5, utr_3 = read_features(args.feature)
	write_genes(args.genes, genes)
	write_alias(args.alias, alias)
	write_exons(args.exons, exons, utr_5, utr_3)
	write_introns(args.introns, args.target_out, introns, targets)

def read_targets(infile):
	targets = set()
	for line in open(infile):
		cur = line.rstrip().split(';')
		targets.add((cur[0], int(cur[1])))
	return targets

def read_features(infile):

	genes = []
	alias = {}
	exons = defaultdict(list)
	introns = defaultdict(list)
	utr_5 = {}
	utr_3 = {}

	for line in open(infile):
		if line[0] != '#':
			cur = line.rstrip().split('\t')
			if cur[2] == 'gene':
				info = cur[8]
				ID_loc = info.find('ID=')
				ID = info[ID_loc+3: info.find(';', ID_loc+3)]
				alias_loc = info.find('Name=')
				if alias_loc == -1:
					alias[ID] = ID
				else:
					alias[ID] = info[alias_loc+5: info.find(';', alias_loc+5)]
				genes.append('%s\t%s\t%s\t%s\t.\t%s' % (cur[0], cur[3], cur[4], ID, cur[6]))
			elif cur[2] == 'CDS':
				info = cur[8]
				parent_loc = info.find('Parent=')
				parent = info[parent_loc+7:]
				if parent[-2:] == '.1':
					ID = parent[:-2]
					exon_loc = info.find('exon:')
					exon = info[exon_loc+5:info.find(';', exon_loc+5)]
					exons[ID].append((cur[0], cur[3], cur[4], ID, exon, cur[6]))
			elif cur[2] == 'intron':
				info = cur[8]
				parent_loc = info.find('Parent=')
				parent = info[parent_loc+7:]
				if parent[-2:] == '.1':
					ID = parent[:-2]
					intron_loc = info.find('intron:')
					intron = info[intron_loc+7:info.find(';', intron_loc+7)]
					introns[ID].append((cur[0], cur[3], cur[4], ID, exon, cur[6]))
			elif cur[2] == 'five_prime_UTR':
				info = cur[8]
				parent_loc = info.find('Parent=')
				parent = info[parent_loc+7:]
				if parent[-2:] == '.1':
					ID = parent[:-2]
					utr_5[ID] = (cur[3], cur[4])
			elif cur[2] == 'three_prime_UTR':
				info = cur[8]
				parent_loc = info.find('Parent=')
				parent = info[parent_loc+7:]
				if parent[-2:] == '.1':
					ID = parent[:-2]
					utr_3[ID] = (cur[3], cur[4])

	return genes, alias, exons, introns, utr_5, utr_3

def write_genes(outfile, genes):
	with open(outfile, 'w') as out:
		out.write('\n'.join(genes))

def write_alias(outfile, alias):
	aliases = []
	for gene in alias:
		aliases.append('%s\t%s' % (gene, alias[gene]))
	with open(outfile, 'w') as out:
		out.write('\n'.join(sorted(aliases)))

def write_exons(outfile, exons, utr_5, utr_3):
	final_exons = []
	for gene in exons:
		orientation = exons[gene][0][5]
		exons[gene].sort(key=lambda x: x[1])
		if orientation == '+':
			if gene in utr_5:
				left = utr_5[gene][0]
			else:
				left = exons[gene][0][1]
			if gene in utr_3:
				right = utr_3[gene][1]
			else:
				right = exons[gene][-1][2]
		else:
			if gene in utr_3:
				left = utr_3[gene][0]
			else:
				left = exons[gene][0][1]
			if gene in utr_5:
				right = utr_5[gene][1]
			else:
				right = exons[gene][-1][2]
		if len(exons[gene]) == 1:
			final_exons.append('%s\t%s\t%s\t%s;E1\t.\t%s' % (exons[gene][0][0], left, right, gene, orientation))
		else:
			if orientation == '+':
				final_exons.append('%s\t%s\t%s\t%s;E1\t.\t%s' % (exons[gene][0][0], left, exons[gene][0][2], gene, orientation))
				for i in range(2,len(exons[gene])):
					final_exons.append('%s\t%s\t%s\t%s;E%d\t.\t%s' % (exons[gene][0][0], exons[gene][i-1][1], exons[gene][i-1][2], gene, i, orientation))
				final_exons.append('%s\t%s\t%s\t%s;E%d\t.\t%s' % (exons[gene][0][0], exons[gene][-1][1], right, gene, len(exons[gene]), orientation))
			else:
				final_exons.append('%s\t%s\t%s\t%s;E%d\t.\t%s' % (exons[gene][0][0], left, exons[gene][0][2], gene, len(exons[gene]), orientation))
				for i in range(2,len(exons[gene])):
					final_exons.append('%s\t%s\t%s\t%s;E%d\t.\t%s' % (exons[gene][0][0], exons[gene][i-1][1], exons[gene][i-1][2], gene, len(exons[gene]) - i + 1, orientation))
				final_exons.append('%s\t%s\t%s\t%s;E1\t.\t%s' % (exons[gene][0][0], exons[gene][-1][1], right, gene, orientation))
	with open(outfile, 'w') as out:
		out.write('\n'.join(final_exons))	

def write_introns(all_outfile, target_outfile, introns, targets):
	final_introns = []
	target_introns = []
	for gene in introns:
		introns[gene].sort(key=lambda x: x[1])
		orientation = introns[gene][0][5]
		for i in range(0, len(introns[gene])):
			if orientation == '+':
				final_introns.append('%s\t%s\t%s\t%s;I%d\t.\t%s' % (introns[gene][0][0], introns[gene][i][1], introns[gene][i][2], gene, i+1, orientation))
				if (gene, i+1) in targets:
					target_introns.append('%s\t%s\t%s\t%s;I%d\t.\t%s' % (introns[gene][0][0], introns[gene][i][1], introns[gene][i][2], gene, i+1, orientation))
			else:
				final_introns.append('%s\t%s\t%s\t%s;I%d\t.\t%s' % (introns[gene][0][0], introns[gene][i][1], introns[gene][i][2], gene, len(introns[gene]) - i, orientation))
				if (gene, len(introns[gene]) - i) in targets:
					target_introns.append('%s\t%s\t%s\t%s;I%d\t.\t%s' % (introns[gene][0][0], introns[gene][i][1], introns[gene][i][2], gene, len(introns[gene]) - i, orientation))
	with open(all_outfile, 'w') as out:
		out.write('\n'.join(final_introns))
	with open(target_outfile, 'w') as out:
		out.write('\n'.join(target_introns))

def parseArguments():
	parser = argparse.ArgumentParser(prog="", description='', usage='%(prog)s')
	input = parser.add_argument_group('input arguments')
	input.add_argument('-f', '--feature', required=True, help='', metavar='', dest='feature')
	input.add_argument('-t', '--targets', required=True, help='', metavar='', dest='targets')
	output = parser.add_argument_group('output arguments')
	output.add_argument('-g', '--genes', required=True, help='', metavar='', dest='genes')
	output.add_argument('-a', '--alias', required=True, help='', metavar='', dest='alias')
	output.add_argument('-e', '--exons', required=True, help='', metavar='', dest='exons')
	output.add_argument('-i', '--introns', required=True, help='', metavar='', dest='introns')
	output.add_argument('-u', '--target_introns', required=True, help='', metavar='', dest='target_out')
	return parser.parse_args()

args = parseArguments()
main(args)
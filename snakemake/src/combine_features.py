import os
from collections import defaultdict

def main():

	premature = defaultdict(int)
	mature = defaultdict(int)
	features = set() 
	samples = set()
	
	for infile in snakemake.input:
		sample = os.path.splitext(os.path.basename(infile))[0]
		samples.add(sample)
		for line in open(infile):
			cur = line.rstrip().split('\t')
			if cur[0] != "Intron":
				features.add(cur[0])
				mature[(sample, cur[0])] = int(cur[1])
				premature[(sample, cur[0])] = int(cur[2])

	with open(snakemake.output[0], 'w') as out:
		out.write('Strain\tReplicate\tIntron\tMature\tPremature\n')
		for sample in sorted(samples):
			info = sample.split('_')
			for feature in sorted(features):
				out.write('%s\t%s\t%s\t%d\t%d\n' % (info[0], info[1], feature, mature[(sample, feature)], premature[(sample, feature)]))

main()
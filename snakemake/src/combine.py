import os
from collections import defaultdict

def main():

	counts = defaultdict(int) # Repository to store counts per feature per sample. Key: (Sample, Feature), Value: Count
	features = set() # Holds the features encountered in all files.
	samples = set() # Holds all the samples run on.
	
	# Iterate through each line of each input file. Store the value in the second column of the input file in the count repository.
	
	for infile in snakemake.input:
		sample = os.path.splitext(os.path.basename(infile))[0]
		samples.add(sample)
		samples.add(sample)
		for line in open(infile):
			cur = line.rstrip().split('\t')
			features.add(cur[0])
			counts[(sample, cur[0])] = int(cur[1])

	# Output the combined counts for each file.
	with open(snakemake.output[0], 'w') as out:
		out.write('feature\t%s\n' % ('\t'.join(sample for sample in sorted(samples))))
		for feature in sorted(features):
			out.write('%s\t%s\n' % (feature, '\t'.join(str(counts[(sample, feature)]) for sample in sorted(samples))))

main()
import os, argparse
from collections import defaultdict

def main(args):

	counts = defaultdict(int) # Repository to store counts per feature per sample. Key: (Sample, Feature), Value: Count
	features = set() # Holds the features encountered in all files.
	samples = set() # Holds all the samples run on.
	
	# Iterate through each line of each input file. Store the value in the second column of the input file in the count repository.
	for infile in os.listdir(args.input_directory):
		sample = os.path.splitext(infile)[0]
		samples.add(sample)
		for line in open('%s%s' % (args.input_directory, infile)):
			cur = line.rstrip().split('\t')
			features.add(cur[0])
			counts[(sample, cur[0])] = int(cur[1])

	# Output the combined counts for each file.
	with open(args.output_file, 'w') as out:
		out.write('feature\t%s\n' % ('\t'.join(sample for sample in sorted(samples))))
		for feature in sorted(features):
			out.write('%s\t%s\n' % (feature, '\t'.join(str(counts[(sample, feature)]) for sample in sorted(samples))))

def parseArguments():
	parser = argparse.ArgumentParser(prog="combine.py", description='Combines all files in a directory into single file using the first column as a key and second column as a value.', usage='%(prog)s [options]')
	required = parser.add_argument_group('required arguments')
	required.add_argument('-d', '--directory', required=True, help=' Path to directory of files to combine.', metavar='', dest='input_directory')
	required.add_argument('-o', '--output', required=True, help=' Output file.', metavar='', dest='output_file')
	return parser.parse_args()

args = parseArguments()
main(args)
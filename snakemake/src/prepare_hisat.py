import sys

exons = []
for line in open(snakemake.input[0]):
    cur = line.rstrip().split('\t')
    exons.append('%s\t%d\t%d\t%s' % (cur[0], int(cur[1])-1, int(cur[2])-1, cur[5]))

introns = []
for line in open(snakemake.input[1]):
    cur = line.rstrip().split('\t')
    introns.append('%s\t%d\t%d\t%s' % (cur[0], int(cur[1])-2, int(cur[2]), cur[5]))

with open(snakemake.output[0], 'w') as out:
    out.write('\n'.join(exons))

with open(snakemake.output[1], 'w') as out:
    out.write('\n'.join(introns))

import sys

exons = []
for line in open(sys.argv[1]):
    cur = line.rstrip().split('\t')
    exons.append('%s\t%d\t%d\t%s' % (cur[0], int(cur[1])-1, int(cur[2])-1, cur[5]))

introns = []
for line in open(sys.argv[2]):
    cur = line.rstrip().split('\t')
    introns.append('%s\t%d\t%d\t%s' % (cur[0], int(cur[1])-2, int(cur[2]), cur[5]))

with open(sys.argv[3], 'w') as out:
    out.write('\n'.join(exons))

with open(sys.argv[4], 'w') as out:
    out.write('\n'.join(introns))

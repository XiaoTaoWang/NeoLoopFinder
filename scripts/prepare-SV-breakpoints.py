import sys

infil = sys.argv[1]
outfil = sys.argv[2]

data = [line.rstrip().split() for line in open(infil)]
writeout = []
with open(outfil, 'w') as out:
    for d in data:
        if d[-1]!='10kb':
            continue
        chrom1 = d[1]
        strand1 = d[4]
        chrom2 = d[5]
        strand2 = d[8]
        if strand1=='+':
            pos1 = int(d[3])
        else:
            pos1 = int(d[2])
        
        if strand2=='+':
            pos2 = int(d[7])
        else:
            pos2 = int(d[6])
        
        label = 'unknown'
        if chrom1 == chrom2:
            if pos1 > pos2:
                pos1, pos2 = pos2, pos1
                strand1, strand2 = strand2, strand1
            if strand1=='+' and strand2=='-':
                label = 'deletion'
            elif strand1=='-' and strand2=='+':
                label = 'duplication'
            else:
                label = 'inversion'
        else:
            if chrom1 > chrom2:
                chrom1, chrom2 = chrom2, chrom1
                strand1, strand2 = strand2, strand1
                pos1, pos2 = pos2, pos1
            label = 'translocation'
        
        new = [chrom1, chrom2, strand1+strand2, str(pos1), str(pos2), label]
        writeout.append(new)

with open(outfil, 'w') as out:
    for line in writeout:
        out.write('\t'.join(line)+'\n')

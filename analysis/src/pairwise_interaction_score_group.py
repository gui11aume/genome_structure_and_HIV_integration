import numpy as np
import cooler
import sys

if len(sys.argv) != 6:
    print "usage: {} ranges1.txt ranges2.txt cool_file.cool hic_bin_size out_file".format(sys.argv[0])
    sys.exit(1)
    
# Reps 3 arguments, 2 arxius amb format que conte els gens (sample, background) i el hi-c binsize
r1_file   = sys.argv[1]
r2_file   = sys.argv[2]
cool_file = sys.argv[3]
bin_size  = int(sys.argv[4])
out_file  = sys.argv[5]

# Load Cooler file.
cf = cooler.Cooler(cool_file)

# File format: 1st line is header
# Columns: gene.name, chr, beg, end, param1, param2
range1 = {}
range2 = {}
pcs  = []

with open(r1_file) as fg:
    header = fg.readline()
    for line in fg:
        m = line.rstrip().split('\t')
        if not range1.has_key(m[1]):
            range1[m[1]] = []
        range1[m[1]].append((m[0],int(m[2])/bin_size,int(m[3])/bin_size,float(m[4]),m[5]))

with open(r2_file) as fg:
    header = fg.readline()
    for line in fg:
        m = line.rstrip().split('\t')
        if not range2.has_key(m[1]):
            range2[m[1]] = []
        range2[m[1]].append((m[0],int(m[2])/bin_size,int(m[3])/bin_size,float(m[4]),m[5]))

# Compute all-pairs signal
print 'computing pairwise contact scores..'
for a in xrange(0,len(range1)):
    chr_a  = range1.keys()[a]
    print "{}".format(chr_a)
    for b in xrange(0,len(range2)):
        chr_b  = range2.keys()[b]
        if chr_a == chr_b: continue
        # Load Hi-C interchromosomal matrix.
        m = cf.matrix(balance=False).fetch(chr_a,chr_b)
        # Compute pair-wise contact score.
        for [name_a,beg_a,end_a,p1_a,p2_a] in range1[chr_a]:
            for [name_b,beg_b,end_b,p1_b,p2_b] in range2[chr_b]:
                # Pair-wise contact score (contacts/kb^2)
                reads   = np.sum(m[beg_a:(end_a+1),:][:,beg_b:(end_b+1)])
                pcscore = float(reads)/((end_a-beg_a+1)*(end_b-beg_b+1))/(float(bin_size*bin_size)/1e6)
                pcs.append((name_a,p1_a,p2_a,name_b,p1_b,p2_b,pcscore,reads))

# Print output in .hcs file
with open(out_file,'w') as fo:
    for score in pcs:
        fo.write('\t'.join([str(x) for x in score]) + '\n')

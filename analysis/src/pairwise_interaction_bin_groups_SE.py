import numpy as np
import cooler
import sys

if len(sys.argv) != 5:
    print "usage: {} ranges1.txt cool_file.cool hic_bin_size out_file".format(sys.argv[0])
    sys.exit(1)

# We receive a single file with the following columns:
#    1. Chromosome
#    2. beg
#    3. end
#    4. hiv cnt
#    5. Group (active/silent)

# Diccionari per cromosoma
    # Diccionari per actiu
    # Diccionari per VIH/no VIH
    
# Reps 3 arguments, 2 arxius amb format que conte els gens (sample, background) i el hi-c binsize
r1_file   = sys.argv[1]
cool_file = sys.argv[2]
bin_size  = int(sys.argv[3])
out_file  = sys.argv[4]

# Load Cooler file.
cf = cooler.Cooler(cool_file)

# File format: 1st line is header
# Columns: gene.name, chr, beg, end, param1, param2
ranges = {}
pcs = []
act_value = ['No SE','No SE','SE','SE']
hiv_value = ['No HIV','HIV','No HIV','HIV']

print 'reading data..'
with open(r1_file) as fg:
    header = fg.readline()
    for line in fg:
        m = line.rstrip().split('\t')
        if not ranges.has_key(m[0]):
            ranges[m[0]] = [[],[],[],[]]
        ranges[m[0]][(1 if m[4] == 'SE' else 0)*2 + (1 if int(m[3]) > 0 else 0)].extend(range(int(m[1])/bin_size, int(m[2])/bin_size+1))

# Compute all-pairs signal
print 'computing pairwise contact scores..'
for a in xrange(0,len(ranges)):
    chr_a  = ranges.keys()[a]
    sys.stdout.write("{}\t".format(chr_a))
    for b in xrange(a+1,len(ranges)):
        chr_b  = ranges.keys()[b]
        sys.stdout.write("{},".format(chr_b))
        sys.stdout.flush()
        # Load Hi-C interchromosomal matrix.
        m = cf.matrix(balance=False).fetch(chr_a,chr_b)
        # Compute pair-wise contact scores.
        sys.stdout.write(" (")
        for i in xrange(0,4):
            for j in xrange(i,4):
                if not ranges[chr_a][i] or not ranges[chr_b][j]:
                    continue
                vec_a = np.array(ranges[chr_a][i])
                vec_b = np.array(ranges[chr_b][j])
                # Pair-wise contact score (contacts/kb^2)
                reads   = np.sum(m[vec_a,:][:,vec_b])
                pcscore = float(reads)/(len(vec_a)*len(vec_b))/(float(bin_size*bin_size)/1e6)
                if chr_a == 'chrX' or chr_b == 'chrX':
                    pcscore *= 2
                pcs.append((chr_a,act_value[i],hiv_value[i],len(vec_a),chr_b,act_value[j],hiv_value[j],len(vec_b),pcscore,reads))
    sys.stdout.write("\n".format(chr_b))

# Print output in .hcs file
with open(out_file,'w') as fo:
    for score in pcs:
        fo.write('\t'.join([str(x) for x in score]) + '\n')

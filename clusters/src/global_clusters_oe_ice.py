import hic_basics
import cooler
import os

cf = cooler.Cooler(os.sys.argv[1])
clu_file = os.sys.argv[2]

hic_basics.interchromosomal_clusters(cf=cf, k=5, cluster_file=clu_file, out_file='global_clusters.txt')

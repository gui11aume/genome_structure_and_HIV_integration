import hic_basics
import cooler
import os

cf = cooler.Cooler(os.sys.argv[1])

hic_basics.cluster_compartments(cf=cf, k=15, chrlist=cf.chromnames, eig_dim=30, use_oe=True, use_ice=True, out_allchr='clusters_all.txt')

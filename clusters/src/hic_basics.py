from scipy import sparse
import itertools
import scipy
import scipy.linalg
import scipy.cluster
import numpy as np
import cooler
import matplotlib
import matplotlib.pyplot as plt
import os.path
from sklearn import mixture, manifold
from scipy.cluster.hierarchy import linkage,fcluster
from sklearn.cluster import SpectralClustering, KMeans
import pdb, traceback, sys

def full_matrix(chrlist, cf):
    pad_width = 25
    mat = None
    for c in chrlist:
        # Pad horizontal spacer
        if mat is not None:
            hpad = np.zeros((pad_width, mat.shape[1]))
            mat = np.concatenate([mat, hpad], axis=0)
            
        # Get first matrix and vertical spacer
        vmat = cf.matrix(balance=False).fetch(c,chrlist[0])
        vpad = np.zeros((vmat.shape[0], pad_width))
        
        # Get other matrix combinations
        for chr in chrlist[1:]:
            cm = cf.matrix(balance=False).fetch(c,chr)
            vmat = np.concatenate([vmat,vpad],axis=1)
            vmat = np.concatenate([vmat,cm],axis=1)

        # Concatenate matrices
        if mat is None:
            mat = vmat
        else:
            mat = np.concatenate([mat,vmat], axis=0)
            
    return mat
            
            
def scale_resolution(matrix, factor):
    new_row = matrix.row // factor
    new_col = matrix.col // factor
    bin = (matrix.shape[0] // factor) * new_col + new_row
    unique, bin = np.unique(bin, return_inverse=True)
    sum = np.bincount(bin, weights=matrix.data)
    new_col = unique // (matrix.shape[0] // factor)
    new_row = unique - new_col * (matrix.shape[0] // factor)

    return scipy.sparse.coo_matrix((sum, (new_row, new_col)))

def center_diag(mat, offset):
    d = np.diagonal(mat,offset)
    if (np.isnan(d) == True).all():
        return np.nan, np.nan
    td = d[~np.isnan(d)]
    if (td == 0).all():
        return np.nan, np.nan
    mean_val = np.mean(td)
    td -= mean_val
    span = np.max(td) - np.min(td)
    td[np.where(td > 0)] /= np.max(td)
    td[np.where(td < 0)] /= np.abs(np.min(td))
    d.setflags(write=True)
    d[~np.isnan(d)] = td
    return mean_val, span

def interchromosomal_clusters(cf,k,cluster_file,algorithm='eigh-kmeans',interchr_mat=None,use_builtin_ice=False,use_ice=False,use_cpb=False,mat_save_dir='.',out_file='global_clusters.out'):
    if algorithm not in ['eigh-gmix','eigh-kmeans','spec-kmeans']:
        print "error: algorithm must be either 'eigh-gmix', 'eigh-kmeans' or 'spec-kmeans'"
        
    clusters = {}
    print "[interchromosomal_clusters] k={}, cluster_file={}, out_file={}, algorithm={}".format(k,cluster_file,out_file,algorithm)
    # Read and parse intrachromosomal clusters.
    fc = open(cluster_file,"r")
    for line in fc:
        m = line.rstrip().split('\t')
        chr = m[0]
        clust = m[1].split(',')
        if not clusters.has_key(chr):
            clusters[chr] = []
        clusters[chr].append(np.array([int(v) for v in clust]))

    # Compute interchromosomal sums.
    if interchr_mat is None:
        print "computing interchromosomal matrix..."
        mat_shape = 0
        offset = {}
        for chr in clusters:
            offset[chr] = mat_shape
            mat_shape += len(clusters[chr])

        mat    = np.zeros((mat_shape,mat_shape))
        normat = np.zeros((mat_shape,mat_shape))

        # Iterate over all interchromosomal combinations
        for i in xrange(0,len(clusters)):
            chr_i = clusters.keys()[i]
            ploid_i = (1 if chr_i != 'chrX' else 2)
            sys.stdout.write("{} -> ".format(chr_i));
            for j in xrange(i+1,len(clusters)):
                chr_j = clusters.keys()[j]
                sys.stdout.write("{},".format(chr_j))
                ploid_j = (1 if chr_j != 'chrX' else 2)
                hic_m = cf.matrix(balance=use_builtin_ice).fetch(chr_i,chr_j)
                if use_ice:
                   hic_m = inter_ice(hic_m,20)
                elif use_cpb:
                   hic_m = inter_cpb(hic_m,20)
                   
                # Cluster combinations.
                for ki,clu_i in enumerate(clusters[chr_i]):
                    for kj,clu_j in enumerate(clusters[chr_j]):
                        sumv = np.sum(hic_m[clu_i,:][:,clu_j]) * ploid_i * ploid_j
                        mat[offset[chr_i]+ki,offset[chr_j]+kj] = mat[offset[chr_j]+kj,offset[chr_i]+ki] = sumv
                        normat[offset[chr_i]+ki,offset[chr_j]+kj] = normat[offset[chr_j]+kj,offset[chr_i]+ki] = sumv/float(len(clu_i))/float(len(clu_j))
            # Newline.
            sys.stdout.write('\n')
        # Store matrices for future processing
        np.save('{}/interchr_normat.npy'.format(mat_save_dir),normat)
        np.save('{}/interchr_sums.npy'.format(mat_save_dir),mat)
        
    else:
        print "Using user-provided interchromosomal matrix for cluster computation."
        normat = interchr_mat
        
    N = normat.shape[0]
    print "computing clusters, algorithm {}...".format(algorithm)
    if algorithm == 'spec-kmeans':
        spect_clu = SpectralClustering(n_clusters=k, eigen_solver='arpack', affinity='rbf', assign_labels='kmeans', n_jobs=8)
        hic_clust = spect_clu.fit_predict(normat)
    else:
        w, v = scipy.linalg.eigh(normat, eigvals=(N-k,N-1))
        if algorithm == 'eigh-gmix':
            gmix = mixture.GaussianMixture(n_components=k, covariance_type='full', tol=1e-4, max_iter=1000)
            gmix.fit(v)
            hic_clust = gmix.predict(v)
        elif algorithm == 'eigh-kmeans':
            km = KMeans(n_clusters=k,n_jobs=8)
            hic_clust = km.fit_predict(np.sqrt(w)*v)
            with open('{}'.format(out_file+'.weig'),'w+') as outdata:
                for (c,ev) in zip(hic_clust,w*v):
                    outdata.write(str(c)+'\t'+'\t'.join([str(x) for x in ev[::-1]]) + '\n')


    clu_idx = np.argsort(hic_clust)
    P = np.zeros(normat.shape)
    P[np.arange(0,len(clu_idx)),clu_idx] = 1
    # Permute rows and columns.
    W_clust = np.dot(np.dot(P,normat),np.linalg.inv(P))
    plt.matshow(W_clust,cmap=plt.cm.bwr)

    clust_cnt = [(g[0], len(list(g[1]))) for g in itertools.groupby(sorted(hic_clust))]
    # Compute cluster limits
    cnt = np.zeros(k+1,dtype=int)
    for i in xrange(k):
        cnt[i] += clust_cnt[i][1]
    cmcnt = np.cumsum(cnt)
    l = W_clust.shape[0]-1
    for x in cmcnt:
        plt.plot([0,l], [x,x], color='k', linestyle='-', linewidth=1)
        plt.plot([x,x], [0,l], color='k', linestyle='-', linewidth=1)

    print "writing results to {}...".format(out_file)
    with open(out_file,'w+') as fo:
        c_id = 0
        for i in xrange(0,len(clusters)):
            chr_i = clusters.keys()[i]
            for ki,clu_i in enumerate(clusters[chr_i]):
                fo.write("{}\t{}\t{}\t{}\n".format(clusters.keys()[i],ki,hic_clust[c_id],','.join([str(n) for n in clu_i])))
                c_id += 1
    print "[interchromosomal_clusters] Done."

def cluster_compartments(cf,k,chrlist,eig_dim=None,contact_thr=1,max_sample_size=50000,outlier_pctl=90,corr_outlier_pctl=[5,95],balance_corr_median=False,rm_diags=0,coeffs=None,coeffs_gw=None,seed=None,max_resampling_attempts=10,rearrange_clusters=False,use_oe=True,use_builtin_ice=False,use_ice=False,use_cpb=False,compute_abscore=True,algorithm='eigh-kmeans',outdir='.',out_allchr='clusters_all.txt', return_chr19=True):
    if algorithm not in ['eigh-gmix','eigh-kmeans','spec-kmeans']:
        print "error: algorithm must be either 'eigh-gmix', 'eigh-kmeans' or 'spec-kmeans'"
        return

    print "[intrachromosomal_clusters] k={}, outdir={}, algorithm={}".format(k,outdir,algorithm)

    if eig_dim == None:
        eig_dim = k

    if coeffs == None:
       coeffs = {}

    clusters = {}
    sample_idx = {}
    clusters_idx = {}
    snapshot_mat = None
    snapshot_cor = None

    for chr in chrlist:
        if os.path.isfile('{}/clusters_{}.txt'.format(outdir,chr)):
            print "Warning: {} clusters ({}/clusters_{}.txt) already exist. Skipping chromosome.".format(chr,outdir,chr)
            continue
        print "[{}] balancing matrix...".format(chr)
        
        m = cf.matrix(balance=use_builtin_ice).fetch(chr)

        # Threshold contacts.
        m[np.where(m<contact_thr)] = 0

        print "[{}] computing O/E...".format(chr)
        if use_oe is True:
            if coeffs_gw is not None:
                m_oe = oe_apply(m,coeffs_gw).toarray()
            else:
                if not coeffs.has_key(chr):
                    coeffs[chr] = oe_coeffs_mask(cf,chr,use_builtin_ice=use_builtin_ice)[chr]
                m_oe = oe_apply(m,coeffs[chr]).toarray()
        else:
            m_oe = m

        # Remove diagonals
        if rm_diags > 0:
            m_oe = np.triu(m_oe, k=rm_diags) + np.tril(m_oe,k=-rm_diags)

        # Get idx of high quality regions (measured in raw matrix).
        samp_idx = matrix_mask_idx(m_oe)
        sample_idx[chr] = samp_idx
        print "[{}] removing low-quality regions (matrix rows: {}, sample rows: {})...".format(chr,m.shape[0],samp_idx.shape[0])
        # High-quality matrix size
        l = len(samp_idx)
        ssize = min(l,max_sample_size)
        # Sample iteration (keep sampling while clustering fails).
        np.random.seed(seed)
        successful = False
        cnt = 0
        while not successful and cnt < max_resampling_attempts:
            cnt += 1
            # Get sample
            if ssize < l:
                s = np.sort(np.random.choice(samp_idx,ssize,replace=False))
            else:
                s = np.array(samp_idx)
            m_samp = m_oe[s,:][:,s]

            # Apply CPB.
            if use_cpb:
               print "[{}] computing CPB...".format(chr)
               m_samp = cpb(m_samp,20)
            elif use_ice:
               print "[{}] computing ICE...".format(chr)
               m_samp = ice(m_samp,20)

            # Relax outliers
            m_max = np.percentile(m_samp[np.where(m_samp>0)],outlier_pctl)
            m_samp[np.where(m_samp > m_max)] = m_max
            if (~m_samp.any(axis=1)).any():
                print "[{}] sample contains empty rows (singular matrix). resampling ({})...".format(chr,cnt)
                continue
                
            # Compute correlation and remove diagonal
            print "[{}] computing correlation matrix and balancing...".format(chr)
            m_cor = np.corrcoef(m_samp)
            np.fill_diagonal(m_cor,0)
            
            if return_chr19 and chr == "chr19":
               snapshot_mat = np.matrix(m_samp)
               snapshot_cor = np.matrix(m_cor)

            
            # Increase correlation contrast (5-95 percentiles by default)
            if balance_corr_median:
                m_cor = m_cor - np.median(m_cor[np.triu_indices(ssize,1)])
            min_cor_val = np.percentile(m_cor[np.triu_indices(ssize,1)],corr_outlier_pctl[0])
            max_cor_val = np.percentile(m_cor[np.triu_indices(ssize,1)],corr_outlier_pctl[1])
            m_cor[np.where(m_cor < min_cor_val)] = min_cor_val
            m_cor[np.where(m_cor > max_cor_val)] = max_cor_val

            N = m_cor.shape[0]
            eig_dim = min(N,eig_dim)

            if compute_abscore:
                   print "[{}] computing a/b scores...".format(chr)
                   eigv, abscore = eigs_ab_score(m_cor)

                   # Write weighted eigenvectors
                   with open('{}/clusters_{}.abscore'.format(outdir,chr),'w') as outdata:
                      for i in xrange(0,len(abscore)):
                         outdata.write(str(sample_idx[chr][i])+'\t'+str(abscore[i])+'\n')

            try:
                print "[{}] computing clusters, algorithm {}...".format(chr,algorithm)
                if algorithm == 'spec-kmeans':
                    # some chromosomes crash when using precomputed similarity matrices.
                    # however using RBF seems to give meaningful clustering.
                    spect_clu = SpectralClustering(n_clusters=k, eigen_solver='arpack', affinity='precomputed', assign_labels='kmeans', n_jobs=8)
                    hic_clust = spect_clu.fit_predict(m_cor)
                else:
                    print "[{}] computing eigh...".format(chr)
                    w, v = scipy.linalg.eigh(m_cor, eigvals=(N-eig_dim,N-1))

                    if algorithm == 'eigh-gmix':
                        # Cluster eigenvectors using Gaussian Mixture
                        gmix = mixture.GaussianMixture(n_components=k, covariance_type='full', tol=1e-4, max_iter=1000)
                        gmix.fit(v)
                        hic_clust = gmix.predict(v)
                    elif algorithm == 'eigh-kmeans':
                        # Cluster eigenvalue/eigenvector products with kmeans.
                        print "[{}] computing clusters (k-means)...".format(chr)
                        km = KMeans(n_clusters=k,n_jobs=8)
                        weig = np.sqrt(w)*v
                        hic_clust = km.fit_predict(weig)
                        # Write weighted eigenvectors
                        with open('{}/clusters_{}.weig'.format(outdir,chr),'w') as outdata:
                            for i in xrange(0,len(hic_clust)):
                                outdata.write(str(sample_idx[chr][i])+'\t'+str(hic_clust[i])+'\t'+'\t'.join([str(x) for x in weig[i][::-1]])+'\n')
            except Exception, e:
                print "[{}] error while clustering: {}".format(chr,cnt,str(e))
                cnt = max_resampling_attempts
                break
            successful = True

        if cnt >= max_resampling_attempts:
            print "[{}] max reampling attempts reached, skipping chromosome.".format(chr)
            continue

        # Rearrange clusters for visualization
        # Make cluster index list
        clu_idx = [list() for _ in xrange(k)]
        for i,c in enumerate(hic_clust):
            clu_idx[c].append(i)

        if not rearrange_clusters:
            # Map again to matrix indices
            clusters_idx[chr] = [sample_idx[chr][x] for x in clu_idx]

        else:
            print "[{}] rearranging clusters by similarity...".format(chr)
            for i in xrange(k):
                clu_idx[i] = np.array(clu_idx[i])

            clusters[chr] = list()

            # Find most distant blocks
            l_r = (0,0)
            val = np.inf
            d_sum = np.zeros((k,k))
            for i in xrange(k):
                l_i = len(clu_idx[i])
                for j in xrange(i+1,k):
                    l_j = len(clu_idx[j])
                    d_sum[i,j] = np.sum(m_cor[clu_idx[i],:][:,clu_idx[j]])
                    d = float(d_sum[i,j])/(l_i*l_j)
                    if d < val:
                        l_r = (i,j)
                        val = d

            # Pop left and right blocks (important to do it in this order for index consistency).
            r_idx = clu_idx.pop(l_r[1])
            l_idx = clu_idx.pop(l_r[0])
            r_clusters = [r_idx.copy(),]
            l_clusters = [l_idx.copy(),]

            iters = len(clu_idx)/2 + len(clu_idx)%2
            for i in xrange(iters):
                # Find nearest blocks to L/R.
                len_l = len(l_idx)
                len_r = len(r_idx)
                min_d = np.inf
                max_d = -np.inf
                min_idx = 0
                max_idx = 0
                for i in xrange(len(clu_idx)):
                    len_block = len(clu_idx[i])
                    d_l = float(np.sum(m_cor[l_idx,:][:,clu_idx[i]]))/(len_l*len_block) - val
                    d_r = float(np.sum(m_cor[r_idx,:][:,clu_idx[i]]))/(len_r*len_block) - val
                    r = d_l/d_r
                    if r < min_d:
                        min_idx = i
                        min_d = r
                    if r >= max_d:
                        max_idx = i
                        max_d = r
                # Pop from idx and add to L/R.
                if min_idx > max_idx:
                    r_clusters.append(clu_idx[min_idx].copy())
                    l_clusters.append(clu_idx[max_idx].copy())
                    r_idx = np.append(clu_idx.pop(min_idx),r_idx)
                    l_idx = np.append(l_idx,clu_idx.pop(max_idx))
                elif min_idx < max_idx:
                    r_clusters.append(clu_idx[min_idx].copy())
                    l_clusters.append(clu_idx[max_idx].copy())
                    l_idx = np.append(l_idx,clu_idx.pop(max_idx))
                    r_idx = np.append(clu_idx.pop(min_idx),r_idx)
                else:
                    l_clusters.append(clu_idx[max_idx].copy())
                    l_idx = np.append(l_idx,clu_idx.pop(max_idx))
            # Make final index list.
            clu_idx = np.append(l_idx,r_idx)

            # Make final cluster index list.
            clusters[chr] = l_clusters + list(reversed(r_clusters))

            # Map again to matrix indices
            clusters_idx[chr] = [sample_idx[chr][x] for x in clusters[chr]]

        # Store in disk
        print "[{}] writing clusters to {}/clusters_{}.txt...".format(chr,outdir,chr)
        fout = open('{}/clusters_{}.txt'.format(outdir,chr),'w+')
        for c in clusters_idx[chr]:
            fout.write("{}\t".format(chr))
            fout.write(','.join([str(i) for i in c]))
            fout.write('\n')
        fout.close()
        fall = open('{}/{}'.format(outdir,out_allchr),"a")
        for c in clusters_idx[chr]:
            fall.write("{}\t".format(chr))
            fall.write(','.join([str(i) for i in c]))
            fall.write('\n')
        fall.close()

        '''
        # Make permutation matrix
        P = np.zeros(m_cor.shape)
        P[np.arange(0,len(clu_idx)),clu_idx] = 1

        # Permute rows and columns.
        W_clust = np.dot(np.dot(P,m_cor),np.linalg.inv(P))
        plt.matshow(W_clust,cmap=plt.cm.bwr)

        # Compute cluster limits
        cnt = np.zeros(k+1,dtype=int)
        for i in xrange(len(clusters[chr])):
            cnt[i] += len(clusters[chr][i])
        cnt = np.cumsum(cnt)
        l = W_clust.shape[0]-1
        for x in cnt:
            plot([0,l], [x,x], color='k', linestyle='-', linewidth=1)
            plot([x,x], [0,l], color='k', linestyle='-', linewidth=1)
        title('{} cluster k={}'.format(chr,k))
        '''
        
    print "[intrachromosomal_clusters] Done."

    return snapshot_mat, snapshot_cor, clusters_idx

        
def interchromosomal_cluster_matrix(cf,clusters_idx,c1,c2,use_builtin_ice=False):
    # Interchromosomal cluster matrix
    # clusters_idx['chr'] contains the cluster indices.
    # Read matrix
    m = cf.matrix(balance=use_builtin_ice).fetch(c1,c2)
    # Compute cluster lengths
    len_i = [len(x) for x in clusters_idx[c1]]
    len_j = [len(x) for x in clusters_idx[c2]]

    # Compute interchromosomal matrix with cluster delimiters
    allidx_i = np.array([idxset for clust in clusters_idx[c1] for idxset in clust])
    allidx_j = np.array([idxset for clust in clusters_idx[c2] for idxset in clust])
    m_ic = m[allidx_i,:][:,allidx_j]

    # Plot matrix with cluster delimiters
    plt.matshow(m_ic,cmap=plt.cm.Reds)

    cnt_i = np.cumsum(len_i)
    cnt_j = np.cumsum(len_j)
    m_rows = m_ic.shape[0]
    m_cols = m_ic.shape[1]
    for x in cnt_i:
        plot([0,m_cols], [x,x], color='k', linestyle='-', linewidth=1)
    for x in cnt_j:
        plot([x,x], [0,m_rows], color='k', linestyle='-', linewidth=1)

    # Compute interchromosomal cluster matrix
    clusts_i = len(clusters_idx[c1])
    clusts_j = len(clusters_idx[c2])
    icm = np.zeros((clusts_i,clusts_j))
    for i,idx_i in enumerate(clusters_idx[c1]):
        for j,idx_j in enumerate(clusters_idx[c2]):
            icm[i,j] = float(sum(m[idx_i,:][:,idx_j]))/(len(idx_i)*len(idx_j))

    plt.matshow(icm)

def inter_ice_step(m):
   rs = np.asarray(np.sum(m,1)).reshape(-1)
   cs = np.asarray(np.sum(m,0)).reshape(-1)
   rnz = rs[rs>0]
   cnz = cs[cs>0]
   rf = np.ones(m.shape[0])
   cf = np.ones(m.shape[1])
   rf[rs>0] = np.sqrt(rnz/np.mean(rnz))
   cf[cs>0] = np.sqrt(cnz/np.mean(cnz))
   m = m/np.outer(np.sqrt(rf),np.sqrt(cf))
   return m, np.var(np.append(rnz,cnz))/np.mean(np.append(rnz,cnz))

def inter_ice(m, n_iter):
   for i in xrange(0,n_iter):
      m,v = inter_cpb_step(m)
      print "iCPB variance[{}]: {}".format(i+1, v)
   return m

def inter_cpb_step(m):
   rs = np.asarray(np.sum(m,1)).reshape(-1)
   cs = np.asarray(np.sum(m,0)).reshape(-1)
   rnz = np.asarray(np.sum(m>0,1)).reshape(-1)[rs>0]
   cnz = np.asarray(np.sum(m>0,0)).reshape(-1)[cs>0]
   mrs = rs[rs>0]/rnz.astype('f')
   mcs = cs[cs>0]/cnz.astype('f')
   rf = np.ones(m.shape[0])
   cf = np.ones(m.shape[1])
   rf[rs>0] = mrs/np.mean(mrs)
   cf[cs>0] = mcs/np.mean(mcs)
   m = m/np.outer(np.sqrt(rf),np.sqrt(cf))
   return m, np.var(np.append(rs,cs))/np.mean(np.append(rs,cs))

def inter_cpb(m, n_iter):
   for i in xrange(0,n_iter):
      m,v = inter_cpb_step(m)
      print "iCPB variance[{}]: {}".format(i+1, v)
   frcs = np.append(np.sum(m,0),np.sum(m,1))
   print "iCPB final variance: {}".format(np.var(frcs)/np.mean(frcs))
   
   return m

def ice_step(m):
    # Matrix ice balancing.
    rs = np.asarray(np.sum(m,0)).reshape(-1)        # Total interaction signal
    rsnz = rs[rs>0]
    rf = np.ones(m.shape[0])
    rf[rs>0] = np.sqrt(rsnz/np.mean(rsnz))
    m = m/np.outer(rf,rf)
    return m, np.var(rs)/np.mean(rs)

def ice(m, n_iter):
   '''
    ICE algorithm.
    Returns ICE-balanced version of matrix after applying n_iter iterations.
    '''
   for i in xrange(1,n_iter+1):
      m,v = ice_step(m)
      print "ICE variance/mean[{}]: {}".format(i,v)
   print "ICE final variance/mean {}".format(np.var(np.sum(m,0))/np.mean(np.sum(m,0)))
       
   return m

def cpb_step(m):
   rs = np.asarray(np.sum(m,1)).reshape(-1)
   rnz = np.asarray(np.sum(m>0,1)).reshape(-1)[rs>0]
   mrs = rs[rs>0]/rnz.astype('f')
   rf = np.ones(m.shape[0])
   rf[rs>0] = np.sqrt(mrs/np.mean(mrs))
   m = m/np.outer(rf,rf)
   return m, np.var(rs)/np.mean(rs)

def cpb(m, n_iter):
   '''
    Contact Probability Balancing algorithm.
    Returns CPB-balanced version of matrix after applying n_iter iterations.
    '''
   for i in xrange(1,n_iter+1):
      m,v = cpb_step(m)
      print "CPB variance [{}]: {}".format(i, v)
   print "CPB final variance {}".format(np.var(np.sum(m,0))/np.mean(np.sum(m,0)))
        
   return m

def oe_apply(m, coeffs):
    '''
    Applies O/E normalization given a sparse/coo_matrix and previously computed coefficients.
    '''
    matrix = sparse.coo_matrix(m, dtype=float)
    for i in range(0, len(matrix.data)):
        matrix.data[i] = matrix.data[i] / coeffs[abs(matrix.col[i]-matrix.row[i])]

    return matrix

def matrix_mask_idx(m, bsm_lambda = lambda x: np.median(x)/2, cpb_lambda = lambda x: np.median(x)/2):
        # Mask low quality regions.
        rsm = np.sum(m,1)
        bsm = np.sum(np.sign(m),1)
        cpb = rsm/bsm

        bsm_idx = np.where(bsm >= bsm_lambda(bsm))[0]
        cpb_idx = np.where(cpb >= cpb_lambda(cpb))[0]

        # Mask low quality 
        return np.unique(np.sort(np.append(bsm_idx,cpb_idx)))

def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

def oe_coeffs_mask(cfile, ranges, use_builtin_ice=False, bsm_lambda = lambda x: np.median(x)/2, cpb_lambda = lambda x: np.median(x)/2):
    '''
    Returns a dictionary with O/E normalization coefficients for each chromosome.
    Before computing O/E, some regions of the matrix are masked and excluded.
    Matrix row/columns are masked following the lambda functions:
    bsm_lambda: excludes all rows with bin sum (bins with at least 1 contact) smaller than the value returned by bsm_labmda.
    cpb_lambda: excludes all rows with contact per bin (average contact per row) smaller than the value returned by cpb_labmda.
    Dictionary keys are chromosome names.
    '''
    if type(ranges) is not list:
        ranges = [ranges,]

    pr_int = {}
    # Compute observed interactions for selected chromosomes
    chrlist = []
    chrbins = []
    for cname in ranges:
        chrname = cname.split(':')[0]
        chrlist.append(chrname)
        chrbins.append(cfile.extent(str(chrname)))

    for c,bins in zip(chrlist,chrbins):
        # Compute observed interactions for each chromosomes
        m = cfile.matrix(balance=use_builtin_ice).fetch(c)

        # Mask low quality regions.
        rsm = sum(m,1)
        bsm = sum(np.sign(m),1)
        cpb = rsm/bsm

        bsm_idx = np.where(bsm < bsm_lambda(bsm))[0]
        cpb_idx = np.where(cpb < cpb_lambda(cpb))[0]

        # Mask low quality 
        idx = np.unique(np.sort(np.append(bsm_idx,cpb_idx)))
        m[idx,:] = 0
        m[:,idx] = 0

        # Convert matrix to sparse
        m = scipy.sparse.coo_matrix(m)

        # ObI: Observed Interactions
        # PrI: Probability of Interaction
        max_int_dist = m.shape[0]
        obs_int = [0]*(max_int_dist)
        pr_int[c] = np.zeros(max_int_dist,dtype=float)

        for i,j,v in zip(m.row, m.col, m.data):
            obs_int[abs(i-j)] += v

        for int_dist in range(0,max_int_dist):
            if not obs_int[int_dist]:
                pr_int[c][int_dist] = float('NaN')
                continue
            # PsI: Possible interactions
            possible_int = max_int_dist-int_dist
            # Correct masked rows
            tmp = np.append(idx[np.where(idx >= int_dist)], idx + int_dist)
            tmp = tmp[np.where(tmp < max_int_dist)]
            possible_int -= len(np.unique(tmp))
            # Probability:
            if possible_int != 0:
                pr_int[c][int_dist] = obs_int[int_dist] / float(possible_int)
            else:
                pr_int[c][int_dist] = float('NaN')

        # Interpolate NaN values.
        nans, x= nan_helper(pr_int[c])
        if len(nans) > 0:
            pr_int[c][nans]= np.interp(x(nans), x(~nans), pr_int[c][~nans])

    return pr_int

    
def oe_coeffs_chrom(cfile, ranges, use_builtin_ice=False):
    '''
    Returns a dictionary with O/E normalization coefficients for each chromosome.
    Dictionary keys are chromosome names.
    '''
    if type(ranges) is not list:
        ranges = [ranges,]
    pr_int = {}
    # Compute observed interactions for selected chromosomes
    chrlist = []
    chrbins = []
    for cname in ranges:
        chrname = cname.split(':')[0]
        chrlist.append(chrname)
        chrbins.append(cfile.extent(str(chrname)))
 
    for c,bins in zip(chrlist,chrbins):
        max_int_dist = bins[1] - bins[0]
        # ObI: Observed Interactions
        # PrI: Probability of Interaction
        obs_int = [0]*(max_int_dist+1)
        pr_int[c] = [0.0]*(max_int_dist+1)

        # Compute observed interactions for each chromosomes
        m = cfile.matrix(balance=use_builtin_ice, sparse=True).fetch(c)
        for i,j,v in zip(m.row, m.col, m.data):
            obs_int[abs(i-j)] += v

        for int_dist in range(0,max_int_dist+1):
            if not obs_int[int_dist]: continue
            # PsI: Possible interactions
            possible_int = max_int_dist-int_dist
            # Probability:
            pr_int[c][int_dist] = obs_int[int_dist] / float(possible_int)

    return pr_int

def oe_coeffs_gw(cfile, use_builtin_ice=False):
    '''
    Returns an array containing the O/E normalization coefficient for all
    interaction distances.
    '''
    # Compute observed interactions for selected chromosomes
    chrlist = []
    chrbins = []
    for cname in cfile.chroms()[:]['name']:
        chrname = cname.split(':')[0]
        chrlist.append(chrname)
        chrbins.append(cfile.extent(str(chrname)))
            
    # Max interaction distance
    chrsz = sorted([x[1]-x[0] for x in chrbins], reverse = 1)
    max_int_dist = chrsz[0]-1
    
    # ObI: Observed Interactions
    # PrI: Probability of Interaction
    obs_int = [0]*(max_int_dist+1)
    pr_int = [0.0]*(max_int_dist+1)

    # Compute observed interactions for all chromosomes
    for c in chrlist:
        m = cfile.matrix(balance=use_builtin_ice, sparse=True).fetch(c)
        for i,j,v in zip(m.row, m.col, m.data):
            obs_int[abs(i-j)] += v

    # Compute probability of interaction
    nchr = len(chrsz)
    cumsz = sum(chrsz)

    for int_dist in xrange(0,max_int_dist+1):
        if obs_int[int_dist] == 0: continue
        # Update chromosome count
        while chrsz[nchr-1] <= int_dist:
            nchr -= 1
            cumsz -= chrsz[nchr]
        # PsI: Possible interactions
        possible_int = cumsz-nchr*int_dist
        # Probability:
        pr_int[int_dist] = obs_int[int_dist] / float(possible_int)

    pr_int = np.array(pr_int)
    pr_int[np.where(pr_int == 0)] = np.min(pr_int[np.where(pr_int > 0)])
            
    return pr_int

def eigs_ab_score(cor_mat, ab_percentile=20):
    '''
    Inputs:
    matrix cor_mat : correlation matrix
    ab_percentile: percentile of reference A/B compartments
    Output:
    returns AB score for all bins in cor_mat
    '''
    N = cor_mat.shape[0]
    # Block eigenvector computation (of correlation matrix):
    ab_score = np.empty(N, dtype='f')
    ab_score[:] = 0
        
    # Compute eigenvectors
    w, v = scipy.linalg.eigh(cor_mat, eigvals=(N-1,N-1))
    v = v[:,np.argmax(w)]
    v = v/np.max(np.abs(v))

    # Compute A/B scores
    a_set = np.where(v > np.percentile(v,100-ab_percentile))[0]
    b_set = np.where(v < np.percentile(v,ab_percentile))[0]
    a_score = np.sum(cor_mat[:,a_set],axis=1)/float(len(a_set))
    b_score = np.sum(cor_mat[:,b_set],axis=1)/float(len(b_set))
    a_mod   = np.sum(np.abs(cor_mat[:,a_set]),axis=1)/float(len(a_set))
    b_mod   = np.sum(np.abs(cor_mat[:,b_set]),axis=1)/float(len(b_set))
    ab_score = (a_score-b_score)/(a_mod+b_mod)*100.0

    # Store data
    return v, ab_score

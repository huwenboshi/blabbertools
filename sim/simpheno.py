#!/usr/bin/python

"""
python script to simulate phenotypes
input data: haplotypes, legend
output data: simulated phenotypes

author: huwenbo shi
"""


from optparse import OptionParser
import sys
import math
import os
import random
import numpy as np

# get command line
parser = OptionParser(add_help_option=False)
parser.add_option("-h", "--hap", dest="hapfile")
parser.add_option("-m", "--marker", dest="markerfile")
parser.add_option("-f", "--freqth", dest="freqth_str")
parser.add_option("-n", "--nsnps", dest="nsnps_str")
parser.add_option("-p", "--pint", dest="pint_str")
parser.add_option("-q", "--ncau", dest="ncau_str")
parser.add_option("-o", "--outprefix")

(options, args) = parser.parse_args()

hapfile_nm = options.hapfile
markerfile_nm = options.markerfile
freqth_str = options.freqth_str
nsnps_str = options.nsnps_str
pint_str = options.pint_str
ncau_str = options.ncau_str
outprefix = options.outprefix

# check command line
if(hapfile_nm == None or markerfile_nm == None or
   freqth_str == None or nsnps_str == None or
   pint_str == None or outprefix == None or
   ncau_str == None):
    sys.stderr.write("Usage:\n")
    sys.stderr.write("\tUse -h to specify hap file\n")
    sys.stderr.write("\tUse -m to specify marker file\n")
    sys.stderr.write("\tUse -f to specify MAF threshold\n")
    sys.stderr.write("\tUse -n to specify the number of SNPs\n")
    sys.stderr.write("\tUse -p to specify the number of interactions\n")
    sys.stderr.write("\tUse -o to specify output prefix\n")
    sys.stderr.write("\tUse -q to specify number of causal SNPs\n")
    sys.exit()

print '#################################\n'
print '** simulation started **'
print '\n#################################\n'

# parse numeric values
freqth = float(freqth_str)
nsnps = int(nsnps_str)
pint = int(pint_str)
ncau = int(ncau_str)

# constants - heritability
h2 = 0.2
h2_int = 0.02
korder = 2

# read in legend file
markers = []
flr = False
markerfile = open(markerfile_nm, 'r')
for line in markerfile:
    if(flr == False):
        flr = True
        continue
    line = line.strip()
    markers.append(line)
markerfile.close()

# read in haplotype file
# each col of haps stores a haplotype
haps = []
hapfile = open(hapfile_nm, 'r')
for line in hapfile:
    line = line.strip()
    cols = line.split()
    for i in xrange(len(cols)):
        cols[i] = int(cols[i])
    haps.append(cols)
hapfile.close()
haps = np.matrix(haps)

# convert haps to gens
nrow = haps.shape[0]
ncol = haps.shape[1]
gens = np.zeros((nrow, ncol/2))
for i in xrange(nrow):
    for j in xrange(0, ncol, 2):
        gens[i,j/2] = haps[i,j]+haps[i,j+1]

# obtain allele frequency
# filter out snps with maf less than freqth
sim_snps = []
sim_snps_idx = []
freqs = gens.mean(1)/2
nrow = freqs.shape[0]
for i in xrange(nrow):
    if(len(sim_snps) >= nsnps):
        break
    if(freqs[i] > freqth and freqs[i] < 1-freqth):
        sim_snps.append(markers[i])
        sim_snps_idx.append(i)

# select the snps, haps, and freqs indexed in sim_snps_idx
# normalize snps
gens_norm = gens[sim_snps_idx,:]
gens = gens[sim_snps_idx,:]
haps = haps[sim_snps_idx,:]
freqs = freqs[sim_snps_idx]
nrow = gens_norm.shape[0]
ncol = gens_norm.shape[1]
for i in xrange(nrow):
    for j in xrange(ncol):
        var = 2*freqs[i]*(1-freqs[i])
        gens_norm[i,j] = (gens_norm[i,j]-2*freqs[i])/math.sqrt(var)

print '** normalized all snps **'
print '-- mean of gi'
print gens_norm.mean(1)
print '-- var of gi'
print np.var(gens_norm)
print '\n#################################\n'

# draw single snps
pool_size = gens_norm.shape[0]
idx_pool = range(pool_size)
cau_idx = []
for i in xrange(ncau):
    random.shuffle(idx_pool)
    cau_idx.append(idx_pool[0])
    del idx_pool[0]

print '** selected %d causal snps, %d interactions **' % (ncau, pint)
print '\n#################################\n'

# draw interaction snps
snp_idx_pairs = []
for i in xrange(pint):
    idx_pair = [-1]*korder
    for j in xrange(korder):
        random.shuffle(idx_pool)
        idx_pair[j] = idx_pool[0]
        del idx_pool[0]
    snp_idx_pairs.append(idx_pair)

# create interaction genotypes
nsnp = gens_norm.shape[0]
nind = gens_norm.shape[1]
nint = len(snp_idx_pairs)
gens_int = np.ones((nint, nind))
for i in xrange(nint):
    for j in xrange(nind):
        for k in xrange(korder):
            gens_int[i,j] = gens_int[i,j]*gens[snp_idx_pairs[i][k],j]

# normalize interaction genotypes
gens_int_norm = gens_int
for i in xrange(nint):
    idx0 = snp_idx_pairs[i][0]
    idx1 = snp_idx_pairs[i][1]
    gens0 = gens[idx0,:]
    gens1 = gens[idx1,:]
    mu0 = gens0.mean(0)
    mu1 = gens1.mean(0)
    cov01 = np.cov(gens0, gens1)[0,1]
    gens0_sq = np.square(gens0)
    gens1_sq = np.square(gens1)
    mu0_sq = gens0_sq.mean(0)
    mu1_sq = gens1_sq.mean(0)
    cov01_sq = np.cov(gens0_sq, gens1_sq)[0,1]
    mu = mu0*mu1+cov01
    var = mu0_sq*mu1_sq+cov01_sq-mu*mu
    for j in xrange(nind):
        gens_int_norm[i,j] = (gens_int_norm[i,j]-mu)/math.sqrt(var)

print '** normalized all interactions **'
print '-- mean of gi*gj'
print gens_int_norm.mean(1)
print '-- var of gi*gj'
print np.var(gens_int_norm)
print '\n#################################\n'

# simulate phenotypes
# initialization
pheno = np.zeros((nind, 1))
if(ncau > 0):
    betas = np.random.normal(0.0, math.sqrt(h2/ncau), ncau)
betas_int = np.random.normal(0.0, math.sqrt(h2_int/nint), nint)

# add single snp effect
gens_norm_cau = gens_norm[cau_idx,:]
for i in xrange(nind):
    for j in xrange(ncau):
        gen = gens_norm_cau[j,i]
        b = betas[j]
        pheno[i] = pheno[i]+b*gen

# add snp interaction effect
for i in xrange(nind):
    for j in xrange(nint):
        gen_int = gens_int_norm[j,i]
        b_int = betas_int[j]
        pheno[i] = pheno[i]+b_int*gen_int

# add environment effect
env_effect = np.random.normal(0.0, math.sqrt(1-h2-h2_int), (nind,1))
pheno = pheno+env_effect

print '** simulation ended **'
print '\n#################################\n'

# write out result
hapfile = open(outprefix+'.sim.hap', 'w')
nrow = haps.shape[0]
ncol = haps.shape[1]
for i in xrange(nrow):
    line = ''
    for j in xrange(ncol):
        line += str(haps[i,j])+' '
    hapfile.write(line+'\n')
hapfile.close()

snpfile = open(outprefix+'.sim.snp', 'w')
snpfile.write('snp_id pos x0 x1\n')
for m in markers:
    snpfile.write(m+'\n')
snpfile.close()

phefile = open(outprefix+'.sim.phe', 'w')
nind = pheno.shape[0]
for i in xrange(nind):
    phefile.write(str(pheno[i,0])+'\n')
phefile.close()

from StrainFinder import *
import cPickle
import sys
import numpy
import parse_midas_data
import parse_timecourse_data

species_name = sys.argv[1]
Nmax = 7
Lmax = 1000
output_directory = os.path.expanduser("~/strainfinder_output/")

filename_prefix = "%s/%s" % (output_directory, species_name)

# Get best (min) AIC version
# Get filenames of EM objects
fns = [filename_prefix+("%d.em.cpickle" % N) for N in range(2,Nmax+1)]    

# Load EM objects
ems = [cPickle.load(open(fn, 'rb')) for fn in fns]

# Get the best AIC in each EM object
aics = [em.select_best_estimates(1)[0].aic for em in ems]

# Select EM with the minimum AIC
#em = ems[numpy.argmin(aics)]
em = ems[1]

input_alignment = em.data.x # alignment data, dim = (M x L x 4)

best_estimate = em.select_best_estimates(1)[0]

strain_genotypes = best_estimate.p # genotypes of first estimate, dim = (N x L x 4)
strain_freqs = best_estimate.z # (M x N)

#output_alignment = numpy.dot(strain_freqs, strain_genotypes)

output_alignment = numpy.einsum('ij,jkl', strain_freqs, strain_genotypes)

print input_alignment.shape
print output_alignment.shape

print "Best # of strains:", em.select_best_estimates(1)[0].z.shape[1]
print em.select_best_estimates(1)[0].z # frequencies of first estimate, dim = (M x N)

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts = numpy.array([sample_time_map[sample] for sample in parse_timecourse_data.morteza_samples])

import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

pylab.figure(figsize=(7,6))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(3,1, height_ratios=[1,1,1], hspace=0.15)

strain_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(strain_axis)
for i in xrange(0,strain_freqs.shape[1]):
    strain_axis.plot(ts,strain_freqs[:,i])
strain_axis.set_ylim([0,1])

predicted_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(predicted_axis)

for i in xrange(0,output_alignment.shape[1]):
    
    allele_counts = output_alignment[:,i,:]
    depths = allele_counts.sum(axis=1)
    freqs = (allele_counts[:,0]/(depths+(depths==0)))
    
    good_idxs = (depths>0)
    
    if good_idxs.sum() < 2:
        continue
    
    masked_ts = ts[good_idxs]
    masked_freqs = freqs[good_idxs]
    
    if (masked_freqs[0])>0.5:
        masked_freqs = 1-masked_freqs

    predicted_axis.plot(masked_ts, masked_freqs,'-',color='0.7',alpha=0.5)

predicted_axis.set_ylim([0,1])

observed_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(observed_axis)
for i in xrange(0,input_alignment.shape[1]):
    
    allele_counts = input_alignment[:,i,:]
    depths = allele_counts.sum(axis=1)
    freqs = (allele_counts[:,0]/(depths+(depths==0)))
    
    good_idxs = (depths>0)
    
    if good_idxs.sum() < 2:
        continue
    
    masked_ts = ts[good_idxs]
    masked_freqs = freqs[good_idxs]
    
    if (masked_freqs[0])>0.5:
        masked_freqs = 1-masked_freqs

    observed_axis.plot(masked_ts, masked_freqs,'-',color='0.7',alpha=0.5)

observed_axis.set_ylim([0,1])

pylab.savefig('strainfinder_output.pdf',bbox_inches='tight')
###############################
#
# Rest of script begins here
#
################################

import pylab
import numpy
import sys
from math import log10
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from numpy.random import binomial, random_sample
import bz2
import parse_midas_data
import parse_timecourse_data
import matplotlib
import matplotlib.pyplot as plt
import timecourse_utils

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("settings_filename", help="settings file")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

settings_filename = args.settings_filename
debug = args.debug
chunk_size = args.chunk_size
################################################################################

sample_time_map = parse_timecourse_data.parse_sample_time_map()
theory_ts = numpy.array([t for t in sorted(set(sample_time_map.values()))])
theory_ts = theory_ts[theory_ts>0]




PLOT_FMAJOR=None
PLOT_FMINOR=None
additional_titles=None
COLORED_LINEWIDTH=0.5

PLOT_APPEARANCE_TIME=False

# Mininum coverage for frequency estimation vs interpolation 
min_coverage = 20


# load settings
settings_file = open(settings_filename,"r")
settings_string = "\n".join(settings_file.readlines())
settings_file.close()
exec settings_string    

if additional_titles==None:
    additional_titles = ["" for species_name in species_names]

mpl.rcParams['font.size'] = 7.0
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

fig_width = 7
fig_height = 1.7*(len(species_names))

fig, axes = plt.subplots(len(species_names),sharex=True,sharey=True,figsize=(fig_width, fig_height))

if len(species_names)<2:
    axes = [axes]


freq_axis = None
    
for species_idx in xrange(0,len(species_names)):        

    species_name = species_names[species_idx]
    
    
    sys.stderr.write("Processing %s...\n" % species_name)
    
    species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
    sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
    desired_samples = numpy.array(samples)[sample_idxs]
    
    times = []
    alt_matrix = []
    depth_matrix = []
    snp_infos = []

    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=desired_samples, initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        if len(samples)<5:
            break
     
    
        sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

        # Calculate fixation matrix
        sys.stderr.write("Calculating allele freqs...\n")
        chunk_alts, chunk_depths, chunk_snp_infos = timecourse_utils.calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
        sys.stderr.write("Done!\n")
    
        chunk_alts = chunk_alts[:,sample_idxs]
        chunk_depths = chunk_depths[:,sample_idxs]
        # polarize using first timepoint
        chunk_alts += (chunk_depths-2*chunk_alts)*(((chunk_alts[:,0]+chunk_alts[:,1])>((chunk_depths[:,0]+chunk_depths[:,1])/2))[:,None])
    
        desired_sites = ((chunk_alts>(0.1*chunk_depths)).sum(axis=1)>2)*((chunk_depths>0).sum(axis=1)>10)
    
        chunk_alts = chunk_alts[desired_sites,:]
        chunk_depths = chunk_depths[desired_sites,:]
        chunk_allele_freqs = chunk_alts*1.0/(chunk_depths+(chunk_depths==0))
    
        if len(times)==0:
            times = sample_ts
        
                      
        if desired_sites.sum()>0:
            alt_matrix.append(chunk_alts)
            depth_matrix.append(chunk_depths)
            desired_site_idxs = numpy.nonzero(desired_sites)[0]
            for idx in desired_site_idxs:
                snp_infos.append(chunk_snp_infos[idx])
                    
    sys.stderr.write("Done!\n")


    # set up figure axis
    freq_axis = axes[species_idx]
    
    if additional_titles[species_idx]=="":
        title_text = species_name
    else:
        title_text = '%s %s' % (species_name, additional_titles[species_idx])  
    freq_axis.set_title(title_text,loc='right',fontsize=6)
        
    freq_axis.set_ylabel('Allele frequency, $f(t)$')
    freq_axis.spines['top'].set_visible(False)
    freq_axis.spines['right'].set_visible(False)
    freq_axis.get_xaxis().tick_bottom()
    freq_axis.get_yaxis().tick_left()
    
    num_colored_mutations = 0
    num_total_mutations = 0

    
    if len(alt_matrix)==0:
        continue
                
    alt_matrix = numpy.vstack(alt_matrix)
    depth_matrix = numpy.vstack(depth_matrix) 
    
    p = min([1,1000.0/(len(snp_infos)+1)])

    for mutation_idx in xrange(0,len(snp_infos)):
        
        num_total_mutations += 1
     
        chromosome, location, gene_name, variant_type = snp_infos[mutation_idx]
        alts = alt_matrix[mutation_idx,:]
        depths = depth_matrix[mutation_idx,:]
        
        freqs = alts*1.0/(depths+(depths==0))
        
        masked_times = times[depths>0]
        masked_freqs = freqs[depths>0]
        
        
        if masked_freqs[0]>0.5:
            masked_freqs = 1-masked_freqs
            
        if (masked_freqs>0.1).sum() < 2:
            continue
         
        interpolation_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs)
        
        if color_condition(species_idx, chromosome, location, gene_name, variant_type, masked_freqs, interpolation_function):
        
            # One of the colored ones!
            num_colored_mutations+=1
            
            #sys.stderr.write("%s %d %s %s\n" % (gene_name, location, var_type, allele)) 
            
            line, = freq_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=COLORED_LINEWIDTH)
            color = pylab.getp(line,'color')
            
        else:  
            # One of the non-colored ones
            #freq_axis.plot(theory_times, interpolation_function(theory_times), '-', alpha=0.5, color='0.7', markersize=3,linewidth=1,zorder=1)
            if random_sample() < p:
                freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
     
    sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
    
freq_axis.set_xlabel('Time, $t$ (days)')

freq_axis.set_ylim([0,1.02])
freq_axis.set_xlim([0,160])   
 
sys.stderr.write("Saving final PNG image...\t")
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")


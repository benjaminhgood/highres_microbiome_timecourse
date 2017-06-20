###############################
#
# Rest of script begins here
#
################################
import matplotlib  
matplotlib.use('Agg') 
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
import parse_patric

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

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
species_idx_map = {species[i]: i for i in xrange(0,len(species))}
sample_time_map = parse_timecourse_data.parse_sample_time_map()

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]
species_freq_matrix = species_coverage_matrix*1.0/(species_coverage_matrix.sum(axis=0))

desired_samples = numpy.array(samples)[species_time_idxs]
    

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
fig_height = 1.7*(2*len(species_names))

fig, axes = plt.subplots(2*len(species_names),sharex=True,sharey=False,figsize=(fig_width, fig_height))    

#if len(species_names)<2:
#    axes = [axes]

freq_axis = None
    
for species_idx in xrange(0,len(species_names)):        

    species_name = species_names[species_idx]
    
    
    sys.stderr.write("Processing %s...\n" % species_name)
    
    
    antibiotic_resistance_genes = parse_patric.load_antibiotic_resistance_genes(species_name)
    virulence_factors = parse_patric.load_virulence_factors(species_name)
    core_genes = parse_timecourse_data.load_core_timecourse_genes(species_name, min_copynum=0.3, min_prevalence=0.9, min_marker_coverage=min_coverage)
    
    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=desired_samples)
    sys.stderr.write("Done!\n")
    times, marker_coverage_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)

    gene_names = numpy.array(gene_names)
    gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
    marker_coverages = marker_coverages[marker_coverage_idxs]
    gene_copynum_matrix = gene_copynum_matrix[:,marker_coverage_idxs]
    
    gene_idxs = ((gene_copynum_matrix>0.3).sum(axis=1)>0)
    gene_copynum_matrix = gene_copynum_matrix[gene_idxs,:]
    gene_names = gene_names[gene_idxs] 

    gene_copynum_matrix = numpy.clip(gene_copynum_matrix, 1e-02,100)

    # set up figure axis
    
    species_freq_axis = axes[2*species_idx]
    if additional_titles[species_idx]=="":
        title_text = '%s abundance' % species_name
    else:
        title_text = '%s %s coverage' % (species_name, additional_titles[species_idx])  
    species_freq_axis.set_title(title_text,loc='right',fontsize=6)
        
    species_freq_axis.set_ylabel('Species abundance')
    species_freq_axis.spines['top'].set_visible(False)
    species_freq_axis.get_xaxis().tick_bottom()
    
    depth_axis = species_freq_axis.twinx()
    depth_axis.set_ylabel('Marker coverage',color='#007ccd',rotation=270,labelpad=10)
    for tl in depth_axis.get_yticklabels():
        tl.set_color('#007ccd')
    depth_axis.spines['top'].set_visible(False)
    depth_axis.spines['right'].set_color('#007ccd')
    species_freq_axis.spines['right'].set_color('#007ccd')
    depth_axis.tick_params(axis='y', colors='#007ccd')
    depth_axis.tick_params(axis='y', which='minor', colors='#007ccd')
    
    depth_axis.set_ylim([2,2e03])
    depth_axis.set_xlim([0,160])   
    
    depth_axis.set_zorder(0)
    species_freq_axis.set_zorder(1)  
        
    
    copynum_axis = axes[2*species_idx+1]
    
    if additional_titles[species_idx]=="":
        title_text = '%s genes' % species_name
    else:
        title_text = '%s %s genes' % (species_name, additional_titles[species_idx])  
    copynum_axis.set_title(title_text,loc='right',fontsize=6)
        
    copynum_axis.set_ylabel('Gene copynum')
    copynum_axis.spines['top'].set_visible(False)
    copynum_axis.spines['right'].set_visible(False)
    copynum_axis.get_xaxis().tick_bottom()
    copynum_axis.get_yaxis().tick_left()
    
    copynum_axis.set_xlim([0,160])   
    copynum_axis.set_ylim([1e-01,100])
    
    copynum_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[1e-01,1e-01],[100,100],color=parse_timecourse_data.antibiotics_color,linewidth=0)
    copynum_axis.fill_between([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[1e-01,1e-01],[100,100],color=parse_timecourse_data.lyme_color,linewidth=0)
    copynum_axis.plot([parse_timecourse_data.hrv_infection, parse_timecourse_data.hrv_infection], [1e-01,100], 'k-',linewidth=0.25,zorder=-1)
    
    depth_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[2,2],[2e03,2e03],color=parse_timecourse_data.antibiotics_color,linewidth=0,zorder=-1)
    depth_axis.fill_between([parse_timecourse_data.lyme_infection, parse_timecourse_data.antibiotic_start],[2,2],[2e03,2e03],color=parse_timecourse_data.lyme_color,linewidth=0,zorder=-1)
    depth_axis.plot([parse_timecourse_data.hrv_infection, parse_timecourse_data.hrv_infection], [2,2e03], 'k-',linewidth=0.25,zorder=-1)
    
    
    num_colored_mutations = 0
    num_total_mutations = 0

    # Plot species abundance and marker coverage
    species_freqs = species_freq_matrix[species_idx_map[species_name],:]
    
    depth_axis.semilogy(times[marker_coverages>0], marker_coverages[marker_coverages>0],'.-',color='#007ccd',markersize=3)
    species_freq_axis.semilogy(species_times[species_freqs>0], species_freqs[species_freqs>0],'k.-',markersize=3)
    
    if species_freqs[species_freqs>0].min() < 1e-04:
        print species_freqs[species_freqs>0].min()
        species_freq_axis.set_ylim(bottom=1e-04)
        
    species_freq_axis.set_xlim([0,160])   
    

    p = 1 #min([1,1000.0/(len(snp_infos)+1)])

    for gene_idx in xrange(0,gene_copynum_matrix.shape[0]):
        
        num_total_mutations += 1
     
        gene_copynums = gene_copynum_matrix[gene_idx,:]
        gene_name = gene_names[gene_idx]
    
        if not (marker_coverages>=min_coverage).any():
            continue
    
        masked_gene_copynums = gene_copynums[marker_coverages>=min_coverage]
        masked_times = times[marker_coverages>=min_coverage] 
        masked_marker_coverages = marker_coverages[marker_coverages>=min_coverage]
    
        
        
        if color_condition(species_idx, gene_name, masked_times, masked_gene_copynums, masked_marker_coverages):
        
            # One of the colored ones!
            num_colored_mutations+=1
            
            #sys.stderr.write("%s %d %s %s\n" % (gene_name, location, var_type, allele)) 
            
            line, = copynum_axis.semilogy(masked_times, masked_gene_copynums, '-o', alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=COLORED_LINEWIDTH)
            color = pylab.getp(line,'color')
            
        else:  
            # One of the non-colored ones
            #freq_axis.plot(theory_times, interpolation_function(theory_times), '-', alpha=0.5, color='0.7', markersize=3,linewidth=1,zorder=1)
            if random_sample() < p:
                copynum_axis.semilogy(masked_times, masked_gene_copynums, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
     
    sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
    
copynum_axis.set_xlabel('Time, $t$ (days)')

copynum_axis.set_xlim([0,160])   
 
sys.stderr.write("Saving final PNG image...\t")
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")


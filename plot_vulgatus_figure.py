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
import calculate_preexisting_snps
import cluster_utils
################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('settingsfilename', type=str, help="settings-filename")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
settings_filename = args.settingsfilename
################################################################################


       
desired_samples = parse_timecourse_data.morteza_samples
    
# Mininum coverage for frequency estimation vs interpolation 
min_coverage = 10
COLORED_LINEWIDTH=1
    
mpl.rcParams['font.size'] = 7.0
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

# load settings
settings_file = open(settings_filename,"r")
settings_string = "\n".join(settings_file.readlines())
settings_file.close()
exec settings_string    

desired_species = species_names


num_panels = len(desired_species)

fig_width = 7
fig_height = 1.7*(num_panels)

fig, axes = plt.subplots(num_panels,sharex=True,sharey=False,figsize=(fig_width, fig_height))    
axes = [axes]
    
fig_idx = 0
freq_axis = None

output_items = []

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
samples = numpy.array(samples)[sample_idxs]
species_coverage_matrix = species_coverage_matrix[:,sample_idxs]
total_coverage = species_coverage_matrix.sum(axis=0)
species_freq_matrix = numpy.clip(species_coverage_matrix*1.0/total_coverage,0, 2)    

 
for species_idx in xrange(0,len(desired_species)):        

    print species_idx, desired_species[species_idx]

    species_name = desired_species[species_idx]
    
    sys.stderr.write("Processing %s...\n" % species_name)
    
    
    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=desired_samples)
    sys.stderr.write("Done!\n")
    if len(gene_names)==0:
        continue
    
    marker_coverage_times, marker_coverage_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, gene_samples)

    gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
    marker_coverages = marker_coverages[marker_coverage_idxs]
    gene_copynum_matrix = gene_copynum_matrix[:,marker_coverage_idxs]
    
    #preexisting_snps = calculate_preexisting_snps.parse_preexisting_snps(species_name)
    
    #sys.stderr.write("Loaded database of %d snps!\n" % len(preexisting_snps))
    
    times = []
    alt_matrix = []
    depth_matrix = []
    snp_infos = []

    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','2D','3D','4D']),chunk_size=chunk_size,allowed_samples=desired_samples, initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        sample_ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

        # Calculate fixation matrix
        sys.stderr.write("Calculating allele freqs...\n")
        chunk_alts, chunk_depths, chunk_snp_infos = timecourse_utils.calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))    
        sys.stderr.write("Done!\n")
    
        chunk_alts = chunk_alts[:,sample_idxs]
        chunk_depths = chunk_depths[:,sample_idxs]
        
        desired_sites = numpy.logical_not(numpy.logical_or((chunk_alts<=(0.1*chunk_depths)).all(axis=1), (chunk_alts>=(0.9*chunk_depths)).all(axis=1)))*((chunk_depths>0).sum(axis=1)>2)
        
        
            
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
                    
    sys.stderr.write("Done loading SNP freqs!\n")
    
        
    if len(alt_matrix)>0:     
        alt_matrix = numpy.vstack(alt_matrix)
        depth_matrix = numpy.vstack(depth_matrix) 
    else:
        alt_matrix = numpy.array([])
        depth_matrix = numpy.array([])
    

    # set up figure axis
    
    freq_axis = axes[fig_idx]
    fig_idx+=1
    
    title_text = ('%s' % species_name)
    freq_axis.set_title(title_text,loc='right',fontsize=6)
        
    freq_axis.set_ylabel('Allele frequency, $f(t)$')
    freq_axis.spines['top'].set_visible(False)
    freq_axis.spines['right'].set_visible(False)
    freq_axis.get_xaxis().tick_bottom()
    freq_axis.get_yaxis().tick_left()
    
    freq_axis.set_xlim([0,160])   
    freq_axis.set_ylim([0,1.02])
    
    freq_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[0,0],[1.02,1.02],color=parse_timecourse_data.antibiotics_color,linewidth=0)
    
    num_colored_mutations = 0
    num_total_mutations = 0

    if len(alt_matrix)>0:
        # lets plot some SNPs!
        p = min([1,1000.0/(len(snp_infos)+1)])

        cluster_As = []
        cluster_Ds = []

        for mutation_idx in xrange(0,len(snp_infos)):
        
            num_total_mutations += 1
     
            chromosome, location, gene_name, variant_type = snp_infos[mutation_idx]
            alts = alt_matrix[mutation_idx,:]
            depths = depth_matrix[mutation_idx,:]
        
            alts = alts*(depths>=min_coverage)
            depths = depths*(depths>=min_coverage)
        
        
            freqs = alts*1.0/(depths+(depths==0))
        
            masked_times = times[depths>=min_coverage]
            masked_freqs = freqs[depths>=min_coverage]
            masked_depths = depths[depths>=min_coverage]
        
            if len(masked_freqs)<3:
                continue
        
            allele='A'
        
            if masked_freqs[0]>0.5:
                masked_freqs = 1-masked_freqs
                allele='R'
            
            # figure out whether to color it!
            id = (gene_name, chromosome, location, variant_type)
            
            color_condition = False
            #(id in antibiotic_fixation_map[species_name]['snps']) or (id in final_fixation_map[species_name]['snps'])
            
            cluster_As.append(alts)
            cluster_Ds.append(depths)
            
     
        cluster_As = numpy.array(cluster_As)
        cluster_Ds = numpy.array(cluster_Ds)
        sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
        if len(cluster_As)==0:
            continue
        
        cluster_map = cluster_utils.cluster_snps_by_distance(cluster_As, cluster_Ds)
        
        for cluster_id in cluster_map:
            color = None
            avg_fs, total_Ds = cluster_map[cluster_id]['centroid']
            
            good_idxs = (total_Ds>0)
            
            if good_idxs.sum()<1.5:
                continue
            
            cluster_size = len(cluster_map[cluster_id]['snps'])
            fractional_size = cluster_size*1.0/len(cluster_As)
            
            if fractional_size < 0.05:
            #if cluster_size<1000:
                continue
 
            sys.stderr.write("Cluster %d (n=%d SNVs, %g)\n" % (cluster_id, cluster_size, fractional_size)) 
 
            
            masked_times = times[good_idxs]
            masked_avg_fs = avg_fs[good_idxs]
            
            freq_axis.plot(masked_times, masked_avg_fs, 'k-', alpha=1, markersize=3, markeredgecolor='none', zorder=5, linewidth=1)
            line, = freq_axis.plot(masked_times, masked_avg_fs, 'o', alpha=1, markersize=3, zorder=6, linewidth=0.5)
            color = pylab.getp(line,'color')
            
            
            for snp_idx, flip in cluster_map[cluster_id]['snps']:
                good_idxs = cluster_Ds[snp_idx,:]>0
                masked_times = times[good_idxs]
                masked_As = cluster_As[snp_idx,good_idxs]
                masked_Ds = cluster_Ds[snp_idx,good_idxs]
                if flip:
                    masked_As = masked_Ds-masked_As
                masked_freqs = masked_As*1.0/masked_Ds
                    
                freq_axis.plot(masked_times, masked_freqs, '-', color=color, alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=0.5)
           
     
freq_axis.set_ylabel('Allele frequency, $f(t)$')
freq_axis.set_xlabel('Time, $t$ (days)')        
    
filename = parse_midas_data.analysis_directory+('vulgatus_cluster_figure.png')
sys.stderr.write("Saving final PNG image...\t")
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")


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

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fixation_type', type=str, help="antibiotic or final")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--private", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--output-filename', type=str, help='Store colored snps')
args = parser.parse_args()

fixation_type = args.fixation_type
debug = args.debug
chunk_size = args.chunk_size
output_filename = args.output_filename
private_only = args.private
################################################################################


#####
#
# Load precalculated fixation events
#
####
file = open("%s_within_host_output.txt" % fixation_type,"r")
fixation_map = {}
for line in file:
    items = line.split("\t")
    species = items[0]
    change_strs = items[1:]
    snp_changes = []
    gene_changes = []
    
    for change_str in change_strs:
        subitems = change_str[1:-1].split(",")
        if len(subitems)>2:
            # snp change
            snp_changes.append((subitems[0][1:-1], subitems[1].strip()[1:-1], long(subitems[2]), subitems[3].strip()[1:-1]))
        else:
            gene_changes.append(subitems[0][1:-1])
        
    fixation_map[species] = {'snps': set(snp_changes), 'genes': set(gene_changes)}
       
desired_samples = parse_timecourse_data.morteza_samples
desired_species = fixation_map.keys()
if debug:
    desired_species = desired_species[:3]
            
sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.morteza_samples)

species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
species_idx_map = {species[i]: i for i in xrange(0,len(species))}

species_times, species_time_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

species_coverage_matrix = species_coverage_matrix[:,species_time_idxs]
species_freq_matrix = species_coverage_matrix*1.0/(species_coverage_matrix.sum(axis=0))

    
# Mininum coverage for frequency estimation vs interpolation 
min_coverage = 5
COLORED_LINEWIDTH=1
    
mpl.rcParams['font.size'] = 7.0
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

fig_width = 7
fig_height = 1.7*(3*len(desired_species))

fig, axes = plt.subplots(len(desired_species)*3,sharex=True,sharey=False,figsize=(fig_width, fig_height))    

fig_idx = 0
freq_axis = None


output_items = {}
 
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
    
    if private_only:
        preexisting_snps = calculate_preexisting_snps.parse_preexisting_snps(species_name)
    else:
        preexisting_snps = set([])    
    
    sys.stderr.write("Loaded database of %d snps!\n" % len(preexisting_snps))
    
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
    
    species_freq_axis = axes[fig_idx]
    fig_idx += 1
    title_text = '%s abundance' % species_name
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
    
    freq_axis = axes[fig_idx]
    fig_idx+=1
    
    title_text = ('%s SNVs' % species_name)
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

    # Plot species abundance and marker coverage
    species_freqs = species_freq_matrix[species_idx_map[species_name],:]
    
    depth_axis.semilogy(marker_coverage_times[marker_coverages>0], marker_coverages[marker_coverages>0],'.-',color='#007ccd',markersize=3)
    species_freq_axis.semilogy(species_times[species_freqs>0], species_freqs[species_freqs>0],'k.-',markersize=3,zorder=4)
    
    if species_freqs[species_freqs>0].min() < 1e-04:
        print species_freqs[species_freqs>0].min()
        species_freq_axis.set_ylim(bottom=1e-04)
        
    species_freq_axis.set_xlim([0,160])   
    
    if len(alt_matrix)>0:
        # lets plot some SNPs!
        p = min([1,1000.0/(len(snp_infos)+1)])

        for mutation_idx in xrange(0,len(snp_infos)):
        
            num_total_mutations += 1
     
            chromosome, location, gene_name, variant_type = snp_infos[mutation_idx]
            alts = alt_matrix[mutation_idx,:]
            depths = depth_matrix[mutation_idx,:]
        
        
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
            
            color_condition = (id in fixation_map[species_name]['snps'])
            
            if len(preexisting_snps)>0:
                
                #color_condition = color_condition and ((chromosome,location) not in preexisting_snps)
                color_condition = ((chromosome,location) not in preexisting_snps)
            
            if color_condition:
            
            
                output_str = "\t".join([species_name, snp_infos[mutation_idx][0], str(snp_infos[mutation_idx][1]), allele, str(species_idx)])
        
                output_key = (species_idx, chromosome, location)
                output_items[output_key] = output_str
        
                # One of the colored ones!
                num_colored_mutations+=1
            
                #sys.stderr.write("%s %d %s %s\n" % (gene_name, location, var_type, allele)) 
            
                line, = freq_axis.plot(masked_times, masked_freqs, '-o', alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=COLORED_LINEWIDTH)
                color = pylab.getp(line,'color')
            
            else:  
                if random_sample() < p:
                    freq_axis.plot(masked_times, masked_freqs, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
     
        sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
    copynum_axis = axes[fig_idx]
    fig_idx+=1
    
    title_text = '%s genes' % species_name
    copynum_axis.set_title(title_text,loc='right',fontsize=6)
        
    copynum_axis.set_ylabel('Gene copynum')
    copynum_axis.spines['top'].set_visible(False)
    copynum_axis.spines['right'].set_visible(False)
    copynum_axis.get_xaxis().tick_bottom()
    copynum_axis.get_yaxis().tick_left()
    
    copynum_axis.set_xlim([0,160])   
    copynum_axis.set_ylim([3e-02,3])
    
    copynum_axis.fill_between([parse_timecourse_data.antibiotic_start, parse_timecourse_data.antibiotic_end],[1e-01,1e-01],[100,100],color=parse_timecourse_data.antibiotics_color,linewidth=0)
    
    p = 1 
    
    for gene_idx in xrange(0,gene_copynum_matrix.shape[0]):
        
        num_total_mutations += 1
     
        gene_copynums = gene_copynum_matrix[gene_idx,:]
        gene_name = gene_names[gene_idx]
    
        if not (marker_coverages>=min_coverage).any():
            continue
    
        if not ((gene_copynums>=0.3)*(gene_copynums<3)).any():
            continue
    
        masked_gene_copynums = numpy.clip(gene_copynums[marker_coverages>=min_coverage],1e-04,1e02)
        masked_times = marker_coverage_times[marker_coverages>=min_coverage]
        masked_marker_coverages = marker_coverages[marker_coverages>=min_coverage]
    
        if len(masked_gene_copynums)<3:
            continue
    
        masked_gene_copynums = numpy.clip(masked_gene_copynums,5e-02,100)
    
        if gene_name in fixation_map[species_name]['genes']:
        
            output_str = "\t".join([species_name, gene_name])
        
            output_key = (species_idx, gene_name)
            output_items[output_key] = output_str
        
            # One of the colored ones!
            num_colored_mutations+=1
            
            line, = copynum_axis.semilogy(masked_times, masked_gene_copynums, '-o', alpha=0.5, markersize=2, markeredgecolor='none', zorder=4, linewidth=COLORED_LINEWIDTH)
            color = pylab.getp(line,'color')
            
        else:  
            # One of the non-colored ones
            #freq_axis.plot(theory_times, interpolation_function(theory_times), '-', alpha=0.5, color='0.7', markersize=3,linewidth=1,zorder=1)
            if random_sample() < p:
                copynum_axis.semilogy(masked_times, masked_gene_copynums, '-', color='0.7', alpha=0.5, markersize=3,label=gene_name,linewidth=0.25,zorder=1)
     
    sys.stderr.write("Colored=%d, Total=%d\n" % (num_colored_mutations, num_total_mutations))
    
# do this for the last one!     
copynum_axis.set_xlabel('Time, $t$ (days)')
copynum_axis.set_xlim([0,160])   

filename = parse_midas_data.analysis_directory+('%s_fixation_timecourse.png' % fixation_type)
sys.stderr.write("Saving final PNG image...\t")
fig.savefig(filename, bbox_inches='tight', dpi=300, transparent=True)
pylab.close(fig)
sys.stderr.write("Done!\n")


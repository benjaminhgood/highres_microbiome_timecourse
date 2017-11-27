import sys
import parse_timecourse_data
import parse_midas_data
import barcode_utils
import numpy

species_name = "Bacteroides_vulgatus_57955"

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""


 
site_base_map = {} 
snp_groups = {}
desired_sites = []
filename = "Bacteroides_vulgatus_57955_final_snps.txt"
file = open(filename,"r")
file.readline() # header
for line in file:
    items = line.split("\t")
    if items[0].strip() == species_name:
        
        contig = items[1].strip()
        position = long(items[2])
        base = items[3]
        group = long(items[4])
        
        site = (contig, position)
        desired_sites.append(site)
        
        if group not in snp_groups:
            snp_groups[group] = []
        
        snp_groups[group].append(site)

        site_base_map[site] = base

desired_sites_set = set(desired_sites)        
snp_group_sets = {group: set(snp_groups[group]) for group in snp_groups.keys()}

 
num_shared_barcodes_per_site = {}

samples = [parse_timecourse_data.morteza_samples[-2], parse_timecourse_data.highcoverage_end]

#samples = parse_timecourse_data.morteza_samples[-4:]

for sample_name in samples:

    # Load barcodes from disk     
    sys.stderr.write("Loading data...\t")
    allele_barcode_map = barcode_utils.parse_allele_barcodes(species_name, sample_name)
    sys.stderr.write("Done!\n")
    sys.stderr.write("Inverting allele->barcode map...\t")
    barcode_allele_map = barcode_utils.calculate_barcode_allele_map(allele_barcode_map)
    sys.stderr.write("Done!\n")
       
    sys.stderr.write("Calculating shared barcodes per site...\t")        
    barcode_utils.calculate_num_shared_barcodes_per_site(allele_barcode_map, barcode_allele_map, desired_sites=desired_sites, num_shared_barcodes_per_site=num_shared_barcodes_per_site)
    sys.stderr.write("Done!\n")

L00s = []
L11s = []
L10s = []

n00_fourgamete_failed = 0
n11_fourgamete_failed = 0
n10_fourgamete_failed = 0

min_barcodes = 4
num_pairs_with_linkage = 0
num_pairs_fail_fourgamete = 0
for site_1 in desired_sites:
    for site_2 in num_shared_barcodes_per_site[site_1].keys():
        
        if site_1==site_2:
            continue
            
        if site_2 not in desired_sites_set:
            continue
        
        is_linked = (num_shared_barcodes_per_site[site_1][site_2].sum()>=min_barcodes)
        
        if not is_linked:
            continue
        
        base_1 = site_base_map[site_1]
        base_2 = site_base_map[site_2]
        
        linkage_score = barcode_utils.calculate_linkage_score(num_shared_barcodes_per_site[site_1][site_2], base_1, base_2)
        
        if barcode_utils.minimum_gamete_fraction(num_shared_barcodes_per_site[site_1][site_2])>=0.1:
            fails_fourgamete_test = 1
        else:
            fails_fourgamete_test = 0
        
        if site_1 in snp_group_sets[0]:
            if site_2 in snp_group_sets[0]:
                # both in group 0
                L00s.append(linkage_score)
                
                n00_fourgamete_failed += fails_fourgamete_test
            else:
                L10s.append(linkage_score)
                
                n10_fourgamete_failed += fails_fourgamete_test
        
        elif site_1 in snp_group_sets[1]:
            if site_2 in snp_group_sets[1]:
                # both in group 1
                L11s.append(linkage_score)
                n11_fourgamete_failed += fails_fourgamete_test
            else:
                L10s.append(linkage_score)
                n10_fourgamete_failed += fails_fourgamete_test
        else:
            pass
            
L00s = numpy.array(L00s)
L11s = numpy.array(L11s)
L10s = numpy.array(L10s)

print '00:', len(L00s), n00_fourgamete_failed, 'median l =', numpy.median(L00s)
print '11:', len(L11s), n11_fourgamete_failed, 'median l =',numpy.median(L11s)
print '10:', len(L10s), n10_fourgamete_failed, 'median l =',numpy.median(L10s)



# Expect mostly +1s for group 0 -> group 0 comparisons
# Expect mostly +1s for group 1 -> group 1 comparisons
# Expect mostly -1s for in between, except that some will be actually linked!            
        
import stats_utils
import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

pylab.figure(1,figsize=(3.42,2))
pylab.xlabel('Linkage score, $\ell$')
pylab.ylabel('Fraction SNP pairs $\geq \ell$')

pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(L00s)
pylab.step(xs,ns/ns[0],color='b',alpha=0.5,label='within group 1')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(L11s)
pylab.step(xs,ns/ns[0],color='r',alpha=0.5,label='within group 2')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(L10s)
pylab.step(xs,ns/ns[0],color='k',alpha=0.5,label='between groups')

pylab.legend(loc='lower left',frameon=False)

pylab.xlim([-1.05,1.05])
pylab.ylim([0,1.05])
pylab.savefig(parse_midas_data.analysis_directory+'linked_snp_groups.pdf',bbox_inches='tight')

       
       
 
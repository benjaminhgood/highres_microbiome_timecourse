#!/usr/bin/env python 
### This script runs the necessary post-processing of the MIDAS output 
### across all species in serial. 

import os
import sys
import config

if len(sys.argv) > 1:
    argument=sys.argv[1]
else:
    argument = 'all'

# First calculate core genes for each species
os.system('python core_gene_utils.py')

# Make some directories
os.system('mkdir -p %ssnp_prevalences' % config.data_directory)

# Call postprocess_midas_data.py for each species
os.system('python loop_over_species_wrapper.py %s python postprocess_midas_data.py' % argument)

# Calculate substitution rates for the most prevalent species
os.system('python calculate_substitution_rates.py')

# Calculate temporal changes for the most prevalent species
os.system('python calculate_temporal_changes.py')

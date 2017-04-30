###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 

data_directory = os.path.expanduser("~/highres_microbiome_timecourse_data/")
analysis_directory = os.path.expanduser("~/highres_microbiome_timecourse_analysis/")
scripts_directory = os.path.expanduser("~/highres_microbiome_timecourse_scripts/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 20
good_species_min_prevalence = 5
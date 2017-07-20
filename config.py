###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 

data_directory = os.path.expanduser("~/highres_microbiome_timecourse_data/")
analysis_directory = os.path.expanduser("~/highres_microbiome_timecourse_analysis/")
scripts_directory = os.path.expanduser("~/highres_microbiome_timecourse_scripts/")
patric_directory = os.path.expanduser("~/patric_db/")
midas_directory = os.path.expanduser("~/midas_db/")
barcode_directory = os.path.expanduser("~/highres_microbiome_timecourse_barcode_data/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 20
good_species_min_prevalence = 5
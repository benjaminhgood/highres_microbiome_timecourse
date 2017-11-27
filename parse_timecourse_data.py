import parse_midas_data
import config
import numpy
import gene_diversity_utils

focal_patient = 'patient0'

hrv_infection = 36
lyme_infection = 41
antibiotic_start = 53
antibiotic_end = 70

highcoverage_start_1 = "6037"
highcoverage_start_2 = "6038.1"
highcoverage_antibiotic = "1014.2"
highcoverage_end = "6041"

morteza_samples = ['6037',
           '6037.2',
           '6037.3',
           '6038.1',
           '1021',
           '1022',
           '1022.1',
           '1023',
           '1014.2',
           '1025',
           '4021A',
           '4021.1A',
           '4022',
           '4023', #4023.1 absent
           '4024.1', #4025 absent
           '4025.4',
           '4026',
           '4026.2',
           '6041']
           
highcoverage_samples = ['6037',
                        '6038.1',
                        '1021',
                        '1022.1',
                        '1014.2',
                        '4021A',
                        '4023',
                        '4025',
                        '6041']
                        
antibiotics_color = '#bdd7e7'
lyme_color = '#eff3ff'





###############################################################################
#
# Loads metadata for HMP samples 
# Returns map from subject -> map of samples -> set of accession IDs
#
###############################################################################
def parse_subject_sample_map(): 

    subject_sample_map = {}
    
    
    # Then load Kuleshov data 
    file = open(parse_midas_data.scripts_directory+"highres_timecourse_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        
        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}
            
        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()
            
        subject_sample_map[subject_id][sample_id].add(accession_id)
    file.close()
     
    return subject_sample_map 

#####
#
# Loads country metadata for samples
#
#####
def parse_sample_visno_map(): 

    sample_visno_map = {}
    
    file = open(parse_midas_data.scripts_directory+"highres_timecourse_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        visno = float(items[5].strip())
        
        sample_visno_map[sample_id] = visno
    
    file.close()
    return sample_visno_map
    
    
#####
#
# Loads country metadata for samples
#
#####
def parse_sample_time_map(): 

    sample_time_map = {}
    
    file = open(parse_midas_data.scripts_directory+"highres_timecourse_ids.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        subject_id = items[0].strip()
        sample_id = items[1].strip()
        accession_id = items[2].strip()
        country = items[3].strip()
        continent = items[4].strip()
        visno = float(items[5].strip())
        time = float(items[6].strip())
        
        sample_time_map[sample_id] = time
    
    file.close()
    return sample_time_map

    
def calculate_timecourse_idxs(sample_time_map, desired_samples, min_time=1):

    # TODO
    idxs = []
    times = []
    for i in xrange(0,len(desired_samples)):
        sample = desired_samples[i]
        if sample in sample_time_map:
            time = sample_time_map[sample]
            if time>=min_time:
                idxs.append(i)
                times.append(time)
                
    times, idxs = zip(*sorted(zip(times, idxs)))
    times = numpy.array(times)
    idxs = numpy.array(idxs)
    
    return times, idxs
    
###############################################################################
#
# Loads a subset of "core" genes using copynum information in the genes/ folder 
#
###############################################################################   
def load_core_timecourse_genes(desired_species_name, min_copynum=0.3, min_prevalence=0.9, min_marker_coverage=20):

    # Load subject and sample metadata
    subject_sample_map = parse_subject_sample_map()
    sample_time_map = parse_sample_time_map()
    
    desired_samples = set(subject_sample_map[focal_patient].keys())
    
    # Load reference genes
    reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
    
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(desired_species_name)
    
    gene_names = numpy.array(gene_names)
   
    reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])

    sample_idxs = numpy.array([sample_name in desired_samples for sample_name in gene_samples])*(marker_coverages>=min_marker_coverage)
    
    if sample_idxs.sum()>0:   
        prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,sample_idxs], marker_coverages[sample_idxs], min_copynum)
        core_gene_idxs = reference_gene_idxs*(prevalences>=min_prevalence)  
    else:
        sys.stderr.write("Not enough samples for core genome!\n")
        reference_gene_idxs

    return set(gene_names[core_gene_idxs])
                
    
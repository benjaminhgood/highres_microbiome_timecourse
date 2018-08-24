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
highcoverage_hrv = "1021"
highcoverage_lyme = "1022.1"
highcoverage_antibiotic = "1014.2"
highcoverage_postantibiotic = "4023"
highcoverage_end = "6041"

morteza_samples = ['6037',
           '6037.2',
           '6037.3',
           '6038.1', # hc (day 26)
           '1021', # hc, hrv, day 36
           '1022',
           '1022.1', # lyme diagnosed. day 47
           '1023',
           '1014.2', # antibiotic
           '1025',
           '4021A', # hc (5 days post)
           '4021.1A',
           '4022',
           '4023', # hc (11 days post) #4023.1 absent #4025 absent
           '4024.1', 
           '4025.4',
           '4026',
           '4026.2',
           '6041'] # final (73 days post)
           
highcoverage_samples = ['6037', # start 1, day 1
                        '6038.1', # start 2, day 26
                        '1021', # hrv, day 36
                        '1022.1', # lyme, day 47
                        '1014.2', # antibiotic (2 days before finishing), day 68
                        '4021A', # diagnostic 1 (5 days post), day 73
                        '4023', # diagnostic 2 (11 days post), day 81
                        '6041'] # final (84 days post antibiotic), day 154 (day 73) 

alistipes_onderdonkii_gene_gain_samples = ['1014.2', '1025', '4021A', '4021.1A', '4022', '4023']
                        
antibiotics_color = '#bdd7e7'
lyme_color = '#eff3ff'

# total # of reads in the raw fastq files. 
fastq_coverage_map = {} # dividing by 4 because 4 lines per read
fastq_coverage_map['1014.2'] = 1228272236/4
fastq_coverage_map['1021'] = 1276652816/4
fastq_coverage_map['1022'] = 262493376/4
fastq_coverage_map['1022.1'] = 2538752928/4
fastq_coverage_map['1023'] = 217127600/4
fastq_coverage_map['1025'] = 462883992/4
fastq_coverage_map['4021.1A'] = 379311296/4
fastq_coverage_map['4021A'] = 2365861024/4
fastq_coverage_map['4022'] = 328461556/4
fastq_coverage_map['4023'] = 1191171684/4
fastq_coverage_map['4024.1'] = 268695392/4
fastq_coverage_map['4025.4'] = 639037900/4
fastq_coverage_map['4026'] = 582648592/4
fastq_coverage_map['4026.2'] = 323633936/4
fastq_coverage_map['6037'] = 1210346620/4
fastq_coverage_map['6037.2'] = 1208179808/4
fastq_coverage_map['6037.3'] = 130772348/4
fastq_coverage_map['6038.1'] = 1228380944/4
fastq_coverage_map['6041'] = 2354297232/4



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
    
def parse_sample_order_map(): 

    sample_time_map = parse_sample_time_map()
    samples = list(sample_time_map.keys())
    
    sorted_idxs = range(0,len(samples))
    sorted_idxs = list(sorted(sorted_idxs, key = lambda idx: sample_time_map[samples[idx]]))
    
    sample_order_map = {}
    for idx in xrange(0,len(sorted_idxs)):
        sample_order_map[samples[sorted_idxs[idx]]] = 'patient0', idx
        
    return sample_order_map

    
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
    
    if len(idxs)>0:        
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
    
    
def get_initial_antibiotic_final_idxs(samples):
    samples = list(samples)
    if highcoverage_start_2 in samples:
        initial_idx = samples.index(highcoverage_start_2)
    elif highcoverage_start_1 in samples:
        initial_idx = samples.index(highcoverage_start_1)
    elif highcoverage_lyme in samples:
        initial_idx = samples.index(highcoverage_lyme)
    else:
        initial_idx = -1
        
    if highcoverage_antibiotic in samples:
        antibiotic_idx = samples.index(highcoverage_antibiotic)
    elif highcoverage_postantibiotic in samples:
        antibiotic_idx = samples.index(highcoverage_postantibiotic)
    else:
        antibiotic_idx = -1
        
    if highcoverage_end in samples:
        final_idx = samples.index(highcoverage_end)
    else:
        final_idx = -1
        
    return initial_idx, antibiotic_idx, final_idx
        
           
def get_initial_idxs(samples):

    sample_time_map = parse_sample_time_map()
    ts = numpy.array([sample_time_map[sample] for sample in samples])
    
    initial_idxs = (ts<28)
    if initial_idxs.sum() > 0:
        return numpy.nonzero(initial_idxs)[0]
    else:
        return []

def get_antibiotic_idxs(samples):     
    
    
    sample_time_map = parse_sample_time_map()
    ts = numpy.array([sample_time_map[sample] for sample in samples])
    
    initial_idxs = (ts>66)*(ts<83)
    
    if initial_idxs.sum() > 0:
        return numpy.nonzero(initial_idxs)[0]
    else:
        return []

def get_final_idxs(samples):
    
    sample_time_map = parse_sample_time_map()
    ts = numpy.array([sample_time_map[sample] for sample in samples])
    
    initial_idxs = (ts>118)
    
    if initial_idxs.sum() > 0:
        return numpy.nonzero(initial_idxs)[0]
    else:
        return []

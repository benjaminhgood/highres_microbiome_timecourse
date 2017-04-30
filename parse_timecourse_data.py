import parse_midas_data
import numpy

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
        visno = long(items[5].strip())
        
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
        visno = long(items[5].strip())
        time = long(items[6].strip())
        
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
    
    return numpy.array(times), numpy.array(idxs)
                
    
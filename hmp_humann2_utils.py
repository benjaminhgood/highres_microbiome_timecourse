import os
import sys
import gzip
import numpy

hmp_humann2_directory = os.path.expanduser("~/hmp_humann2_data/")

real_hmp_humann2_filename = hmp_humann2_directory+"hmp1-II_genestrata_Stool-mtd-qcd.tsv.gz"

postprocessed_hmp_humann2_filename = hmp_humann2_directory+"hmp_humann2_uniref50.tsv.gz"


def pipe_uniref50_table():
    
    file = gzip.open(real_hmp_humann2_filename,"r")
    header_1_line = file.readline() # header 1
    header_items = numpy.array(header_1_line.strip().split("\t"))
    file.readline() # RANDSID
    visno_line = file.readline() # VISNO
    visno_items = visno_line.strip().split('\t')
    good_idxs = [item.strip()=='1' for item in visno_items]
    good_idxs[0]=True
    good_idxs = numpy.array(good_idxs)
    
    file.readline() # Sample Area (gut)
    file.readline() # Site (stool)
    file.readline() # SNPRNT
    file.readline() # Male/Female
    file.readline() # WMSPhase
    file.readline() # SRS
    
    print "\t".join(header_items[good_idxs])
    
    for line in file:
        items = numpy.array(line.strip().split('\t'))
        
        if '|' in items[0]:
            continue
        
        print "\t".join(items[good_idxs])
        
        

def parse_uniref50_table():
    
    file = gzip.open(postprocessed_hmp_humann2_filename,"r")
    file.readline() # header 1
    
    abundance_matrix = []
    uniref50_names = []
    for line in file:
        items = line.split()
        uniref50_name = items[0]
        if not uniref50_name.startswith('UniRef50'):
            print uniref50_name
        abundances = numpy.array([float(item) for item in items[1:]])
        abundance_matrix.append(abundances)
        uniref50_names.append(uniref50_name)
        
    abundance_matrix = numpy.array(abundance_matrix)
        
    return uniref50_names, abundance_matrix
    
if __name__=='__main__':
    
    #pipe_uniref50_table()
    
    uniref50_names, abundance_matrix = parse_uniref50_table()
    
    
    
    print len(uniref50_names)
    print abundance_matrix.sum(axis=0)
    
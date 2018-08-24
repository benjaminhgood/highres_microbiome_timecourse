import config
import numpy
import sys
import parse_uniref_data

def parse_kegg_module_abundance_matrix(samples, level=3, uniref_abundance_data=None):

    if uniref_abundance_data==None:
        uniref_names, abundance_matrix = parse_gene_family_abundance_matrix(parse_timecourse_data.morteza_samples, include_unmapped=True, include_unknown=True)
    else:
        uniref_names, abundance_matrix = uniref_abundance_data
        
    uniref_kegg_map, kegg_modules = parse_uniref_data.load_uniref_kegg_map(level)
    
    module_idx_map = {kegg_modules[idx] : idx for idx in xrange(0,len(kegg_modules))}
    module_abundance_matrix = numpy.zeros((len(kegg_modules), abundance_matrix.shape[1]))*1.0
    
    for uniref_idx in xrange(0,len(uniref_names)):
        uniref_name = uniref_names[uniref_idx]
        
        if uniref_name not in uniref_kegg_map:
            continue
        
        freqs = abundance_matrix[uniref_idx,:]
        
        target_modules = uniref_kegg_map[uniref_name]
        for module_name in target_modules:
            module_idx = module_idx_map[module_name]
            module_abundance_matrix[module_idx,:] += freqs
            
    # renormalize
    module_abundance_matrix /= module_abundance_matrix.sum(axis=0)
    return kegg_modules, module_abundance_matrix   

def parse_gene_family_abundance_matrix(samples, include_unmapped=True, include_unknown=True):
    
    name_idx_map = {}
    idx_name_map = []
    idx_rpk_map = []
    
    for sample_idx in xrange(0,len(samples)):
        
        sample = samples[sample_idx]
        
        sys.stderr.write("Processing sample %s...\n" % sample)
        
        filename = "%s%s.join_genefamilies.tsv" % (config.humann2_directory, sample)
        file = open(filename,"r")
        if not include_unmapped:
            file.readline() # UNMAPPED READS
        file.readline() # header
    
        for line in file:
    
            if not include_unknown:
                if line.startswith('UniRef90_unknown'):
                    continue
        
            items = line.split()
            name = items[0]
            if "|" in items[0]:
                continue
            
            rpk = float(items[1])
        
            if name not in name_idx_map:
                name_idx_map[name] = len(idx_name_map)
                idx_name_map.append(name)
                idx_rpk_map.append(numpy.zeros(len(samples))*1.0)
            
            name_idx = name_idx_map[name]
            idx_rpk_map[name_idx][sample_idx] = rpk
            
        file.close()
        
    # create giant matrix
    idx_rpk_map = numpy.array(idx_rpk_map)
    # normalize
    idx_rpk_map /= idx_rpk_map.sum(axis=0)
    
    return idx_name_map, idx_rpk_map

def parse_species_specific_gene_family_abundance_matrix(samples, desired_species_names, include_unknown=False):
    
    include_unmapped=False
    
    name_idx_map = {}
    idx_name_map = []
    idx_rpk_map = []
    
    for sample_idx in xrange(0,len(samples)):
        
        sample = samples[sample_idx]
        
        sys.stderr.write("Processing sample %s...\n" % sample)
        
        filename = "%s%s.join_genefamilies.tsv" % (config.humann2_directory, sample)
        file = open(filename,"r")
        if not include_unmapped:
            file.readline() # UNMAPPED READS
        file.readline() # header
    
        for line in file:
    
            if not include_unknown:
                if line.startswith('UniRef90_unknown'):
                    continue
        
            items = line.split()
            name_items = items[0].split("|")
            
            if len(name_items)==1:
                continue # ignore values summed over all species
                
            uniref_name = name_items[0]
            species_name = name_items[1]
            
            if species_name not in desired_species_names:
                continue
                
            rpk = float(items[1])
        
            if uniref_name not in name_idx_map:
                name_idx_map[uniref_name] = len(idx_name_map)
                idx_name_map.append(uniref_name)
                idx_rpk_map.append(numpy.zeros(len(samples))*1.0)
            
            name_idx = name_idx_map[uniref_name]
            idx_rpk_map[name_idx][sample_idx] = rpk
            
        file.close()
        
    # create giant matrix
    idx_rpk_map = numpy.array(idx_rpk_map)
    # normalize
    idx_rpk_map /= idx_rpk_map.sum(axis=0)
    
    return idx_name_map, idx_rpk_map


def parse_gene_family_abundances(sample,include_unmapped=True,include_unknown=True):
    
    filename = "%s%s.join_genefamilies.tsv" % (config.humann2_directory, sample)
    file = open(filename,"r")
    if not include_unmapped:
        file.readline() # UNMAPPED READS
    file.readline() # header
    
    uniref90_names = []
    uniref90_rpks = []
    
    for line in file:
    
        if not include_unknown:
            if line.startswith('UniRef90_unknown'):
                continue
        
        items = line.split()
        name = items[0]
        if "|" in items[0]:
            continue
            
        rpk = float(items[1])
        
        uniref90_names.append(name)
        uniref90_rpks.append(rpk)
     
    file.close()
        
    uniref90_names = numpy.array(uniref90_names)
    uniref90_rpks = numpy.array(uniref90_rpks)
    
    # normalize to copy number
    uniref90_rpks /= uniref90_rpks.sum()

    return uniref90_names, uniref90_rpks

    
def merge_uniref_abundances(names1, abundances1, names2, abundances2):

    name_abundance_map = {}
    for idx in xrange(0,len(names1)):
        
        if names1[idx] not in name_abundance_map:
            name_abundance_map[names1[idx]] = [0,0]
        
        name_abundance_map[names1[idx]][0] = abundances1[idx]
    
    for idx in xrange(0,len(names2)):
        
        if names2[idx] not in name_abundance_map:
            name_abundance_map[names2[idx]] = [0,0]
        
        name_abundance_map[names2[idx]][1] = abundances2[idx]
        
    merged_names = []
    merged_abundance_matrix = []
    for name in name_abundance_map:
        merged_names.append(name)
        merged_abundance_matrix.append(name_abundance_map[name])
      
    merged_names = numpy.array(merged_names)  
    merged_abundance_matrix = numpy.array(merged_abundance_matrix)
    
    return merged_names, merged_abundance_matrix
    
if __name__=='__main__':
    
    import parse_midas_data
    import parse_timecourse_data

    #uniref_names, uniref_rpks = parse_gene_family_abundances(parse_timecourse_data.highcoverage_start_1)
    
    uniref_names, abundance_matrix = parse_gene_family_abundance_matrix(parse_timecourse_data.morteza_samples, include_unmapped=True, include_unknown=True)
    print 'UniRef90_Q08425' in uniref_names
    print 'UniRef90_A0A1Q3L398' in uniref_names
    uniref_names = numpy.array(uniref_names)
    tetQ_idxs = (uniref_names == 'UniRef90_Q08425')
    if tetQ_idxs.sum()>0:
        tetQ_idx = numpy.nonzero(tetQ_idxs)[0]
    else:
        tetQ_idx = -1
        
    print tetQ_idx
    
    print uniref_names[-1], uniref_names[10]
    print len(uniref_names)
    print abundance_matrix[tetQ_idx,:]
    
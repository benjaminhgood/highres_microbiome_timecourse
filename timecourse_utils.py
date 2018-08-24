import numpy
from scipy.interpolate import interp1d

def calculate_read_count_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']),allowed_genes=set([])):

    if len(allowed_genes) == 0:
        allowed_genes = set(passed_sites_map.keys())
     
    alt_matrix = []
    depth_matrix = []
    snp_infos = []
     
    for gene_name in sorted(allowed_genes):
    
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue
        
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue
                
            
            depths = allele_counts.sum(axis=2)
            alts = allele_counts[:,:,0]
        
            for site_idx in xrange(0,len(allele_counts_map[gene_name][variant_type]['locations'])):
                chromosome, location = allele_counts_map[gene_name][variant_type]['locations'][site_idx]
                snp_infos.append((chromosome, location, gene_name, variant_type))
                
            alt_matrix.append(alts)
            depth_matrix.append(depths) 
    

    alt_matrix = numpy.vstack(alt_matrix)
    depth_matrix = numpy.vstack(depth_matrix)
    
    return alt_matrix, depth_matrix, snp_infos


###########
#
# Creates a continuous interpolation function for frequency trajectory,
# so that it can be evaluated anywhere on time interval
#
###########    
def create_interpolation_function(times,freqs,kind='linear'):
    # can create it for anything!

    interpolating_function = interp1d(times, freqs, kind=kind,bounds_error=True)
    
    return interpolating_function

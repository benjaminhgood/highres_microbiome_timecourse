import numpy
import diversity_utils
import calculate_substitution_rates
import parse_timecourse_data
import parse_midas_data
import sys
import calculate_temporal_changes

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.morteza_samples)
all_samples = numpy.array(parse_timecourse_data.morteza_samples)[sample_idxs]

desired_samples = all_samples


def calculate_combined_within_species_changes(temporal_change_map, sample_1, sample_2):
    
    snp_perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2)
    
    gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2)
    
    if snp_perr>-0.5 and gene_perr>-0.5:
        
        total_perr = snp_perr + gene_perr
        total_snp_changes = mutations+reversions
        total_gene_changes = gains+losses
        total_changes = []
        for snp_change in total_snp_changes:
            
            id = (snp_change[0], snp_change[1], snp_change[2], snp_change[3]) 
                # gene name,       contig         position       variant_type
            properties = (snp_change[4], snp_change[5], snp_change[6], snp_change[7])
            
            total_changes.append(id)
        
        for  gene_change in total_gene_changes:
            id = (gene_change[0],)
            properties = (gene_change[1], gene_change[2], gene_change[3], gene_change[4])
            total_changes.append(id)
            
        return total_perr, set(total_changes)
        
    else:
        return -1, set([])

good_species_list = parse_midas_data.parse_good_species_list()
for species_name in good_species_list:
    # Only plot samples above a certain depth threshold that are confidently phaseable.

    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")    

    initial_idxs = parse_timecourse_data.get_initial_idxs(desired_samples)
    antibiotic_idxs = parse_timecourse_data.get_antibiotic_idxs(desired_samples)
    final_idxs = parse_timecourse_data.get_final_idxs(desired_samples)
    
    if len(initial_idxs)<1:
        continue
  
    antibiotic_perr=-1    
    # Calculate changes between any intial timepoint and any antibiotic 
    antibiotic_changes = set([])
    for initial_idx in initial_idxs:
        sample_1 = desired_samples[initial_idx]
        for antibiotic_idx in antibiotic_idxs:
            sample_2 = desired_samples[antibiotic_idx]
            chunk_perr, chunk_changes = calculate_combined_within_species_changes(temporal_change_map, sample_1, sample_2)
            if chunk_perr>-0.01:
                antibiotic_changes.update(chunk_changes)
                antibiotic_perr = max([antibiotic_perr,chunk_perr])    
    # Calculate changes between any intial timepoint and any final
     
    final_changes = set([])
    final_perr = -1
    for initial_idx in initial_idxs:
        sample_1 = desired_samples[initial_idx]
        for final_idx in final_idxs:
            sample_2 = desired_samples[final_idx]
            
            chunk_perr, chunk_changes = calculate_combined_within_species_changes(temporal_change_map, sample_1, sample_2)
            if chunk_perr>-0.01:
                final_changes.update(chunk_changes)
                final_perr = max([final_perr,chunk_perr])
    
    #if len(final_changes)>0 and final_perr>-0.5:
    #    print "\t".join([species_name] + [str(change) for change in final_changes])
    #    pass

    if len(antibiotic_changes)>0 and antibiotic_perr>-0.5:
        print "\t".join([species_name] + [str(change) for change in antibiotic_changes])
        pass
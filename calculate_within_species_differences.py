import numpy
import pylab
import diversity_utils
import calculate_substitution_rates
import parse_timecourse_data
import parse_midas_data
import sys
import calculate_temporal_changes

sample_time_map = parse_timecourse_data.parse_sample_time_map()
ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, parse_timecourse_data.morteza_samples)
all_samples = numpy.array(parse_timecourse_data.morteza_samples)[sample_idxs]

desired_samples = numpy.array([parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_end])



desired_samples_with_antibiotic = numpy.array([parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_end])

pylab.figure(1,figsize=(4,3))
fig = pylab.gcf()

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
            
            total_changes.append((id, [properties]))
        
        for  gene_change in total_gene_changes:
            id = (gene_change[0],)
            properties = (gene_change[1], gene_change[2], gene_change[3], gene_change[4])
            total_changes.append((id, [properties]))
            
        return total_perr, total_changes
        
    else:
        return -1, []

def merge_change_lists(perr_1, changes_1, perr_2, changes_2):
    
    if perr_1 < -0.5:
        return perr_2, changes_2
        
    if perr_2 < -0.5:
        return perr_1, changes_1
    
    # first find changes that are in both lists
    total_change_map = {}
    num_hits_map = {}
    for id, property_list in (changes_1+changes_2):
        
        if id not in total_change_map:
            total_change_map[id] = (id, [])
            num_hits_map[id] = 0
        
        total_change_map[id][1].extend(property_list)
        num_hits_map[id] += 1
        
    total_changes = []
    for id in total_change_map:
        if num_hits_map[id] == 2:
            total_changes.append(total_change_map[id])
            
    return max([perr_1, perr_2]), total_changes
        
good_species_list = parse_midas_data.parse_good_species_list()
for species_name in good_species_list:
    
    # Only plot samples above a certain depth threshold that are confidently phaseable.
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name)
    
    #sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    #sys.stderr.write("Done!\n")
    
    # first see if there are changes between begnning and beginning 2
    # (sort of a control)
    start_perr, start_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_start_2)
    
    # first see if there are changes between begnning_1 and antibiotic
    antibiotic_1_perr, antibiotic_1_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_antibiotic)
    
    # first see if there are changes between begnning_2 and antibiotic
    antibiotic_2_perr, antibiotic_2_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_antibiotic)
    
    # first see if there are changes between begnning_2 and antibiotic
    antibiotic_3_perr, antibiotic_3_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_lyme, parse_timecourse_data.highcoverage_antibiotic)
    
    
    if antibiotic_1_perr < -0.5 and antibiotic_2_perr < -0.5 and antibiotic_3_perr<-0.5:
        # species is not present during antibiotic, look for changes
        # at next non-antibiotic timepoint (e.g. recolonization)
        
        antibiotic_1_perr, antibiotic_1_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_postantibiotic)
        antibiotic_2_perr, antibiotic_2_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_postantibiotic)
        antibiotic_3_perr, antibiotic_2_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_lyme, parse_timecourse_data.highcoverage_postantibiotic)
        
        # you can't relax from antibiotic to post-antibiotic if you
        # aren't present during antibiotic!
        relaxation_perr = -1
        relaxation_changes = []
    else:
        # present during antibiotic and after, so 
        # calculate changes during "relaxation phase"
        relaxation_perr, relaxation_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_antibiotic, parse_timecourse_data.highcoverage_postantibiotic)
        
    antibiotic_perr, antibiotic_changes = merge_change_lists(antibiotic_1_perr, antibiotic_1_changes, antibiotic_2_perr, antibiotic_2_changes)  
    antibiotic_perr, antibiotic_changes = merge_change_lists(antibiotic_perr, antibiotic_changes, antibiotic_3_perr, antibiotic_3_changes)
      
    post_perr, post_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_postantibiotic, parse_timecourse_data.highcoverage_end)
    
    # Between S->lyme
    lyme_1_perr, lyme_1_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_lyme)
    lyme_2_perr, lyme_2_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_lyme)
    
    # Between S1->end
    whole_1_perr, whole_1_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_1, parse_timecourse_data.highcoverage_end)
    
    # S2->end
    whole_2_perr, whole_2_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_start_2, parse_timecourse_data.highcoverage_end)
    
    # lyme->end
    whole_3_perr, whole_3_changes = calculate_combined_within_species_changes(temporal_change_map, parse_timecourse_data.highcoverage_lyme, parse_timecourse_data.highcoverage_end)
    
    whole_perr, whole_changes = merge_change_lists(whole_1_perr, whole_1_changes, whole_2_perr, whole_2_changes)  
    whole_perr, whole_changes = merge_change_lists(whole_perr, whole_changes, whole_3_perr, whole_3_changes)  
    
    whole_and_antibiotic_perr, whole_and_antibiotic_changes = merge_change_lists(antibiotic_perr, antibiotic_changes, whole_perr, whole_changes)
    
    if antibiotic_perr > -0.5 or whole_perr > -0.5:
        if False and (len(antibiotic_changes)>0 or len(whole_changes)>0):
            print species_name
            #print "S1->S2", len(start_changes), start_perr
            print "S->antibiotic", len(antibiotic_changes), "(%d, %d, %d)" % (len(antibiotic_1_changes), len(antibiotic_2_changes), len(antibiotic_3_changes)), antibiotic_1_perr, antibiotic_2_perr
            #print "antibiotic->relaxation", len(relaxation_changes), relaxation_perr
            #print "relaxation->end", len(post_changes), post_perr
            print "S->end", len(whole_changes), "(%d, %d, %d)" % (len(whole_1_changes), len(whole_2_changes), len(whole_3_changes)), whole_1_perr, whole_2_perr
            print "S->antibiotic->end", len(whole_and_antibiotic_changes)
            #sys.stderr.write("%s\n" % species_name)
            #sys.stderr.write(" ".join(["antibiotic & end", str(len(whole_and_antibiotic_changes)), str((len(whole_and_antibiotic_changes)+(len(antibiotic_changes)==0))*1.0/(len(antibiotic_changes)-(len(antibiotic_changes)==0)))]))
            #sys.stderr.write("\n")
        
        if len(antibiotic_changes)>0:
            print "\t".join([species_name] + [str(change[0]) for change in whole_changes])
            pass

    continue

    sys.stderr.write("Calculating Fst values...\n")           
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    
    desired_samples = numpy.array(parse_timecourse_data.morteza_samples)

    samples, pi_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_pi', allowed_samples=desired_samples)
    
    ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)

    if parse_timecourse_data.highcoverage_start_1 in samples:
        initial_idx = samples.index(parse_timecourse_data.highcoverage_start_1)
    elif parse_timecourse_data.highcoverage_start_2 in samples:
        initial_idx = samples.index(parse_timecourse_data.highcoverage_start_2)
    elif parse_timecourse_data.highcoverage_lyme in samples:
        initial_idx = samples.index(parse_timecourse_data.highcoverage_lyme)
    else:
        continue
    
    if parse_timecourse_data.highcoverage_antibiotic in samples:
        antibiotic_idx = samples.index(parse_timecourse_data.highcoverage_antibiotic)
    elif parse_timecourse_data.highcoverage_postantibiotic in samples:
        antibiotic_idx = samples.index(parse_timecourse_data.highcoverage_postantibiotic)
    else:
        continue
        
    if parse_timecourse_data.highcoverage_end in samples:
        final_idx = samples.index(parse_timecourse_data.highcoverage_end)
    else:
        continue
    
    # Make sure we have initial, antibiotic, and end    
                 
    normalized_pi_matrix = pi_matrix*1.0/opportunity_matrix
    self_pis = numpy.diag(normalized_pi_matrix)
    avg_pis = (self_pis[:,None]+self_pis[None,:])/2.0
        
    #fst_matrix = 1.0 -  avg_pis/(normalized_pi_matrix+(normalized_pi_matrix==0))
    #fst_matrix = numpy.clip(fst_matrix,1e-05,1e08)
    fst_matrix = normalized_pi_matrix
    
    fst_antibiotic = fst_matrix[initial_idx, antibiotic_idx]
    fst_final = fst_matrix[initial_idx, final_idx]
    dfst = fst_matrix[antibiotic_idx, final_idx]
    
    if fst_antibiotic > 0.1 or fst_final> 0.1:
        print species_name, fst_antibiotic, fst_final, dfst
    
    fsts = fst_matrix[initial_idx,:][sample_idxs]
    if fst_antibiotic > fst_final:
        fig_idx = 2
    else:
        fig_idx = 3
    
    pylab.figure(fig_idx)
    pylab.plot(ts, fsts,'.-')
    
    continue
    
    fst0 = numpy.fabs(fst_matrix[initial_idx,1])
    fstf = numpy.fabs(fst_matrix[0,2])
    pylab.figure(1)
    pylab.loglog([fst0],[fstf],'.')
    pylab.ylim([1e-05,1e-01])
    pylab.xlim([1e-05,1e-01])
    pylab.loglog([1e-05,1e03],[1e-05,1e03],'k-',color='0.7')
    
    
    
    
    if fstf > 2*fst0 and fstf > 1e-04:
        
        # load full thing:
        samples, pi_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_pi', allowed_samples=all_samples)
    
        ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
        
        normalized_pi_matrix = pi_matrix*1.0/opportunity_matrix
        self_pis = numpy.diag(normalized_pi_matrix)
        
        fst_matrix = normalized_pi_matrix - (self_pis[:,None]+self_pis[None,:])/2.0
        fst_matrix = numpy.fabs(fst_matrix)
        fst_matrix = numpy.clip(fst_matrix,1e-05,1e08)
    
        
    
    
    samples, pi_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=desired_samples_with_antibiotic)
    ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
    
    if len(samples)<4:
        continue
        
    normalized_pi_matrix = pi_matrix*1.0/opportunity_matrix
    self_pis = numpy.diag(normalized_pi_matrix)
        
    fst_matrix = normalized_pi_matrix # - (self_pis[:,None]+self_pis[None,:])/2.0
    fst_matrix = numpy.fabs(fst_matrix)
    fst_matrix = numpy.clip(fst_matrix,1e-05,1e08)
    fst0 = numpy.fabs(fst_matrix[0,1])
    fstf = numpy.fabs(fst_matrix[0,2])
    pylab.figure(3)
    pylab.loglog([fst0],[fstf],'.')
    pylab.ylim([1e-05,1e-01])
    pylab.xlim([1e-05,1e-01])
    pylab.loglog([1e-05,1e03],[1e-05,1e03],'k-',color='0.7')
    
    if fstf > 2*fst0 and fstf > 1e-04:
        
        # load full thing:
        samples, pi_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all_pi', allowed_samples=all_samples)
    
        ts, sample_idxs = parse_timecourse_data.calculate_timecourse_idxs(sample_time_map, samples)
        
        normalized_pi_matrix = pi_matrix*1.0/opportunity_matrix
        self_pis = numpy.diag(normalized_pi_matrix)
        
        fst_matrix = normalized_pi_matrix - (self_pis[:,None]+self_pis[None,:])/2.0
        fst_matrix = numpy.fabs(fst_matrix)
        fst_matrix = numpy.clip(fst_matrix,1e-05,1e08)
    
        fsts = numpy.fabs(fst_matrix[0,:])
    
        pylab.figure(4)
        print species_name
        pylab.plot(ts, fsts,'.-')
    
pylab.figure(1)
fig = pylab.gcf()  
pylab.figure(2)
fig = pylab.gcf() 
fig.savefig('%s/within_species_differences_antibiotic.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
pylab.figure(3)
fig = pylab.gcf()
fig.savefig('%s/within_species_differences_end.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
pylab.figure(4)
fig = pylab.gcf()
    
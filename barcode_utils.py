import config
import gzip 
import numpy
import os.path

# Returns list of barcodes and what alleles they were found on.

base_idx_map = {('R','R'):0, ('A','A'):1, ('R','A'):2, ('A','R'):3}


def calculate_linkage_score(gamete_vector, base_1='A', base_2='A'):
    
    n = gamete_vector.sum()
    
    linkage_score = (gamete_vector[0]+gamete_vector[1])*1.0/n - (gamete_vector[2]+gamete_vector[3])*1.0/n        

    if base_1!=base_2:
        linkage_score *= -1
        
    return linkage_score

def minimum_gamete_fraction(gamete_vector):
    
    n = gamete_vector.sum()
    return gamete_vector.min()*1.0/n


def barcodes_exist(species_name, sample_name):
    barcode_filename = "%s%s/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
        
    if os.path.isfile(barcode_filename):
        return True
    else:
        return False
        

def parse_allele_barcode_tuples(species_name, sample_name):
    
    barcode_filename = "%s%s/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
    
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    barcode_file.readline()
    
    allele_barcode_map = {}
    
    for line in barcode_file:
        line = line.strip()
        items = line.split("\t")
        allele = items[0].strip()
        allele_barcode_map[allele] = []
        
        if len(items)>1:
        
            barcode_items = items[1].split(",")
            for barcode_item in barcode_items:
                barcode_subitems = barcode_item.split(":")
                barcode_id = long(barcode_subitems[0])
                barcode_weight = long(barcode_subitems[1])
                
                allele_barcode_map[allele].append((barcode_id, barcode_weight))
                
    return allele_barcode_map

    
def parse_allele_barcodes(species_name, sample_name):
    
    barcode_filename = "%s%s/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
    
    barcode_file = gzip.GzipFile(barcode_filename,"r")
    barcode_file.readline()
    
    allele_barcode_map = {}
    
    for line in barcode_file:
        line = line.strip()
        items = line.split("\t")
        allele = items[0].strip()
        allele_barcode_map[allele] = []
        
        if len(items)>1:
        
            barcode_items = items[1].split(",")
            for barcode_item in barcode_items:
                barcode_subitems = barcode_item.split(":")
                barcode_id = long(barcode_subitems[0])
                
                allele_barcode_map[allele].append(barcode_id)
                
    return allele_barcode_map
    
def calculate_barcode_allele_map(allele_barcode_map):

    barcode_allele_map = {}
    for allele in allele_barcode_map.keys():
        
        for barcode in allele_barcode_map[allele]:
            
            if barcode not in barcode_allele_map:
                barcode_allele_map[barcode] = []
            
            barcode_allele_map[barcode].append(allele)
    
    return barcode_allele_map

def calculate_num_shared_barcodes(allele_barcode_map, barcode_allele_map, desired_alleles):
    
    num_shared_barcodes = {allele: {} for allele in desired_alleles}
    for allele in desired_alleles:
        
        if allele not in allele_barcode_map:
            continue
            
        barcodes = allele_barcode_map[allele]
        
        if len(barcodes) == 0:
            continue
            
        for barcode in barcodes:
            for other_allele in barcode_allele_map[barcode]:
                if other_allele not in num_shared_barcodes[allele]:
                    num_shared_barcodes[allele][other_allele] = 0
                    
                num_shared_barcodes[allele][other_allele] += 1
                
                
    return num_shared_barcodes
    
def calculate_num_shared_barcodes_per_site(allele_barcode_map, barcode_allele_map, desired_sites=set([]), num_shared_barcodes_per_site = {}):
    
    desired_sites_set = set(desired_sites)
    
    allele_items_map = {}
    
    for allele_1 in allele_barcode_map.keys():
        
        if not ((allele_1[-1]=='A') or (allele_1[-1]=='R')):
            # not a SNP
            continue    
        
        if not allele_1 in allele_items_map:
            
            allele_items_1 = allele_1.split("|")
            contig_1 = allele_items_1[0]
            position_1 = long(allele_items_1[1]) 
            base_1 = allele_items_1[2].strip()
            site_1 = (contig_1, position_1)
            allele_items_map[allele_1] = (site_1, base_1)
            
        site_1, base_1 = allele_items_map[allele_1]
        
        if site_1 not in desired_sites_set:
            continue
        
        if site_1 not in num_shared_barcodes_per_site:
            num_shared_barcodes_per_site[site_1] = {}
    
        barcodes = allele_barcode_map[allele_1]
        
        if len(barcodes)==0:
            # only linked to itself!
            continue
            
        for barcode in barcodes:
        
            if len(barcode_allele_map[barcode]) < 2:
                # only linked to itself!
                continue
                
            for allele_2 in barcode_allele_map[barcode]:
                
                if not ((allele_2[-1]=='A') or (allele_2[-1]=='R')):
                    # not a SNP
                    continue
            
                if not allele_2 in allele_items_map:
        
                    allele_items_2 = allele_2.split("|")
                    contig_2 = allele_items_2[0]
                    position_2 = long(allele_items_2[1]) 
                    base_2 = allele_items_2[2].strip()
                    site_2 = (contig_2, position_2)
                    allele_items_map[allele_2] = (site_2, base_2)
            
                site_2, base_2 = allele_items_map[allele_2]
         
                if site_2 not in num_shared_barcodes_per_site[site_1]:
                    num_shared_barcodes_per_site[site_1][site_2] = numpy.zeros(len(base_idx_map))
        
                num_shared_barcodes_per_site[site_1][site_2][base_idx_map[(base_1, base_2)]] += 1
                
    return num_shared_barcodes_per_site

    
def calculate_linked_set(num_shared_barcodes, min_shared_barcodes=1):
    
    linked_set = set([])
    linked_set.update(num_shared_barcodes) 
    
    for allele in num_shared_barcodes:
        for other_allele in num_shared_barcodes[allele].keys():
            if num_shared_barcodes[allele][other_allele] >= min_shared_barcodes:
                linked_set.add(other_allele)
                
    return linked_set

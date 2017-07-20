import config
import gzip 

# Returns list of barcodes and what alleles they were found on.

def parse_allele_barcodes(species_name, sample_name):
    
    barcode_filename = "%s/%s/barcodes/output/%s.barcodes.gz" % (config.barcode_directory, sample_name, species_name)
    
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
    
def calculate_linked_set(num_shared_barcodes, min_shared_barcodes=1):
    
    linked_set = set([])
    linked_set.update(num_shared_barcodes) 
    
    for allele in num_shared_barcodes:
        for other_allele in num_shared_barcodes[allele].keys():
            if num_shared_barcodes[allele][other_allele] >= min_shared_barcodes:
                linked_set.add(other_allele)
                
    return linked_set

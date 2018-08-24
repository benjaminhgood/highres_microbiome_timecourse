import parse_midas_data

def get_genus_name(species_name):
    return species_name.split("_")[0]
    
def get_genus_family_map():

    genus_family_map = {}
    
    file = open("%sgenome_taxonomy.txt" % parse_midas_data.midas_directory,"r")
    file.readline() # header
    for line in file:
        items =  line.split("\t")
        family = items[7].strip()
        genus = items[8].strip()
        genus_family_map[genus] = family
        
    return genus_family_map
    
def get_family_name(species_name, genus_family_map={}):
    
    if len(genus_family_map)==0:
        genus_family_map = get_genus_family_map()
        
    genus_name = get_genus_name(species_name)
    if genus_name in genus_family_map:
        return genus_family_map[genus_name]
    else:
        return genus_name
        

def get_genus_phylum_map():

    genus_phylum_map = {}
    
    file = open("%sgenome_taxonomy.txt" % parse_midas_data.midas_directory,"r")
    file.readline() # header
    for line in file:
        items =  line.split("\t")
        phylum = items[4].strip()
        genus = items[8].strip()
        genus_phylum_map[genus] = phylum
        
    return genus_phylum_map
    
def get_phylum_name(species_name, genus_phylum_map={}):
    
    if len(genus_phylum_map)==0:
        genus_phylum_map = get_genus_phylum_map()
        
    genus_name = get_genus_name(species_name)
    if genus_name in genus_phylum_map:
        return genus_phylum_map[genus_name]
    else:
        return genus_name
    
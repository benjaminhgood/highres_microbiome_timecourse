import config 
import json
import gzip

def load_uniref_kegg_map(level=3):
    # returns a map from uniref90 family to KOXXXXX gene
    module_ko_map, ko_module_map = load_kegg_map(level)

    uniref_kegg_map = {}
    kegg_modules = set()
    
    filename = filename = "%shumann_maps/map_ko_uniref90.txt.gz" % config.uniref_directory
    file = gzip.GzipFile(filename,"r")
    for line in file:
        items = line.split()
        ko_name = items[0]
        
        if ko_name not in ko_module_map:
            continue
            
        module_names = ko_module_map[ko_name]
        
        for uniref_name in items[1:]:
            if uniref_name not in uniref_kegg_map:
                uniref_kegg_map[uniref_name] = []
            uniref_kegg_map[uniref_name].extend(module_names)
            kegg_modules.update(module_names)
    file.close()
    
    return uniref_kegg_map, list(sorted(kegg_modules))
    

def load_kegg_map(level=4):

    # possible levels:
    # level 1: pathway module, structural complex, etc.
    # level 2: Energy metabolism, Carbohydrate and Lipid metabolism, etc...
    # level 3: Central carbohydrate metabolism, Other carbohydrate metabolism, etc..
    # level 4: MXXXXX module, contains KOXXXXX genes
    
    filename = "%skegg_maps/ko00002.json" % config.uniref_directory
    file = open(filename,"r")
    data = json.load(file)
    file.close()
    data = data['children']


    #data = data['children'][0]

    # data object has 4 layers of hierarchy
    # layer 1: pathway module, structural complex, etc.
    # layer 2: Energy metabolism, Carbohydrate and Lipid metabolism, etc...
    # layer 3: Central carbohydrate metabolism, Other carbohydrate metabolism, etc..
    # layer 4: MXXXXX module, contains KOXXXXX genes

    kegg_map = {}
    
    
    for layer1_idx in xrange(0,len(data)):
        
        layer1_name = data[layer1_idx]['name']
        layer1_data = data[layer1_idx]['children']
        
        #if layer1_name != "Pathway module":
        #    continue
        
        #print layer1_name
        
        for layer2_idx in xrange(0,len(layer1_data)):
            
            layer2_name = layer1_data[layer2_idx]['name']
            layer2_data = layer1_data[layer2_idx]['children']
            
            #print layer2_name
            
            for layer3_idx in xrange(0,len(layer2_data)):
            
                layer3_name = layer2_data[layer3_idx]['name']
                layer3_data = layer2_data[layer3_idx]['children']
                
                #print layer3_name
                
                for layer4_idx in xrange(0,len(layer3_data)):
            
                    layer4_name = layer3_data[layer4_idx]['name'].split()[0]
                    layer4_data = layer3_data[layer4_idx]['children']
                    
                    #print layer4_name
                    
                    if level==1:
                        module_name = layer1_name
                    elif level==2:
                        module_name = layer2_name
                    elif level==3:
                        module_name = layer3_name
                    elif level==4:
                        module_name = layer4_name
                    
                    # promote drug resistance to top:
                    if module_name.startswith('Drug'):
                        module_name = "Drug resistance or efflux"
                    
                    if module_name not in kegg_map:
                        kegg_map[module_name] = set()
                    
                    # Layer 4 is the MXXXXXX module, so get KO names
                    for object in layer4_data:
                        ko_name = object['name'].split()[0]
                        kegg_map[module_name].add(ko_name)
    
    # Invert map, from KOXXXXX label to coarse-grained label
    ko_module_map = {}
    
    for module_name in kegg_map.keys():
        for ko_name in kegg_map[module_name]:
            
            if ko_name not in ko_module_map:
                ko_module_map[ko_name] = []
                
            ko_module_map[ko_name].append(module_name)
                       
    return kegg_map, ko_module_map
            
    
if __name__=='__main__':
    module_ko_map, ko_module_map = load_kegg_map(level=3)
    
    #print len(module_ko_map), "top-level modules"
    #print module_ko_map.keys()[0], module_ko_map[module_ko_map.keys()[0]]
    
    num_doubles = 0
    for ko_name in ko_module_map:
        if len(ko_module_map[ko_name])>1:
            num_doubles+=1
            #print ko_name, ko_module_map[ko_name]
            
    #print  len(ko_module_map), "KO groups, ", num_doubles, "assigned to multiple modules"
    
    uniref_kegg_map, kegg_modules = load_uniref_kegg_map()
    num_doubles = 0
    for uniref_name in uniref_kegg_map:
        if len(uniref_kegg_map[uniref_name])>1:
            num_doubles+=1
            #print ko_name, ko_module_map[ko_name]
    
    print  len(uniref_kegg_map), "Uniref90 families, ", len(kegg_modules), "modules, ",  num_doubles, " families assigned to multiple modules"
    

import config
import gzip
import sys

def parse_uniref_annotations(filename=""):
    uniref_id_map = {}
    annotation_id_map = {}
    id_annotation_map = []
    
    num_lines = 1
    
    if len(filename)==0:
        filename = "%s/idmapping_go.tab.gz" % config.uniref_directory
    file = gzip.GzipFile(filename,"r")
    for line in file:
        if num_lines%1000000==0:
            sys.stderr.write("%dM\n" % (num_lines/1000000))
        num_lines+=1
        
        items = line.split("\t")
        uniref_id = items[0].strip()
        go_strs = items[1].split("; ")
        go_ids = []
        for go_str in go_strs:
            if go_str not in annotation_id_map:
                id = len(id_annotation_map) 
                annotation_id_map[go_str] = id
                id_annotation_map.append(go_str)
            go_ids.append(annotation_id_map[go_str])
            
        uniref_id_map[uniref_id] = go_ids
        
    file.close()
    
    return uniref_id_map, id_annotation_map    

if __name__=='__main__':
    
    # test to see about go annotations
    uniref_annotation_id_map, id_annotation_map = parse_uniref_annotations()
    filename = "%s/idmapping_go.tab.gz" % config.uniref_db
    file = gzip.GzipFile(filename,"r")
    for line in file:
        items = line.split("\t")
        uniref_id = items[0].strip()
        uniref90_id = items[8]
        uniref50_id = items[9]
        if uniref90_id == "":
            continue
            
        if uniref50_id == "":
            continue
        
        uniref90_id = uniref90_id[9:]
        uniref50_id = uniref50_id[9:]
        
        # compare thigns
        print uniref_annotation_id_map[uniref_id], uniref_annotation_id_map[uniref90_id] #, uniref_annotation_id_map[uniref50_id]
        
    file.close()
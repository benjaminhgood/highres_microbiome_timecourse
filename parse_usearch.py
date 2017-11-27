import sys
import os
import os.path
import gzip
import shutil

################################################################################


# read in the usearch files and generate files for midas
'''
##############################################################
# i. create the centroid_functions.txt.gz
# this is a list of gene names in the centroids.ffn.gz file
##############################################################

inFN='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_95.fa'
inFile=open(inFN,'r')

outFN='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/files_for_midas/pan_genomes/centroid_functions.txt.gz'

outFile=gzip.open(outFN, 'wb')

outFile.write('gene_id\tfunction_id\tontology\n')

for line in inFile:
    line=line.strip()
    if line[0]=='>':
        header=line[1:]
        outFile.write(header + '\tNA\tNA\n')
            
##############################################################
# ii. copy over centroids.ffn.gz
# replicates what is in the fasta file output of usearch
############################################################## 
inFN='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_95.fa'
outFN='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/files_for_midas/pan_genomes/centroids.ffn.gz'


inFile=open(inFN, 'r') 
outFile=gzip.open(outFN, 'wb')

shutil.copyfileobj(inFile, outFile)
'''
##############################################################
# iii. gene_info.txt.gz
# List of all genes with their centroid mapping. Could be comprehensive but may not be necessary (i.e. I will just list the centroids only)
# see http://drive5.com/usearch/manual/opt_uc.html for parsing output 
############################################################## 

inFN_95='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_95.uc'
inFile_95=open(inFN_95, 'r')
inFN_99='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_99.uc'
inFile_99=open(inFN_99, 'r')

outFN='/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/files_for_midas/pan_genomes/gene_info.txt.gz'
outFile=gzip.open(outFN, 'wb')

outFile.write('gene_id\tgenome_id\tcentroid_99\tcentroid_95\tcentroid_90\tcentroid_85\tcentroid_80\tcentroid_75\n')

# make a dictionary with the key=query, value=[query_99, query_95]
centroid_dictionary={}
for line in inFile_99:
    items=line.strip().split('\t')
    record_type=items[0]
    query=items[8]
    target=items[9]
    if record_type=='S':
        # this means that the query = target
        target=query
    if record_type == 'S' or record_type=='H':
        centroid_dictionary[query]=[target]


for line in inFile_95:
    items=line.strip().split('\t')
    record_type=items[0]
    query=items[8]
    target=items[9]
    if record_type=='S':
        # this means that the query = target
        target=query
    if record_type == 'S' or record_type=='H':
        centroid_dictionary[query].append(target)

print centroid_dictionary


for query in centroid_dictionary.keys():
    if len(centroid_dictionary[query]) == 2: # this means that the 99% centroid was also a 95% centroid. The key and both values should match. 
        target_99=centroid_dictionary[query][0]
        target_95=centroid_dictionary[query][1]
    else: # this means that the 99% centroid was clustered into a 95% cluster
        target_99=centroid_dictionary[query][0]
        target_95=centroid_dictionary[target_99][1]
    outFile.write(query + '\tnew_species\t' + target_99 + '\t' + target_95 + '\t' + target_95 + '\t' + target_95 + '\t' + target_95 + '\t' + target_95 + '\n' )

'''
for line in inFile:
    items=line.strip().split('\t')
    record_type=items[0]
    query=items[8]
    target=items[9]
    if record_type=='S':
        # this means that the query = target
        target=query
    if record_type == 'S' or record_type=='H':
        # want to exclude 'C' entries, which are duplicates of S entries
        outFile.write(query + '\tnew_species\t' + target + '\t' + target + '\t' + target + '\t' + target + '\t' + target + '\t' + target + '\n' )

'''
##############################################################
# Update species_info.txt ?

# update genome_info.txt ?



import sys
import os
import os.path
import gzip

timept = sys.argv[1]
################################################################################


# read in the blast output and find the set of genes that have a good match to the pangenome:

inFN='/srv/gsfs0/projects/snyder/ngarud/Usearch_BLAST_prodigal/' + timept + '_usearch_blast.b6'
inFile=open(inFN,'r')

gene_matches=[]

for line in inFile:
    items=line.strip().split('\t')
    query=items[0]
    percent_id=float(items[2])
    if percent_id >=0.95:
        gene_matches.append(query)

# next, read through the fasta file. If a gene is not in gene_matches, then output it to the unmatched file. 
inFN_fasta='/srv/gsfs0/projects/snyder/wenyu/LongReads_10X/IDBA/prodigal_NT/' + timept + '.prodigal.fa'
inFile_fasta=open(inFN_fasta,'r')

outFN='/srv/gsfs0/projects/snyder/ngarud/unannotated_prodigal_genes/' + timept + '_unannotated_prodigal_genes.fa'
outFile=open(outFN, 'w')

output_dna=False
for line in inFile_fasta:
    line=line.strip()
    #check whether the gene is in the gene_matches. 
    if line[0]=='>':
        header=line[1:]
        if header not in gene_matches:
            output_dna=True
            #reformat header so that only the gene ID + timept is outputted
            outFile.write('>' + header.split(' ')[0] + '_'+timept +'\n')
        else:
            output_dna=False
    elif output_dna==True:
         outFile.write(line+'\n')
            

import os
import pandas as pd
import pysam

def get_read(samfile, fastafile, header, start):
    start = int(start)
    end = start+1
    read = []

    for pileupcolumn in samfile.pileup(header, start, end):
        for pileupread in pileupcolumn.pileups:
            if (pileupcolumn.pos == start):
               read.append(str(pileupread.alignment.query_name))

    read = str(read)
    return read

cwd = os.getcwd()

#unzip snp files
os.system('gunzip '+cwd+'/output/*.gz')

#generate bam index files
os.system('samtools index '+cwd+'/temp/genomes.bam '+cwd+'/temp/genomes.bam.bai')

sampath = cwd+'/temp/genomes.bam'
fastapath = cwd+'/temp/genomes.fa'

samfile = pysam.AlignmentFile(sampath, 'rb')
fastafile = pysam.FastaFile(fastapath)

flist = os.listdir(cwd+'/output')

for f in flist:
    spec_name = os.path.splitext(f)[0]
    df = pd.read_csv('~/test_path/output/'+f, delimiter='\t')
    df = df.dropna()
    df['reads'] = df.apply(lambda row: get_read(samfile, fastafile, row['ref_id'], row['ref_pos']), axis=1)
    df.to_csv(spec_name+'_snps_reads.csv', sep='\t')

samfile.close()
fastafile.close()

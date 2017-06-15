import os
import pandas as pd
import pickle
import pysam

def get_read(samfile, fastafile, header, start):
    end = start+1

    samfile = pysam.AlignmentFile(args.bam, "rb")
    fastafile = pysam.FastaFile(args.fasta)
    
    read = []
    
    for pileupcolumn in samfile.pileup(header, start, end):
        for pileupread in pileupcolumn.pileups:
            if (pileupcolumn.pos == args.start):
               read.extend(pileupread.alignment.query_name)

    read = str(read)
    samfile.close()
    fastafile.close()

    return read


#unzip snp files
os.system('gunzip output/*.gz')

#generate bam index files
os.system('samtools index temp/genomes.bam temp/genomes.bam.bai')

flist = os.listdir('output/')

for f in flist:
    spec_name = os.path.splitext(f)[0]
    df = pd.read_csv(args.path+f, delimiter='\t')
    df = df.dropna()
    df['reads'] = df.apply(lambda row: get_read(samfile, fastafile, row['ref_id'], row['ref_pos']), axis=1)
    df.to_csv(spec_name+'_snps_reads.csv', sep='\t')

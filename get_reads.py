import os
import pandas as pd
import pysam

def get_read(samfile, fastafile, header, start, ref_base):
    start = int(start)
    end = start+1
    ref_read = []
    alt_read = []
    err_read = []

    for pileupcolumn in samfile.pileup(header, start, end):
        for pileupread in pileupcolumn.pileups:
            if (pileupcolumn.pos == start):
               if pileupread.query_position == None:
                   err_read.append(str(pileupread.alignment.query_name))
               else:
                   if (pileupread.alignment.query_sequence[pileupread.query_position]==ref_base):
                       ref_read.append(str(pileupread.alignment.query_name))
                   else:
                       alt_read.append(str(pileupread.alignment.query_name))
    ref_read = str(ref_read)
    alt_read = str(alt_read)
    err_read = str(err_read)
    return ref_read, alt_read, err_read

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

lambdafunc = lambda row: pd.Series(get_read(samfile, fastafile, row['ref_id'], row['ref_pos'], row['ref_allele']))

for f in flist:
    spec_name = os.path.splitext(f)[0]
    df = pd.read_csv('~/test_path/output/'+f, delimiter='\t')
    df = df.dropna()
    df[['ref_reads', 'alt_reads', 'err_reads']] = df.apply(lambdafunc, axis=1)
    df.to_csv(spec_name+'_snps_reads.csv', sep='\t')

samfile.close()
fastafile.close()

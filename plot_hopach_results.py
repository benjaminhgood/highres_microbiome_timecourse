import pylab
import numpy
import os
import sys

def calculate_covariance(f1,f2):
    
    f1_avg = f1.mean()
    f2_avg = f2.mean()
    
    return ((f1-f1_avg)*(f2-f2_avg)).sum()
    
    
snp_filename = sys.argv[1]
snp_file = open(snp_filename,"r")
snp_fs = []
for line in snp_file:
    if line.strip()=="":
        continue
        
    items = line.split("\t")
    snp_fs.append( [float(item) for item in items] )
snp_file.close()
snp_fs = numpy.array(snp_fs)

cluster_labels = set()

hopach_filename = sys.argv[2]
hopach_file = open(hopach_filename,"r")
hopach_file.readline()
snp_cluster_map = {}
cluster_snp_map = {}
for line in hopach_file:
    items = line.split()
    snp_idx = long(items[0])-1
    cluster_label = long(items[2])
    snp_cluster_map[snp_idx] = cluster_label
    
    if cluster_label not in cluster_snp_map:
        cluster_snp_map[cluster_label] = []
        
    cluster_snp_map[cluster_label].append(snp_idx)
hopach_file.close()

nodes = numpy.array([snp_cluster_map[snp_idx] for snp_idx in xrange(0,snp_fs.shape[0])])

print len(nodes)

cluster_fs_map = {}
for cluster_label in cluster_snp_map.keys():
    
    anchor_fs = snp_fs[cluster_snp_map[cluster_label][0]]
    if anchor_fs[0]>0.5:
        anchor_fs = 1-anchor_fs

    cluster_fs_map[cluster_label] = [anchor_fs]

    if len(cluster_snp_map[cluster_label]) > 1:
    
        for snp_idx in cluster_snp_map[cluster_label][1:]:
            
            target_fs = snp_fs[snp_idx]
            
            cov = calculate_covariance(anchor_fs, target_fs)
            
            if cov<0:
                
                target_fs = 1-target_fs
                
            cluster_fs_map[cluster_label].append(target_fs)
            
    cluster_fs_map[cluster_label] = numpy.array(cluster_fs_map[cluster_label])
    
print len(cluster_fs_map.keys())


# now plot!
import pylab

pylab.figure(1,figsize=(7,2))
pylab.figure(2) # dummy

ts = numpy.arange(0,snp_fs.shape[1])

for cluster_label in cluster_fs_map.keys():

    cluster_fs = cluster_fs_map[cluster_label]
    print cluster_fs.shape[0]    
    pylab.figure(2)
    line, = pylab.plot([1],[1],'.')
    color = pylab.getp(line,'color')
    
    avg_fs = numpy.median(cluster_fs, axis=0)
    
    pylab.figure(1)
    for snp_idx in xrange(0,cluster_fs.shape[0]):
        fs = cluster_fs[snp_idx]
        
        #if fs[8]>0.7 and fs[7]<0.3:
        #if numpy.fabs(fs-avg_fs).mean() > 0.2:
        if True:
            pylab.plot(ts,fs,'-',linewidth=0.5,alpha=0.5,color=color,zorder=0)
    
    pylab.plot(ts, avg_fs, 'o-', linewidth=2, color=color, zorder=1,markeredgecolor='k')
pylab.figure(1)
pylab.xlabel('Timepoint')
pylab.ylabel('Frequency')
pylab.ylim([-0.05,1.05])
pylab.savefig('cluster_example.png',bbox_inches='tight',dpi=300)
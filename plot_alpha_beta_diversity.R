# goal: plot alpha diversity vs time 

# set the working directory
setwd('~/highres_microbiome_timecourse_data/species')

##################################
# alpha diversity per acccesion: #
##################################  

#read in the data 
data=read.table('relative_abundance.txt.bz2', header=T)

# iterate through each accession (columns) to compute alpha (shannon) diversity for each accession

H_vector=c() #vector to story every accession's shannon diveristy (H=SUM(-p*ln(p)), summing over all species

for (accession in c(2:length(data[1,]))){
    p=as.numeric(data[,accession])
    H=sum((-1)*p*log(p), na.rm=T)  # natural log
    H_vector=c(H_vector,H)
    }

######################################################
# plot the values in H_vector as a function of days  #
###################################################### 

meta=read.table('~/highres_microbiome_timecourse_scripts/highres_timecourse_ids.txt', header=T, sep='\t', row.names=NULL)

names=c("700024086", "700037453", "4026.2", "6037", "4026", "6037.3", "4022","4023", "4024.1","1025", "1022.1", "4021.1A" ,  "1014.2", "1021", "6041" ,"6037.2", "1022", "1023", "SRR2822459", "6038.1","4025.4","4021A")

days=c()
for (name in names){
   time=meta$day[match(name, meta$subject_id)]
   days=c(days,time)
}

tmp2=data.frame(names,days,H_vector)
tmp3=tmp2[order(tmp2$days),]

out_file="~/highres_microbiome_timecourse_analysis/alpha_diversity.pdf"
pdf(out_file,width=8.25,height=4.86, title = "Figure1",paper="special")   
plot(tmp3$days[4:length(tmp3$days)],tmp3$H_vector[4:length(tmp3$days)], col='blue',  main='Alpha diversity vs time', ylab='alpha (Shannon diversity)', xlab='days', pch=20, type='o', ylim=c(2.5,3.6))
points(tmp3$days[1:2],tmp3$H_vector[1:2], col='red', pch=20)
points(tmp3$days[3],tmp3$H_vector[3], col='green', pch=20)  
legend('topright',col=c('blue','red','green'), c('patient','HMP controls','Kuleshov'), pch=20)
dev.off()


## try again rescaling coverage
#coverage=read.table('coverage.txt.bz2', header=T)
#H_vector_2=c() #vector to story every accession's shannon diveristy (H=SUM(-p*ln(p)), summing over all species

#for (accession in c(2:length(data[1,]))){
#    cov=as.numeric(coverage[,accession])  
#    p=cov/sum(cov)
#    H=sum((-1)*p*log(p), na.rm=T)
#    H_vector_2=c(H_vector_2,H)        
#    }




###########################
# compute gamma diversity #
###########################

# first, we want the coverage across all samples (total poplation)
total_coverage=rowSums(sample_coverage[,2:length(sample_coverage[1,])])

# now compute the gamma, which is alpha of the total pop. 
p=total_coverage/sum(total_coverage)
gamma=sum((-1)*p*log(p), na.rm=T)


##########################
# Compute beta diversity
# Beta=gamma/alpha
#########################

beta=gamma/H_vector_sample


###### 
#plot#
######
# Question: do the alpha diveristy computed per accession vs per sample comapare?
# Answer: YES
out_file="~/ben_nandita_hmp_analysis/alpha_diversity.pdf"
pdf(out_file,width=4.25,height=4.86, title = "Figure1",paper="special")   
boxplot(H_vector, H_vector_sample, col='blue', names=c('Per accession', 'Per sample'), main='Alpha values in HMP data', ylab='alpha (Shannon diversity)')
dev.off()

# Question: what does alpha vs beta diversity look like for the samples?
# Answer: beta diversity is a little higher than 1, meaning than within-patient diversity at the species' level is less than the global pop. species' diversity.

out_file="~/ben_nandita_hmp_analysis/alpha_beta_diversity.pdf"
pdf(out_file,width=4.25,height=4.86, title = "Figure1",paper="special")   
boxplot(H_vector_sample, beta, col='blue', names=c('alpha', 'beta'), main='Alpha and beta values in HMP data', ylab='value')
dev.off()

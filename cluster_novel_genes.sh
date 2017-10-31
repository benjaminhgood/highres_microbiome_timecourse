########################################################## 
# September 5 2017
# identify novel gene clusters from Morteza's gene annotation
##########################################################  

# morteza's genes:
/srv/gsfs0/projects/snyder/wenyu/LongReads_10X/IDBA


#######################################
# First we want to find genes that do not match any centroid in the pangenome.
# usearch for blast:
# This does a global alignment
# preferable, but slower than BLAST
####################################

# create a databse with all fasta files for all centroids in the midas db. We will blast new genes against this. 

rm ~/highres_microbiome_timecourse_analysis/novel_centroids/centroids_all_species.ffn.gz
while read species; do
    echo $species
    cat /home/ngarud/midas_db/pan_genomes/${species}/centroids.ffn.gz >> ~/highres_microbiome_timecourse_analysis/novel_centroids/centroids_all_species.ffn.gz     
done < ~/highres_microbiome_timecourse_analysis/novel_centroids/species_list.txt

gunzip ~/highres_microbiome_timecourse_analysis/novel_centroids/centroids_all_species.ffn.gz

# run usearch to blast new genes against the db created above. 

/srv/gsfs0/projects/snyder/wenyu/usearch61 -usearch_global /srv/gsfs0/projects/snyder/wenyu/LongReads_10X/IDBA/prodigal_NT/${timept}.prodigal.fa -db ~/highres_microbiome_timecourse_analysis/novel_centroids/centroids_all_species.ffn -id 0.95 -strand both -blast6out /srv/gsfs0/projects/snyder/ngarud/Usearch_BLAST_prodigal/${timept}_usearch_blast.b6

qsub ~/highres_microbiome_timecourse_scripts/qsub_scripts/Usearch_BLAST_prodigal

##############################################################
# parse blast output
# retain only those genes that have poor matches/no matches
##############################################################  

# python script to parse blast output (run in parallel for all time points)

python parse_prodigal_blast.py $timept
qsub ~/highres_microbiome_timecourse_scripts/qsub_scripts/parse_prodigal_blast

# concatenate all fasta sequences for unknown genes into 1 file accross all time points.# The reason for this is in the next step we will cluster these genes to find non-redundant genes. 
 
rm /srv/gsfs0/projects/snyder/ngarud/unannotated_prodigal_genes/all_timepts_unannotated_prodigal_genes.fa
cat /srv/gsfs0/projects/snyder/ngarud/unannotated_prodigal_genes/* > /srv/gsfs0/projects/snyder/ngarud/unannotated_prodigal_genes/all_timepts_unannotated_prodigal_genes.fa

#####################################################
# cluster the new genes identified in previous step
# usearch
#######################################################

module load usearch/7.0.1090

# test fasta sequence:
#fasta=/home/ngarud/midas_db/marker_genes/phyeco.fa
#usearch -cluster_fast $fasta -id 0.95 -centroids nr.fasta -uc clusters.uc

# run usearch on the genes that have poor/no blast matches
fasta=/srv/gsfs0/projects/snyder/ngarud/unannotated_prodigal_genes/all_timepts_unannotated_prodigal_genes.fa

/srv/gsfs0/projects/snyder/wenyu/usearch61 -cluster_fast $fasta -id 0.99 -centroids /srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_99.fa -uc /srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_99.uc

qsub ~/highres_microbiome_timecourse_scripts/qsub_scripts/Usearch_cluster_prodigal_99

# then cluster again at the 95% level
fasta=/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_99.fa

/srv/gsfs0/projects/snyder/wenyu/usearch61 -cluster_fast $fasta -id 0.95 -centroids /srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_95.fa -uc /srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids_95.uc

qsub ~/highres_microbiome_timecourse_scripts/qsub_scripts/Usearch_cluster_prodigal_95


# python script to parse the uclust output and generate pangenome MIDAS files. 

# files to make:
# i.centroid_functions.txt.gz  -- list of gene names in the centroids.ffn.gz file
# ii.centroids.ffn.gz  -- replicates what is in the fasta file output of usearch
# iii.gene_info.txt.gz – List of all genes with their centroid mapping. Could be comprehensive but may not be necessary (i.e. I could just list the centroids only). 

# data lives here:
/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids.fa
/srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/unique_prodigal_centroids.uc

# create a directory with the new species information for running MIDAS
mkdir /srv/gsfs0/projects/snyder/ngarud/unique_prodigal_genes/files_for_midas/pan_genomes

python ~/highres_microbiome_timecourse_scripts/parse_usearch.py


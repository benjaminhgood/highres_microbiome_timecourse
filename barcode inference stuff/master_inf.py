import os


d = os.listdir('highres_microbiome_timecourse_code/snps_to_track')

for f in d:
    os.system('sbatch run_inf.sh highres_microbiome_timecourse_code/snps_to_track/'+f)

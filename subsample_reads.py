import sys
from random import random

if len(sys.argv)<4:
    print "usage: subsample_reads.py R1.fastq R2.fastq subsampling_fraction"

r1_fastq_filename = sys.argv[1]
r2_fastq_filename = sys.argv[2]
subsampling_fraction = float(sys.argv[3])

r1_filename_items = r1_fastq_filename.split(".")
r2_filename_items = r2_fastq_filename.split(".")

r1_output_filename = ".".join(r1_filename_items[:-1]+["subsampled"]+[r1_filename_items[-1]])

r2_output_filename = ".".join(r2_filename_items[:-1]+["subsampled"]+[r2_filename_items[-1]])

r1_input_file = open(r1_fastq_filename,"r")
r2_input_file = open(r2_fastq_filename,"r")

r1_output_file = open(r1_output_filename,"w")
r2_output_file = open(r2_output_filename,"w")

while True:
    
    # read in next read in both files
    r1_name = r1_input_file.readline()
    r1_sequence = r1_input_file.readline()    
    r1_plus = r1_input_file.readline()
    r1_quality_score = r1_input_file.readline()

    r2_name = r2_input_file.readline()
    r2_sequence = r2_input_file.readline()
    r2_plus = r2_input_file.readline()
    r2_quality_score = r2_input_file.readline()

    if r1_name.strip()=="":
        break
    
    if random() < subsampling_fraction:
    
        # print read
        r1_output_file.write(r1_name)
        r1_output_file.write(r1_sequence)
        r1_output_file.write(r1_plus)
        r1_output_file.write(r1_quality_score)

        r2_output_file.write(r2_name)
        r2_output_file.write(r2_sequence)
        r2_output_file.write(r2_plus)
        r2_output_file.write(r2_quality_score) 

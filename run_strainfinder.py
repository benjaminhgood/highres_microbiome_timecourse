import sys
import os

species_name = sys.argv[1]
Nmax = 7
Lmax = 1000
output_directory = os.path.expanduser("~/strainfinder_output/")

filename_prefix = "%s/%s" % (output_directory, species_name)

os.system('python create_StrainFinderInput.py -o %s --species %s -L %d' % (output_directory, species_name, Lmax))

for N in xrange(2,Nmax+1):

    alignment_filename = filename_prefix+".strainfinder.p" 
    em_filename = filename_prefix+("%d.em.cpickle" % N)
    log_filename = filename_prefix+("%d.log.txt" % N)

    os.system('python StrainFinder.py --aln %s -N %d --max_reps 10 --dtol 1 --ntol 2 --max_time 3600 --converge --em %s --em_out %s --log %s --n_keep 3 --force_update --merge_out --msg' % (alignment_filename,N,em_filename, em_filename,log_filename))
    
    
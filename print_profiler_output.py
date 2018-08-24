import pstats
p = pstats.Stats('profile_output')
p.sort_stats('cumulative').print_stats(10)
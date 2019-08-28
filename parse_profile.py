import pstats

p = pstats.Stats('c_test_run.prof')
p.sort_stats('tottime').print_stats(10)
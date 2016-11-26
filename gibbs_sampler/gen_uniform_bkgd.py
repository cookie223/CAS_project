from gibbs_sampler_lib import *
import sys

if sys.argv[1] == '--help':
	print('gen_uniform_bkgd.py [file] [col]')
	sys.exit()

m = sys.argv[1]
col = int(sys.argv[2])
m = alphabet_to_index(m, col)
n = float(len(m[0]))
for i in m[0].keys():
	print('\t'.join([i, str(1 / n)]))

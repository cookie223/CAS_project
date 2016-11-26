from gibbs_sampler_lib import *
import matplotlib.pyplot as plt
import pickle
import sys

if sys.argv[1] == '--help':
	print('visualize_motif.py [pssm.pckl]')
	sys.exit()

pssm = sys.argv[1]

f = open(pssm, 'rb')
object = pickle.load(f)
score = object['score']
# width = object['width']
# alphabet_dic = object['alphabet_dic']
f.close()

plt.imshow(score.T, cmap='hot',interpolation='nearest', vmax = math.log(2), vmin = 0)
plt.savefig(outname + '_width' + str(width) + '_motif.png')
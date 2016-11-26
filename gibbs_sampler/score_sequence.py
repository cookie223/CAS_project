import sys
import re
import pickle
from gibbs_sampler_lib import *
import matplotlib.pyplot as plt

if sys.argv[1] == '--help':
	print('score_sequence.py [seq_pool] [pssm.pckl]')
	sys.exit()

seqfile = sys.argv[1]
pssm = sys.argv[2]

f = open(pssm, 'rb')
object = pickle.load(f)
score = object['score']
width = object['width']
alphabet_dic = object['alphabet_dic']
f.close()

(seq_pool, lengths, names) = read_fasta(seqfile)

scores = []

for i in seq_pool:
	temp = scoring_seq_by_pssm(i, width, score, alphabet_dic)
	plt.plot(temp)
	plt.hold(True)
	scores.append(temp)
# plt.show()
plt.savefig(seqfile + '_' + pssm + '.png')
import sys
import re
from gibbs_sampler_lib import *
import os
import numpy as np

if sys.argv[1] == '--help':
	print('sequence_generator.py [seed_region] [background_dist] [background_percentage] [alphabet] [mutation_rate]')
	sys.exit()

seeds = sys.argv[1]
background_distf = sys.argv[2]
background_percentage = float(sys.argv[3])
alphabet = sys.argv[4]
mutation_rate = float(sys.argv[5])
temp = alphabet.split('-')
alphabet_file = temp[0]
col = int(temp[1])
alphabet_dic = alphabet_to_index(alphabet_file, col, 'q')
background_dist = read_as_pssm(background_distf, alphabet_dic)

(seed_pool, lengths, names) = read_fasta(seeds)

build_seqs(10, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate)
build_pair_seqs(10, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate)
build_triple_seqs(10, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate)





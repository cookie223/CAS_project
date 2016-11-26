import sys
import re
import numpy as np
from gibbs_sampler_lib import *
import matplotlib.pyplot as plt

if sys.argv[1] == '--help':
	print('gibbs_sampler.py [seqs_fasta] [length_of_motif] [max_step] [alphabet_file-col] [background_dist] [pseudo_count]')
	sys.exit()

fasta = sys.argv[1]
width = int(sys.argv[2])
max_step = int(sys.argv[3])
alphabet = sys.argv[4]
pseudo_count = float(sys.argv[6])
temp = alphabet.split('-')
alphabet_file = temp[0]
col = int(temp[1])
alphabet_dic = alphabet_to_index(alphabet_file, col, 'q')
background_dist = read_as_pssm(sys.argv[5], alphabet_dic)

(seq_pool, lengths, names) = read_fasta(fasta)
k = len(lengths)

## initialization 
pos = []
for i in range(0, len(lengths)):
	# check if width is valid
	if lengths[i] < width:
		print('input width is too big...\nexiting...')
		sys.exit()
	# random draw
	index = np.random.randint(0, lengths[i] - width + 1)
	pos.append([index, i])

index_star = np.random.randint(0, k)
index_star_pre = index_star
current_pssm = construct_pssm(seq_pool, pos, width, index_star, alphabet_dic, pseudo_count, background_dist)
counter = 0


## main body
objs = []
flag = 0
best_pos = pos
best_score = 0
while max_step >= counter:
	counter += 1
	current_pssm = update_pssm(seq_pool, pos, width, index_star, index_star_pre, alphabet_dic, current_pssm, pseudo_count, background_dist)
	# score_across_seq_star = scroing_seq(seq_pool[index_star], current_pssm, alphabet_dic)
	# print([index_star, pos[index_star][0]])
	# print([counter, counter - (counter % 10) * 10])
	if counter % 10 == 1:
		(score_across_seq_star, obj) = compute_objective(seq_pool, pos, width, current_pssm, index_star, alphabet_dic, pseudo_count, background_dist, 'full')
		print(obj)
		objs.append(obj)
		if best_score > obj or flag == 0:
			best_pos = pos
			flag = 1
	else:
		(score_across_seq_star, obj) = compute_objective(seq_pool, pos, width, current_pssm, index_star, alphabet_dic, pseudo_count, background_dist, '')
	
	pos[index_star][0] = np.argmax(score_across_seq_star) # np.argmax(np.random.multinomial(1, score_across_seq_star, size = 1))
	index_star_pre = index_star
	index_star = draw_from(index_star, k)

	
	# print([index_star_pre, pos[index_star_pre][0]])

for i in range(len(best_pos)):
	print(names[i] + ' ----- ' + str(best_pos[i][0]))





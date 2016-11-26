import sys
import re
import numpy as np
from gibbs_sampler_lib import *
import matplotlib.pyplot as plt
import pickle

if sys.argv[1] == '--help':
	print('gibbs_sampler.py [seqs_fasta] [length_of_motif] [max_step] [alphabet_file-col] [background_dist] [pseudo_count] [repeat_times] [outname]')
	sys.exit()

fasta = sys.argv[1]
width = int(sys.argv[2])
max_step = int(sys.argv[3])
alphabet = sys.argv[4]
pseudo_count = float(sys.argv[6])
temp = alphabet.split('-')
alphabet_file = temp[0]
col = int(temp[1])
repeat_times = int(sys.argv[7])
outname = sys.argv[8]
alphabet_dic = alphabet_to_index(alphabet_file, col, 'q')
background_dist = read_as_pssm(sys.argv[5], alphabet_dic)

(seq_pool, lengths, names) = read_fasta(fasta)
k = len(lengths)


my_plot = []


## initialization 
for e in range(repeat_times): 
	flag = 0
	best_pos = []
	best_score = 0
	
	objs = []
	print('--- repeat ' + str(e + 1) + ' ----')
	print(sys.argv)
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
	# print('construct' + str(sum(current_pssm[0][0])))

	## main body
	
	while max_step >= counter:
		counter += 1
		current_pssm = update_pssm(seq_pool, pos, width, index_star, index_star_pre, alphabet_dic, current_pssm, pseudo_count, background_dist)
		# print(sum(current_pssm[0][0]))

		if counter % 1 == 0:
			(score_across_seq_star, obj, model) = compute_objective(seq_pool, pos, width, current_pssm, index_star, alphabet_dic, pseudo_count, background_dist, 'full')
			# print(obj)
			objs.append(obj)
			if best_score < obj or flag == 0:
				best_score = obj
				best_pos = pos
				best_model = model
				flag = 1
		else:
			(score_across_seq_star, obj) = compute_objective(seq_pool, pos, width, current_pssm, index_star, alphabet_dic, pseudo_count, background_dist, '')
		# print(sum(current_pssm[0][0]))

		pos[index_star][0] = np.argmax(score_across_seq_star) # np.argmax(np.random.multinomial(1, score_across_seq_star, size = 1))
		index_star_pre = index_star
		index_star = draw_from(index_star, k)
		# index_star = weighted_draw_from(index_star, k, current_pssm)

		
		# print([index_star_pre, pos[index_star_pre][0]])
	if len(my_plot) == 0:
		my_plot = np.array(objs)
	else:
		my_plot = np.vstack((my_plot, objs))
	f = open(outname + '_width' + str(width) + '_' + str(e) + '.pckl', 'wb')
	(bpssm_score, bpssm_prob) = count_to_prob(model, pseudo_count, alphabet_dic, background_dist)
	forsave = {}
	forsave['count'] = model
	forsave['score'] = bpssm_score
	forsave['prob'] = bpssm_prob
	forsave['background_dist'] = background_dist
	forsave['pseudo_count'] = pseudo_count
	forsave['alphabet_dic'] = alphabet_dic
	forsave['width'] = width
	pickle.dump(forsave, f)
	f.close()
	

	for i in range(len(best_pos)):
		print(seq_pool[i][best_pos[i][0] : best_pos[i][0] + width] + ' <---> ' + str(best_pos[i][0]))

plt.plot(my_plot.T)
plt.savefig(outname + '_width' + str(width) + '_curve.png')





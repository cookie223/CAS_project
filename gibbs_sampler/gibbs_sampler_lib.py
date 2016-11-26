import re
import numpy as np
import math
import copy

def alphabet_to_index(alphabet_file, col, mode):
	dic_forward = {}
	dic_backward = {}
	t = open(alphabet_file, 'r')
	t = t.readlines()
	if mode != 'q':
		print('reading alphabet')
	counter = 0
	for i in t:
		temp = i.split(' ')
		temp = filter(None, temp)
		now = temp[col - 1]
		now = re.sub('\n', '', now)
		dic_forward[now] = counter
		dic_backward[counter] = now
		if mode != 'q':
			print(' '.join([now, 'is encoded as', str(counter)]))
		counter += 1
	dic = (dic_forward, dic_backward)
	return dic
def construct_pssm(seq_pool, pos, width, index_star, dic, pseudo_count, background_dist):
	temp = []
	for i in range(len(seq_pool)):
		if i != index_star:
			temp.append(seq_pool[i][pos[i][0] : pos[i][0] + width])
	pos_string = []
	for i in range(width):
		col_string = ''
		for j in range(len(temp)):
			col_string += temp[j][i]
		pos_string.append(col_string)
	pssm_count = count_alphabet(pos_string, dic, pseudo_count)
	(pssm_score, pssm_prob) = count_to_prob(pssm_count, pseudo_count, dic, background_dist)
	return (pssm_count, pssm_score, pssm_prob)

def count_alphabet(pos_string, dic, pseudo_count):
	pssm_count = []
	for i in pos_string:
		pos_count = []
		for j in range(len(dic[1].keys())):
			result = re.findall(dic[1][j], i)
			pos_count.append(len(result))
		pssm_count.append(pos_count)
	return pssm_count

def count_to_prob(pssm_count, pseudo_count, dic, background_dist):
	pssm_score = []
	pssm_prob = []
	n = len(dic[0].keys())
	sums = np.sum(pssm_count, 1) + float(n * pseudo_count)
	temp_pssm = np.array(pssm_count) + float(pseudo_count)
	pssm_prob = temp_pssm.transpose() / sums
	# print(pssm_prob)
	# print(background_dist)
	pssm_score = np.log(pssm_prob.transpose() / background_dist)
	# pssm_score = pssm_score
	pssm_prob = np.log(pssm_prob.transpose())
	return (pssm_score, pssm_prob)

def read_as_pssm(filein, dic):
	f = open(filein, 'r')
	f = f.readlines()
	temp = {}
	back = []
	for i in f:
		i = re.sub('\n', '', i)
		i = i.split('\t')
		temp[i[0]] = float(i[1])
	for i in range(len(temp.keys())):
		back.append(temp[dic[1][i]])
	return np.array(back)

def update_pssm(seq_pool, pos, width, index_star, index_star_pre, dic, current_pssm, pseudo_count, background_dist):
	pssm_count = current_pssm[0]
	pre_in = seq_pool[index_star_pre][pos[index_star_pre][0] : pos[index_star_pre][0] + width]
	now_out = seq_pool[index_star][pos[index_star][0] : pos[index_star][0] + width]
	for i in range(width):
		# print('in=' + str(sum(pssm_count[i])))
		pssm_count[i][dic[0][pre_in[i]]] = pssm_count[i][dic[0][pre_in[i]]] + 1
		pssm_count[i][dic[0][now_out[i]]] = pssm_count[i][dic[0][now_out[i]]] - 1
		# print('out=' + str(sum(pssm_count[i])))
	(pssm_score, pssm_prob) = count_to_prob(pssm_count, pseudo_count, dic, background_dist)
	return (pssm_count, pssm_score, pssm_prob)

def scroing_seq(seq, current_pssm, dic):
	w = len(current_pssm[0])
	scores = []
	# print(current_pssm[2])
	for i in range(len(seq) - w):
		substr = seq[i : i + w]
		score = 0
		for j in range(len(substr)):
			score = logadd(score, current_pssm[2][j][dic[0][substr[j]]])
			# print(substr + '-' + substr[j] + '++' + str(current_pssm[2][j][dic[0][substr[j]]]) + '  ' + str(score))
		scores.append(score)

		# print(str(score) + '--' + substr)
	scores = np.exp(np.array(scores))
	scores = scores / np.sum(scores)
	return scores

def draw_from(index_star, k):
	temp = np.random.randint(0, k - 1)
	if temp >= index_star:
		temp += 1
	return temp

def logadd(x, y):
	return max(x, y) + math.log1p(math.exp( - abs(x - y)))

def read_fasta(filename):
	seq_file = open(filename, 'r')
	seqs = seq_file.readlines()
	seq_pool = []
	flag = 0
	now_seq = []
	lengths = []
	names = []
	for i in seqs:
		if i[0] == '>':
			names.append(re.sub('\n', '', i[1 : ]))
			flag = 1
			if len(now_seq) > 0:
				seq_pool.append(now_seq)
				lengths.append(len(now_seq))
			continue
		else:

			if flag == 1:
				now_seq = re.sub('\n', '', i)
				flag = 0
			else:
				now_seq += re.sub('\n', '', i)
	seq_pool.append(now_seq)
	lengths.append(len(now_seq))
	return (seq_pool, lengths, names)

def build_seqs(seqnum, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate):
	for i in range(len(seed_pool)):
		# temp = open(foldername + '/' + str(i), 'w')
		seq = seed_pool[i]
		blength = int(len(seq) * background_percentage)
		for j in range(seqnum):
			back_seq = ramdon_seq(blength, background_dist, alphabet_dic)
			random_pos = np.random.randint(0, blength)
			full_seq = back_seq[ : random_pos] + random_flip(seq, alphabet_dic, mutation_rate) + back_seq[random_pos : ]
			print('> single ' + str(i) + ' : ' + str(j) + ' ' + str(random_pos))
			print(full_seq)
			# temp.write('>' + str(i) + ':' str(j) + str(random_pos) + '\n')
			# temp.write(full_seq + '\n')

def build_pair_seqs(seqnum, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate):
	for i in range(len(seed_pool)):
		for j in range(i + 1, len(seed_pool)):
			seqi = seed_pool[i]
			seqj = seed_pool[j]
			blength = int(len(seqi + seqj) * background_percentage)
			for k in range(seqnum):
				back_seq = ramdon_seq(blength, background_dist, alphabet_dic)
				random_pos1 = np.random.randint(0, blength)
				random_pos2 = np.random.randint(0, blength)
				randomi = min(random_pos1, random_pos2)
				randomj = max(random_pos1, random_pos2)
				full_seq = back_seq[ : randomi] + random_flip(seqi, alphabet_dic, mutation_rate) + back_seq[randomi : randomj] + random_flip(seqi, alphabet_dic, mutation_rate) + back_seq[randomj : ]
				print('> double ' + str(i) + '-' + str(j) + ' : ' + str(k) + ' ' + str(random_pos1) + '-' + str(random_pos2))
				print(full_seq)	

def build_triple_seqs(seqnum, seed_pool, background_percentage, background_dist, alphabet_dic, mutation_rate):
	for i in range(seqnum):
		full_seq = ''
		random_poss = []
		for j in range(3):
			seqj = seed_pool[j]
			blength = int(len(seqj) * background_percentage)
			back_seq = ramdon_seq(blength, background_dist, alphabet_dic)
			random_pos = np.random.randint(0, blength)
			pre_full_length = len(full_seq)
			full_seq = full_seq + back_seq[ : random_pos] + random_flip(seqj, alphabet_dic, mutation_rate) + back_seq[random_pos : ]
			random_poss.append(str(random_pos + pre_full_length))
		print('> triple ' + '0-1-2 : ' + str(i) + ' ' + '-'.join(random_poss))
		print(full_seq)

def ramdon_seq(blength, background_dist, alphabet_dic):
	out = ''
	for i in range(blength):
		index = np.argmax(np.random.multinomial(1, background_dist, size = 1))		
		out += alphabet_dic[1][index]
	return out

def random_flip(seqj, alphabet_dic, mutation_rate):
	for i in range(len(seqj)):
		mut = np.random.rand()
		if mut < mutation_rate:
			newindex = np.random.randint(0, len(alphabet_dic[0].keys()))
			temp = list(seqj)
			temp[i] = alphabet_dic[1][newindex]
			seqj = ''.join(temp)
	return seqj

def compute_objective(seq_pool, pos, width, current_pssm, index_star, dic, pseudo_count, background_dist, mode):
	score_across_seq_star = scroing_seq(seq_pool[index_star], current_pssm, dic)
	if mode != 'full':
		return (score_across_seq_star, 0)
	pos[index_star][0] = np.argmax(score_across_seq_star)
	pssm_count = copy.deepcopy(current_pssm[0])
	pre_in = seq_pool[index_star][pos[index_star][0] : pos[index_star][0] + width]
	for i in range(width):
		pssm_count[i][dic[0][pre_in[i]]] = pssm_count[i][dic[0][pre_in[i]]] + 1
	(pssm_score, pssm_prob) = count_to_prob(pssm_count, pseudo_count, dic, background_dist)
	score = 0
	for i in pos:
		substr = seq_pool[i[1]][i[0] : i[0] + width]
		for j in range(len(substr)):
			score += pssm_score[j][dic[0][substr[j]]]
	return(score_across_seq_star, score, pssm_count)

def gen_true_pos(names, true_id, start):
	true_pos = []
	counter = 0
	for i in names:
		i = i.split(' ')
		i = i[-1]
		i = i.split('-')
		i = int(i[true_id])
		true_pos.append([i + start, counter])
		counter += 1
	return true_pos

def scoring_seq_by_pssm(seq, width, bpssm_score, alphabet_dic):
	score = []
	for i in range(len(seq) - width):
		subseq = seq[i : i + width]
		nowscore = scoring_substr_by_pssm(subseq, bpssm_score, alphabet_dic)
		score.append(nowscore)
	return score

def scoring_substr_by_pssm(substr, bpssm_score, dic):
	score = 0
	for j in range(len(substr)):
		# print(score)
		score += bpssm_score[j][dic[0][substr[j]]]
	return score



# def weighted_draw_from(index_star, k, current_pssm):
	


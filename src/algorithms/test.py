import os
import random
import string

code = ['A', 'T', 'C', 'G']
str_length = [1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000]
range_length = [10, 20, 40, 80, 160, 320, 640]
word_num = [1, 2, 3, 4, 6, 8]
str_num = [100, 200, 400, 800]
word_length = 10
code = "ATCG"
random_codes = lambda x, y: ''.join([random.choice(x) for i in range(y)])

def insert (source_str, insert_str, pos):
    return source_str[:pos]+insert_str+source_str[pos:]


def generate_word_range(word_length, word_num, range_length):
	str_seed = []
	str_seq = []
	counter = 1
	word = random_codes(code, word_length)
	while counter < word_num:
		str_seed.append(word)
		str_seq.append(word)
		str_seed.append(random_codes('U', range_length))
		str_seq.append(random_codes('Y', range_length))
		counter += 1

	str_seed.append(word)
	str_seq.append(word)

	return ''.join(str_seed), ''.join(str_seq)



def test_str_length():

	for i in str_length:
		filename_pos = "./test_files/str_len/" + str(i) + "_pos.fas"
		filename_neg = "./test_files/str_len/" + str(i) + "_neg.fas"
		print filename_pos
		file = open(filename_pos, 'w+')
		insert_str_seed, insert_str_seq = generate_word_range(word_length,word_num[1],range_length[1])

		for j in range(str_num[1]):
			file.write( ">str_len_" + str(i) + "_test_" + str(j) + "\n")
			source_str = random_codes(code, i) + "\n"
			index = random.randint(0, len(source_str)-1)
			if(j==0):
				seq = insert(source_str, insert_str_seed, index)
			else:
				seq = insert(source_str, insert_str_seq, index)
			file.write(seq)
		file.close()
		command = "fasta-shuffle-letters " + filename_pos + ' ' + filename_neg
		os.popen(command)


def test_range():

	for i in range_length:
		filename_pos = "./test_files/range/" + str(i) + "_pos.fas"
		filename_neg = "./test_files/range/" + str(i) + "_neg.fas"
		print filename_pos
		file = open(filename_pos, mode='w')
		insert_str_seed, insert_str_seq = generate_word_range(word_length,word_num[5],i)

		for j in range(str_num[2]):
			file.write( ">range_" + str(i) + "_test_" + str(j) + "\n")
			source_str = random_codes(code, str_length[1]) + "\n"
			index = random.randint(0, len(source_str)-1)
			if(j==0):
				seq = insert(source_str, insert_str_seed, index)
			else:
				seq = insert(source_str, insert_str_seq, index)
			file.write(seq)
		file.close()
		command = "fasta-shuffle-letters " + filename_pos + ' ' + filename_neg
		os.popen(command)		

def test_str_num():

	for i in str_num:
		filename_pos = "./test_files/str_num/" + str(i) + "_pos.fas"
		filename_neg = "./test_files/str_num/" + str(i) + "_neg.fas"
		print filename_pos
		file = open(filename_pos, mode='w')
		insert_str_seed, insert_str_seq = generate_word_range(word_length,word_num[2],range_length[1])

		for j in range(i):
			file.write( ">str_num_" + str(i) + "_test_" + str(j) + "\n")
			source_str = random_codes(code, str_length[1]) + "\n"
			index = random.randint(0, len(source_str)-1)
			if(j==0):
				seq = insert(source_str, insert_str_seed, index)
			else:
				seq = insert(source_str, insert_str_seq, index)
			file.write(seq)
		file.close()
		command = "fasta-shuffle-letters " + filename_pos + ' ' + filename_neg
		os.popen(command)	

def test_word_num():

	for i in word_num:
		filename_pos = "./test_files/word_num/" + str(i) + "_pos.fas"
		filename_neg = "./test_files/word_num/" + str(i) + "_neg.fas"
		print filename_pos
		file = open(filename_pos, mode='w')
		insert_str_seed, insert_str_seq = generate_word_range(word_length, i, range_length[6])

		for j in range(str_num[3]):
			file.write( ">word_num_" + str(i) + "_test_" + str(j) + "\n")
			source_str = random_codes(code, str_length[6]) + "\n"
			index = random.randint(0, len(source_str)-1)
			if(j==0):
				seq = insert(source_str, insert_str_seed, index)
			else:
				seq = insert(source_str, insert_str_seq, index)
			file.write(seq)
		file.close()
		command = "fasta-shuffle-letters " + filename_pos + ' ' + filename_neg
		os.popen(command)	


def limit_test():

	filename_pos = "./test_files/limit_test/" + "test_pos.fas"
	filename_neg = "./test_files/limit_test/" + "test_neg.fas"
	print filename_pos
	file = open(filename_pos, mode='w')
	insert_str_seed, insert_str_seq = generate_word_range(word_length, word_num[5], range_length[6])

	for j in range(str_num[3]):
		file.write( ">limit_test" + "\n")
		source_str = random_codes(code, str_length[8]) + "\n"
		index = random.randint(0, len(source_str)-1)
		if(j==0):
			seq = insert(source_str, insert_str_seed, index)
		else:
			seq = insert(source_str, insert_str_seq, index)
		file.write(seq)
	file.close()
	command = "fasta-shuffle-letters " + filename_pos + ' ' + filename_neg
	os.popen(command)

test_str_length()

test_range()

test_str_num()

test_word_num()

limit_test()
#include "stdio.h"
#include "string.h"
#include <time.h>
#include <stdlib.h>

#define MAX_RANDOM 5000
#define WORD_LENGTH 10
#define MAX_WORD_NUM 100
#define TRUE 1
#define FALSE 0

// seq = (0-MAX_RANDOM)[(WORD_LENGTH + RANDOM_OFFSET) * WORDS_NUM](0-MAX_RANDOM)
char code[] = {'A', 'C', 'G', 'T'};
#define CODE_NUM 4

void generateWord(char *buffer, int length){
	
	for (int i=0; i<length; i++) buffer[i] = code[rand() % CODE_NUM];
	buffer[length] = '\0';
	// printf("buffer = %s", buffer);

}

void generate_random_seq(int pre_num , int suf_num, char *word, int word_num, int offset, char *seq){
	int index = 0;

	for(index; index<pre_num; index++) {
		
		seq[index] = code[rand() % CODE_NUM];
	}
	for(int i=0; i<word_num; i++){
	
	memcpy(seq + index, word, WORD_LENGTH );
	index = index + WORD_LENGTH;
		for(int j=0; j<offset; index++, j++){
			seq[index] = code[rand() % CODE_NUM];
		}
	}

	for(int j=0; j<suf_num; index++, j++) seq[index] = code[rand() % CODE_NUM];
	
	seq[index] = '\0';
	// printf("seq = %s\n", seq );
}

void test(int str_num, int min_support_num, int word_num, char *filename){
//./seed_w test0_pos.fas test0_pos.fas -m save_motifs.xml --min_support 1.0 --word_length 5
  FILE *fh;
  srand(time(NULL));
  int pre_num, suf_num;
  
  char prefix[MAX_RANDOM+1] = {0};
  char suffix[MAX_RANDOM+1] = {0};
  char word[WORD_LENGTH+1] = {0};
  char offset[22] = {0};
  char seq[MAX_RANDOM*2 + (WORD_LENGTH+22)*MAX_WORD_NUM	+1] = {0};
  int offset_num = (rand() % 20) + 1;// its the length between words(range 1 - 21)
	
  fh = fopen( filename, "w" );
  generateWord(word, WORD_LENGTH);
  printf("word = %s\n", word);	

  for(int i=0; i<str_num; i++){
  	fprintf(fh, ">test0\n" );

  	pre_num = rand() % MAX_RANDOM;//random length of prefix
	suf_num = rand() % MAX_RANDOM;//random length of suffix
	// printf("pre_num = %d\n", pre_num);
	// printf("suf_num = %d\n", suf_num);
  	// printf("seq_ori_length = %d, word_insert_index = %d\n", seq_ori_length, word_insert_index);

	if(i<min_support_num)generate_random_seq(pre_num, suf_num, word, word_num, offset_num, seq);
	else generate_random_seq(pre_num, suf_num, word, 0, offset_num, seq);
	// printf("seq  =%s\n", seq);
	fprintf(fh, "%s", seq );
  	fprintf(fh, "\n" );	

  }

  fclose( fh );
}



void main(int argc, char *argv[]){
	// parameters =  str_num, min_support_num, word_num, filename

	// test(10000, 10000, 1, "test0_pos.fas");

	// test(10000, 6000, 1, "test1_pos.fas");

	// test(10000, 5000, 2, "test2_pos.fas");

	// test(10000, 5000, 3, "test3_pos.fas");

	// test(10000, 5000, 3, "test4_pos.fas");

	// test(10000, 100, 3, "test4_neg_1.fas");

	// test(10000, 1000, 3, "test4_neg_2.fas");

	// test(10000, 2000, 3, "test4_neg_3.fas");

	// test(10000, 4000, 3, "test4_neg_4.fas");

	// test(10000, 8000, 3, "test4_pos_5.fas");
	
	// test(10000, 8000, 3, "test4_neg_5.fas");

	// test(10000, 5000, 3, "test5_pos.fas");

	// test(10000, 100, 3, "test5_neg.fas");
}
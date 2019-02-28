/*
* 
* @Author: Weipeng Xu
* @Date:   2019-02-26 16:05:39
* @Last Modified by:   Weipeng Xu
* @Last Modified time: 2019-02-27 21:19:12
*
* This copyrighted source code is freely distributed under the terms
* of the GNU General Public License. 
* See the files COPYRIGHT and LICENSE for details.
*/


#include "libdev.h"
#include "seq.h"
#include "vector.h"
#include "libvtree.h"
#include "match_w.h"
#include "seed_w.h"
#include "kfunc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <stdbool.h>

int NUM_INPUT_STRING;     // Number of Input Strings under Current Mode(POS/NEG)
int NUM_INPUT_STRING_POS; // Number of Positive Input Strings
int NUM_INPUT_STRING_NEG; // Number of Negative Input Strings

/*****************************************************************
 * display_banner - displays the version of the program          *
 *****************************************************************/

static void
display_banner(parameters_t *params)
{
  printf("%s - Identification of consensus word-based motifs inference\n\n", params->version);

  printf("Copyright (C) 2019 University of Ottawa\n");
  printf("All Rights Reserved\n");
  printf("\n");
  printf("This program is distributed under the terms of the\n");
  printf("GNU General Public License. See the source code\n");
  printf("for details.\n\n");
}

/*****************************************************************
 * display_usage - general instructions for running the program  *
 *****************************************************************/

static char usage[] = "\
Usage: seed_w [options] file_pos file_neg\n\
where file is a FASTA file that contains k input DNA/RNA sequences.\n\
\n\
Options:\n\
     --min_support <n>         (default 0.70)\n\
     --word_length <n>         (default 3)\n\
  -m --match_file <file>       (no default)\n\
  -v --version\n\
  -h --help\n\
\n\
";


static void
display_usage_and_exit()
{
  printf("%s", usage);
  exit(EXIT_SUCCESS);
}

/*****************************************************************
 * free_expression -                                             *
 *****************************************************************/

void dev_free_expression(expression *e)
{
  expression *tmp;
  
  while (e != NULL)
  {
    tmp = e;

    e = e->next;

    dev_free(tmp);
  }
}

/*****************************************************************
 * free_motif -                                                  *
 *****************************************************************/

void dev_free_motif(motif *m)
{
  // print_motif(m);
  dev_free_expression(m->head);
  dev_free(m);
}

/*****************************************************************
 * free_array -                                             *
 *****************************************************************/
void dev_free_vtree_array( void**p, int size, int func)
{
  int i;
  //free vtrees 
  for (i = 0; i < size; i++)
    vtree_free(p[i]);

  //free the array
  if(func == TRUE)dev_free(p);
}


/*****************************************************************
 * param_init -                                                  *
 *****************************************************************/

static void
param_init(parameters_t *params)
{
  params->filename_pos = FILENAME_POS;
  params->filename_neg = FILENAME_NEG;
  params->match_file = MATCH_FILE;
  params->min_support = MIN_SUPPORT;
  params->word_length = WORD_LENGTH;
  params->version = VERSION;

  params->match_count = 0;
}

/*****************************************************************
 * process_argv - process the command line arguments             *
 *****************************************************************/

static void
process_argv(int argc, char *argv[], parameters_t *params)
{
  int i = 1;

  while (i < argc)
  {

    if (argv[i][0] != '-')
    {

      if (params->filename_pos != NULL)
        dev_die("not a valid argument %s", argv[i]);

      params->filename_pos = argv[i];

      i++;

      if (params->filename_neg != NULL)
        dev_die("not a valid argument %s", argv[i]);

      params->filename_neg = argv[i];
    }

    else if (strcmp("--min_support", argv[i]) == 0)
    {

      params->min_support = dev_parse_float(argv[++i]);
    }

    else if (strcmp("--word_length", argv[i]) == 0)
    {

      params->word_length = dev_parse_int(argv[++i]);
    }
    
    else if ( strcmp("-v", argv[ i ] ) == 0 || strcmp( "--version", argv[ i ] ) == 0 ) {

      display_banner( params );
      exit( EXIT_SUCCESS );
    }
    
    else if (strcmp("-h", argv[i]) == 0 || strcmp("--help", argv[i]) == 0)
    {
      display_usage_and_exit();
    }

    else if (strcmp("-m", argv[i]) == 0 || strcmp("--match_file", argv[i]) == 0)
    {
      params->match_file = argv[++i];
    }
    else
    {
      dev_die("not a valid argument %s", argv[i]);
    }

    i++;
  }

  /* validating */
  if (params->filename_pos == NULL)
    dev_die("please specify an input pos file");

  if (params->filename_neg == NULL)
    dev_die("please specify an input neg file");
}



/*****************************************************************
 * display2 - prints out the various tables                      *
 *****************************************************************/

static void

display2(char *s, dstring_t *ds, vtree_t *v)
{
  printf("\n");

  printf("      s  ds suf lcp  \n");

  for (int i = 0; i < ds->length; i++)
  {

    char c;

    if (dev_isspecial(&lowercase, ds->text[i]))
      c = '$';
    else
      c = dev_decode(&lowercase, ds->text[i]);

    printf(" %2d ", i);
    printf("[%2c]", c);
    printf("[%2d]", ds->text[i]);
    printf("[%2d]", v->suftab[i]);
    printf("[%2d]", v->lcptab[i]);
    printf(" ");

    for (int j = v->suftab[i]; j < v->length; j++)
    {

      if (dev_isspecial(&lowercase, ds->text[j]))
        c = '$';
      else
        c = dev_decode(&lowercase, ds->text[j]);

      printf("%c", c);
    }

    printf("\n");
  }
}

/*****************************************************************
 * print_motif - print all elements inside motif m               *
 *****************************************************************/

static void
print_motif(char *seed, motif *m, parameters_t *params)
{
  
  expression *e = (expression *) m->head;
  int index = 0;
  char dest[params->word_length + 1];
  printf("num_match_pos = %d\n", m->num_match_pos);
  printf("num_match_neg = %d\n", m->num_match_neg);
  
  if(e == NULL)dev_die("expression unknow");
  
  while (e != NULL)
  {
    if (e->type == enum_word)
    {
      printf("word->type = %d (0 = word,  1 = range)\n", e->type);
      strncpy(dest, &seed[e->begin_index], params->word_length);
      dest[params->word_length] = '\0';
      printf("params->word_length = %d\n", params->word_length);
      printf("word = %s\n", dest);
      printf("expression index = %d\n", ++index);
      printf("word->index = %d\n", e->begin_index);
      printf("word->length = %d\n", e->length);
      e = (expression *) e->next;
    }
    else if (e->type == enum_range)
    {
      printf("range->type = %d (0 = word,  1 = range)\n", e->type);
      printf("expression index = %d\n", ++index);
      printf("range->index = %d\n", e->begin_index);
      printf("range->length = %d\n", e->length);
      e = (expression *) e->next;
    }
    else dev_die("expression unknow");
  }

  printf("\n");
}

/*****************************************************************
 * print_seq - print all sequences                               *
 *****************************************************************/

void print_seq(vtree_t *vtree[])
{
  for (int i = 0; i < NUM_INPUT_STRING; i++)
  {
    //print all input strings
    printf("input string %d = ", i);
    for (int j = 0; j < vtree[i]->length - 1; j++)
    {

      printf("%c", dev_decode(&uppercase, vtree[i]->text[j]));
    }
    printf("\n");
  }
}

/*****************************************************************
 * words_enumerate - enumerate all words that has enough support *
 *****************************************************************/

vector_t *
words_enumerate(char seed[], vtree_t *vtree[], parameters_t *params)
{
  //Enumerate all possible motifs made of one word
  //Return a list of motifs(one word motifs)

  printf("length of Seed Sequence = %ld\n", strlen(seed));

  int L = (int)strlen(seed); //SEED_LENGTH

  vector_t *motif_list = dev_new_vector2(VECTOR_INI_CAP, INCREMENT); //(initial capacity, increment)
  dstring_t *pattern;
  
  int cur_support;

  for (int seed_index = 0; seed_index < L - params->word_length + 1; seed_index++)
  {
    char dest[params->word_length + 1];
    cur_support = 1;

    strncpy(dest, &seed[seed_index], params->word_length);
    dest[params->word_length] = '\0';
    pattern = dev_digitalize(&uppercase, dest);
    pattern->length = pattern->length - 1; //This is necessary

    for (int str_index = 1; str_index < NUM_INPUT_STRING; str_index++)
    {
      if (vtree_find_exact_match_new(vtree[str_index], pattern))
        {
          cur_support++;
        }
    }
    
    dev_free_dstring(pattern);
    
    if (cur_support >= params->min_support_num)
    {
      // Found motifs with minimun supports, add to list
      motif *motif_w = (motif *)dev_malloc(sizeof(motif));
      expression *word = (expression *)dev_malloc(sizeof(expression));
      word->begin_index = seed_index;
      word->length = params->word_length;
      word->next = NULL;
      word->type = enum_word;

      motif_w->head = word;
      //Since this motif only contains one word, so the end == head == word
      
      motif_w->num_match_pos = cur_support;
      motif_w->num_match_neg = 1;
      motif_w->end = word;
      motif_w->words_count = 1;//in here is all one word motifs

      dev_vector_add(motif_list, motif_w); //add motif to list
    }
  }

  return motif_list;
}

/*****************************************************************
 * motif_concatenate - concatenate two motifs into one           *
 *****************************************************************/

motif *motif_concatenate(motif *motif_a, motif *motif_b)
{
  // Concatenate two motifs into one, a range will be added in between
  motif *new_motif = (motif *)dev_malloc(sizeof(motif));
  new_motif->head = NULL;
  new_motif->end = NULL;
  new_motif->num_match_pos = 0;
  new_motif->num_match_neg = 0;

  if (motif_a != NULL && motif_b != NULL){
  expression *temp_expression = motif_a->head;
  expression *pre_expression = NULL;
  expression *new_expression = NULL;
  

  if(temp_expression == NULL)dev_die("motif_a unknow");
  while (temp_expression != NULL)
  {
    //Copy info from motfi_a
    new_expression = (expression *)dev_malloc(sizeof(expression));
    new_expression->begin_index = temp_expression->begin_index;
    new_expression->length = temp_expression->length;
    new_expression->type = temp_expression->type;
    //Previous new_expression->next parameter was added at next iteration
    if (pre_expression != NULL)
      pre_expression->next = new_expression;

    if(new_motif->head == NULL)new_motif->head = new_expression;

    pre_expression = new_expression;
    temp_expression = temp_expression->next;
  }

  temp_expression = motif_b->head;

  //insert a range in between
  expression *range = (expression *)dev_malloc(sizeof(expression));
  range->begin_index = pre_expression->begin_index + pre_expression->length;
  range->length = temp_expression->begin_index - range->begin_index;
  range->type = enum_range;
  pre_expression->next = range;
  pre_expression = range;

  while (temp_expression != NULL)
  {
    //Copy info from motfi_b
    new_expression = (expression *)dev_malloc(sizeof(expression));
    new_expression->begin_index = temp_expression->begin_index;
    new_expression->length = temp_expression->length;
    new_expression->type = temp_expression->type;
    if (pre_expression != NULL)
      pre_expression->next = new_expression;

    pre_expression = new_expression;
    temp_expression = temp_expression->next;
  }
  pre_expression->next = NULL;

  new_motif->end = new_expression;
  }
  else dev_die("motif_a or motif_b is NULL");
  return new_motif;
}

int has_enough_support(dstring_t *ds_seed, vtree_t *vtree[], motif *query_motif, int *support, parameters_t *params)
{
  //Determine whether a motif has enough support among all the input strings
  int cur_support = 1; 
  interval2_t *i0;

  for (int i = 1; i < NUM_INPUT_STRING; i++)
  {
    i0 = new_interval2( 0, vtree[i]->length );

    if (match_motif(ds_seed, i0, vtree[i], query_motif) == 1)
    {
      cur_support++;
    }

     dev_free( i0 );
  }
  *support = cur_support;

  if (cur_support >= params->min_support_num)
  {
    return TRUE;
  }

  else
  {
    return FALSE;
  }
}

void get_neg_support(dstring_t *ds_seed, vtree_t *vtree[], vector_t *motif_list, parameters_t *params)
{
  //Travel through negative set and record support for each motifs
  for (int i = 0; i < dev_vector_size(motif_list); i++)
  {
    //Record num_match_neg for each motifs
    motif *query_motif = (motif *)dev_vector_get(motif_list, i);
    int cur_support = 0;
    
    has_enough_support(ds_seed, vtree, query_motif, &cur_support, params);

    query_motif->num_match_neg = cur_support;

  }
}

/*****************************************************************
 * motif_index_com - compare the index of two motifs             *
 *****************************************************************/

int motif_index_com(motif *m1, motif *m2)
{
  //The last expression of any motif should always be a word
  //Since the range at the end is meaningless if there are no words after it
  //Return TRUE if the index of m1's last expression
  //is before the index of m2's first expression, otherwise return FALSE
  if (m1 == m2)
    return FALSE;
  else
  {
    // m1 != m2
    int m1_end_index = ((expression *)(m1->end))->begin_index;
    int m2_head_index = ((expression *)(m2->head))->begin_index;
    if (m1_end_index <= m2_head_index)
      return TRUE;
    else
      return FALSE;
  }
}

/*****************************************************************
 * motif_discovery - discover all motifs with enough supports    *
 *****************************************************************/

vector_t *
motif_discovery(dstring_t *ds_seed, vector_t *motif_list, vtree_t *vtree[], int control, parameters_t *params)
{
  //Discover all motifs that has enough supports

  int one_word_begin = 0; //Begin index of the word length 1
  int one_word_end = dev_vector_size(motif_list); //End index of the word length 1

  int con_begin = one_word_begin; //Begin index of the concatenation motif
  int con_end = one_word_end;     //End index of the concatenation motif

  int has_new_motif = TRUE;
  int is_first_one = TRUE;
  int cur_support;
  while (has_new_motif)
  {
    //while there are motifs that haven't been concatenated
    motif *new_motif;
    has_new_motif = FALSE; 
    is_first_one = TRUE;

    for (int i = con_begin; i < con_end; i++)
    {
      //From begin to the end of motifs made of same number of words
      for (int j = one_word_begin; j < one_word_end; j++)
      {
        //From the begin to the end of motif contain only one word
        motif *motif_i = (motif *)dev_vector_get(motif_list, i);
        motif *motif_j = (motif *)dev_vector_get(motif_list, j);

        // determine if j's index in the seed sequence is greater than i's last word's index
        if (motif_index_com(motif_i, motif_j))
        {

          new_motif = motif_concatenate(motif_i, motif_j);
          
          cur_support = 0;
          if (has_enough_support(ds_seed, vtree, new_motif, &cur_support, params) == TRUE)
          {
            new_motif->num_match_pos = cur_support;
            
            dev_vector_add(motif_list, new_motif);

            if (is_first_one)
            {
              has_new_motif = TRUE;
              con_begin = dev_vector_size(motif_list) - 1;
              is_first_one = FALSE;
            }
          }
          else
          {
            //Don't have enough support, dump the motif
            dev_free_motif(new_motif);
          }
        }
      }
    }
    con_end = dev_vector_size(motif_list); //while loop
  }
  return motif_list;
}

/*****************************************************************
 * motif_ranking - get a p-value for each motifs                 *
 *****************************************************************/

vector_t *
motif_ranking(vector_t *motif_list)
{
  double fisher_left_p, fisher_right_p, fisher_twosided_p;
  motif *m1;
  //Give each motif a p_value
  for (int i = 0; i < dev_vector_size(motif_list); i++)
  {
    m1 = dev_vector_get(motif_list, i);
    kt_fisher_exact(m1->num_match_pos, NUM_INPUT_STRING_POS - m1->num_match_pos,
                    m1->num_match_neg, NUM_INPUT_STRING_NEG - m1->num_match_neg,
                    &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
    m1->p_value = fisher_right_p;

    // printf("11 = %d, 12 = %d\n", m1->num_match_pos, NUM_INPUT_STRING_POS - m1->num_match_pos);
    // printf("21 = %d, 22 = %d\n", m1->num_match_neg, NUM_INPUT_STRING_NEG - m1->num_match_neg);
    // printf("p_val = %.4e \n", fisher_right_p);
  }
}


/*****************************************************************
 * save_motifs - output all motifs to a designated file          *
 *****************************************************************/

void save_motifs(char seed[], vector_t *motif_list, parameters_t *params)
{
  dev_log(1, "[ save_motifs ]");

  FILE *fh;

  fh = dev_fopen(params->match_file, "w");

  fprintf(fh, "<motifs>\n");

  for (int i = 0; i < dev_vector_size(motif_list); i++)
  {

    motif *m = (motif *)dev_vector_get(motif_list, i);

    fprintf(fh, "  <motif id=\"%d\" offset=\"%d\" pvalue=\"%.4e\">\n", i, ((expression *)(m->head))->begin_index, m->p_value);

    char *word;

    expression *e = m->head;
    int index = 0;
    char word_temp[params->word_length+1];

    while (e != NULL)
    {
      if (e->type == enum_word)
      {
        strncpy(word_temp, &seed[e->begin_index], e->length);
        word_temp[params->word_length] = '\0';
        fprintf(fh, "    <word>%s</word>\n", word_temp);
      }
      else if (e->type == enum_range)
      {
        fprintf(fh, "    <range length=\"%d\"></range>\n", e->length);
      }
      e = e->next;
    }
    fprintf(fh, "  </motif>\n");
  }

  fprintf(fh, "</motifs>\n");

  fclose(fh);
}


/*****************************************************************
 * main - main program                                           *
 *****************************************************************/

int main(int argc, char *argv[])
{
  int num_seqs;
  char **seqs_pos, **seqs_neg;
  char **descs_pos, **descs_neg;
  parameters_t params;
  dstring_t *ds, *ds_seed;

  if (argc == 1)
    display_usage_and_exit();

  /* initialisations */
  dev_init();
  param_init(&params);

  /* processing arguments */

  process_argv(argc, argv, &params);

  /* reading data */

  NUM_INPUT_STRING = NUM_INPUT_STRING_POS = bio_read_fasta(params.filename_pos, &seqs_pos, &descs_pos, isnuc);
  printf("NUM_INPUT_STRING = %d\n", NUM_INPUT_STRING);
  // assert( params.seed >= 0 && params.seed < num_seqs );

  params.min_support_num = NUM_INPUT_STRING * params.min_support;
  printf("min_support_num = %d\n", params.min_support_num);
  
  vtree_t **vtree ;
  vtree = (vtree_t **) malloc(NUM_INPUT_STRING_POS * sizeof(vtree_t *));
  // print description
  // for(int i=0; i < NUM_INPUT_STRING; i++) printf("%s\n", descs_pos[i]);
  ds_seed =  dev_digitalize(&uppercase, seqs_pos[0]);

  /*words enumerate - POS*/
  for (int i = 0; i < NUM_INPUT_STRING; i++)
  {
    // to keep consistency with Seed, we transform all
    // characters to upper case in here
    ds = dev_digitalize(&uppercase, seqs_pos[i]);
    vtree[i] = vtree_create(ds);
    dev_free_dstring( ds );
  }
  printf("------ Phase 1 : Word Enumeration------ \n");
  vector_t *motif_list = words_enumerate(seqs_pos[0], vtree, &params);
  printf("Vector_size_after_Enumerate = %d\n", dev_vector_size(motif_list));
  
  /*print_motif after words enumerate -- debug*/
  // for(int i=0; i<dev_vector_size(motif_list); i++){
  //   printf("-------PRINT - Motif : %d-------\n", i);
  //   print_motif(seqs_pos[0], (motif*)dev_vector_get(motif_list, i),  &params);
  // }

  /*motif discovery - POS*/
  printf("------ Phase 2 : Motif Discovery ------ \n");
  motif_list = motif_discovery(ds_seed, motif_list, vtree, FALSE, &params);
  printf("Vector_size_after_Discovery = %d\n", dev_vector_size(motif_list));

  /* free vtree list for positive input strings */
  printf("Free Memory for Vtree_list_Pos \n");
  dev_free_vtree_array((void **) vtree, NUM_INPUT_STRING, TRUE);

  /* motif discovery - NEG */
  NUM_INPUT_STRING = NUM_INPUT_STRING_NEG = bio_read_fasta( params.filename_neg, &seqs_neg, &descs_neg, isnuc );
  
  vtree = (vtree_t **) malloc(NUM_INPUT_STRING_NEG * sizeof(vtree_t *));

  for(int i=0; i<NUM_INPUT_STRING; i++){
  // read the negative data set and create vtree list
    ds = dev_digitalize(&uppercase, seqs_neg[i]);
    vtree[i] = vtree_create(ds);
    dev_free_dstring( ds );
  }

  /* travel through negative set and record support for each motifs */
  printf("------ Phase 3 : Get Support for Negative dataset ------ \n");
  
  get_neg_support(ds_seed, vtree, motif_list, &params);
  
  // /* free vtree list for negative input strings */
  printf("Free Memory for Vtree_list_Neg \n");
  dev_free_vtree_array((void **) vtree, NUM_INPUT_STRING, TRUE);

  // /*print_motif after words discovery -- debug*/
  // for(int i=0; i<dev_vector_size(motif_list); i++){
  //   printf("-------PRINT - Motif : %d-------\n", i);
  //   print_motif(seqs_pos[0], (motif*)dev_vector_get(motif_list, i),  &params);
  // }

  /*motif ranking*/
  printf("------ Phase 4 : motif_ranking ------ \n");
  motif_ranking(motif_list);

  /* Save motifs */
  if(params.match_file != NULL)
    {
      printf("------ Phase 5 : save_motifs ------ \n");
      save_motifs(seqs_pos[0], motif_list, &params);
    }
  else printf("match file is NULL, motifs will not be saved;\nIf you wish to save matches, please specify a match file path.\n");

  /* post-processings */
  printf("------ Clean up ------ \n");
  
  dev_free_vector( motif_list,  ( void ( * )( void * ) ) dev_free_motif );
  dev_free_dstring(ds_seed);
  
  exit(EXIT_SUCCESS);
}

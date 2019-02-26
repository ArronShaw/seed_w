/*                               -*- Mode: C -*- 
 * seed.h --- RNA secondary structure motif inference
 * Author          : Truong Nguyen and Marcel Turcotte
 * Created On      : Thu Jul  7 13:31:36 2005
 * Last Modified By: Marcel Turcotte
 * Last Modified On: Fri Aug  5 18:01:08 2005
 *
 * This copyrighted source code is freely distributed under the terms
 * of the GNU General Public License. 
 * See the files COPYRIGHT and LICENSE for details.
 */

// #ifndef SEED_H
// #define SEED_H

#include "libvtree.h"
#include <sys/types.h>

// typedef struct expression expression;

enum type
{
  enum_word = 0, 
  enum_range = 1
};

typedef struct expression_s{
	int type;// 0 for word, and 1 for range
	int begin_index;
	int length;
	struct expression_s *next;
}expression;

typedef struct {
	expression *head;
  expression *end;
	int words_count;
	int num_match_pos;
	int num_match_neg;
	double p_value;
}motif;

/*****************************************************************
 * alphabet - lower case letters alphabet                        *
 *****************************************************************/

static code_t codes_lower[] =
    {{'*', 0},
     {'a', 1},
     {'b', 2},
     {'c', 3},
     {'d', 4},
     {'e', 5},
     {'f', 6},
     {'g', 7},
     {'h', 8},
     {'i', 9},
     {'j', 10},
     {'k', 11},
     {'l', 12},
     {'m', 13},
     {'n', 14},
     {'o', 15},
     {'p', 16},
     {'q', 17},
     {'r', 18},
     {'s', 19},
     {'t', 20},
     {'u', 21},
     {'v', 22},
     {'w', 23},
     {'x', 24},
     {'y', 25},
     {'z', 26}};

static alphabet_t lowercase = {codes_lower, 27, 27};

static code_t codes_upper[] =
    {{'*', 0},
     {'A', 1},
     {'B', 2},
     {'C', 3},
     {'D', 4},
     {'E', 5},
     {'F', 6},
     {'G', 7},
     {'H', 8},
     {'I', 9},
     {'J', 10},
     {'K', 11},
     {'L', 12},
     {'M', 13},
     {'N', 14},
     {'O', 15},
     {'P', 16},
     {'Q', 17},
     {'R', 18},
     {'S', 19},
     {'T', 20},
     {'U', 21},
     {'V', 22},
     {'W', 23},
     {'X', 24},
     {'Y', 25},
     {'Z', 26}};

static alphabet_t uppercase = {codes_upper, 27, 27};

extern vector_t *dev_new_vector();

extern int match_edge(dstring_t *ds_seed, vtree_t *v, interval2_t *interval, expression *e, int pos, int offset);

extern int match_node(dstring_t *ds_seed, vtree_t *v, interval2_t *interval, expression *e, int pos, int offset);

extern int match_motif(dstring_t *ds_seed, interval2_t *i0, vtree_t *v, motif *m);

extern void print_expression(char *seed, expression *e);
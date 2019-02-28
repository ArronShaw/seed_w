/*
* 
* @Author: Weipeng Xu
* @Date:   2019-02-27 21:18:06
* @Last Modified by:   Weipeng Xu
* @Last Modified time: 2019-02-27 21:18:06
*
* This copyrighted source code is freely distributed under the terms
* of the GNU General Public License. 
* See the files COPYRIGHT and LICENSE for details.
*/

/*****************************************************************
 * param_t -                                                     *
 *****************************************************************/
typedef struct {
 int seed;
 float min_support;
 int min_support_num;
 int word_length;
 char *version;
 char *filename_pos;
 char *filename_neg;
 char *match_file;
 long match_count;
} parameters_t;

/*****************************************************************
 * Default values for the parameters                             *
 *****************************************************************/
#define DEFAULT_SEED 0
#define FILENAME_POS NULL
#define FILENAME_NEG NULL
#define MATCH_FILE NULL
#define MIN_SUPPORT 0.70
#define WORD_LENGTH 3
#define VERSION "SEED_W   v1.0"
#define VECTOR_INI_CAP 1000 //vector initial capacity
#define INCREMENT 1000 // increment

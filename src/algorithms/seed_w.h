// #include "libvtree.h"
// #include <sys/types.h>
// #define VERSION "1.0"



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
// extern vector_t *words_enumerate(char seed[], vtree_t *vtree[], int control)
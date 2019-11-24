/*
*
* @Author: Weipeng Xu
* @Date:   2019-02-27 21:28:30
* @Last Modified by:   Weipeng Xu
* @Last Modified time: 2019-02-27 21:28:42
*
* This copyrighted source code is freely distributed under the terms
* of the GNU General Public License.
* See the files COPYRIGHT and LICENSE for details.
*/
#include "libdev.h"
#include "seq.h"
// #include "stems.h"
#include "match_w.h"
#include "seed_w.h"
int match_motif(dstring_t *ds_seed, interval2_t *i0, vtree_t *v, motif *m,
                parameters_t *params, int mis_match_motif);
int match_node(dstring_t *ds_seed, vtree_t *v, interval2_t *interval,
               expression *e, int pos, int offset, parameters_t *params,
               int mis_match_motif);
int match_edge(dstring_t *ds_seed, vtree_t *v, interval2_t *interval,
               expression *e, int pos, int offset, parameters_t *params,
               int mis_match_motif);

int getCharfromWord(dstring_t *ds_seed, expression *e, int offset) {
  return ds_seed->text[e->begin_index + offset];
}

int error_case() { return -1; }

int match_node(dstring_t *ds_seed, vtree_t *v, interval2_t *interval,
               expression *e, int pos, int offset, parameters_t *params,
               int mis_match_motif) {

  int queryFound = FALSE;

  vector_t *childs;

  childs = vtree_getChildIntervals(v, interval);

  while (dev_vector_size(childs) > 0 && !queryFound) {

    // take out childs one by one
    interval2_t *child = (interval2_t *)dev_vector_serve(childs);

    if (match_edge(ds_seed, v, child, e, pos, offset, params, mis_match_motif))
      queryFound = TRUE;

    dev_free(child);
  }

  dev_free_vector(childs, free);

  return queryFound;
}

int match_edge(dstring_t *ds_seed, vtree_t *v, interval2_t *interval,
               expression *e, int pos, int offset, parameters_t *params,
               int mis_match_motif) {

  int result;

  /* end of expression -- the expression was found in the vtree*/
  if (e == NULL) {

    return TRUE;
  }

  /* at an internal node? */
  if (interval->i != interval->j &&
      pos == vtree_getlcp(v, interval->i, interval->j)) { // ???

    return match_node(ds_seed, v, interval, e, pos, offset, params,
                      mis_match_motif);
  }

  /* else, which is at an edge*/
  switch (e->type) {

  case enum_word:

    if (offset >= e->length) {
      return match_edge(ds_seed, v, interval, e->next, pos, 0, params,
                        mis_match_motif);
    }

    int a = v->text[v->suftab[interval->i] + pos];

    if (ister(a))
      return FALSE;

    int b = getCharfromWord(ds_seed, e, offset);

    // if(a != b) return FALSE;

    if (a != b) {
      if (mis_match_motif < params->mis_match_motif) {
        // return match_edge();
        result = match_edge(ds_seed, v, interval, e, pos + 1, offset + 1,
                            params, ++mis_match_motif);
        // printf("Mismatch = %d\n", params->mis_match_motif);
      } else {
        return FALSE;
      }
    }

    else {
      // a == b, continue to compare next character
      result = match_edge(ds_seed, v, interval, e, pos + 1, offset + 1, params,
                          mis_match_motif);
    }

    break;

  case enum_range:

    // if (offset >= e->length)
    //   result = match_edge(ds_seed, v, interval, e->next, pos, 0, params, 0);

    if (offset >= (e->length - params->range_delta)) {
      //
      if (offset < (e->length + params->range_delta)) {
        // If within delta range
        // 1. try match word
        result = match_edge(ds_seed, v, interval, e->next, pos, 0, params, 0);

        // 2. skip current character
        if (result != TRUE) {
          a = v->text[v->suftab[interval->i] + pos];
          if (ister(a))
            return FALSE;

          result = match_edge(ds_seed, v, interval, e, pos + 1, offset + 1,
                              params, mis_match_motif);
        }
      } else {
        // Outside range delta, match word
        result = match_edge(ds_seed, v, interval, e->next, pos, 0, params, 0);
      }
    }

    else {
      a = v->text[v->suftab[interval->i] + pos];
      if (ister(a))
        return FALSE;

      result = match_edge(ds_seed, v, interval, e, pos + 1, offset + 1, params,
                          mis_match_motif);
    }

    break;

  default:
    error_case();
  }
}

int match_motif(dstring_t *ds_seed, interval2_t *i0, vtree_t *v, motif *m,
                parameters_t *params, int mis_match_motif) {
  // expect TRUE/FALSE indeicate if a motif is exist in a vtree or not
  // mismatch are not allowed in this version 1.0, will be improved in 1.1
  // version

  return match_node(ds_seed, v, i0, m->head, 0, 0, params, mis_match_motif);
}

void print_expression(char *seed, expression *e) {
  printf("expression = ");
  while (e != NULL) {
    if (e->type == enum_word) {
      for (int i = e->begin_index; i < e->begin_index + e->length; i++)
        printf("%c", seed[i]);
    } else if (e->type == enum_range) {
      for (int i = e->begin_index; i < e->begin_index + e->length; i++)
        printf("*");
    } else {
      dev_die("internal error, expression type error");
    }
    e = e->next;
  }
  printf("\n");
}

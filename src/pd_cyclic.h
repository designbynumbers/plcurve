/*

  pd_cyclic.h : This is a minimal cyclic group
  implementation which allows one to select a cyclic
  group element.

*/

#ifndef __PD_CYCLIC_H_ 
#define __PD_CYCLIC_H_ 1

typedef struct pd_cyclic_struct {

  pd_idx_t n;     /* number of elts */
  pd_idx_t *map;  /* i -> map[i] */
   
} pd_cyclic_t;

void           *pd_new_cyclic(void *n); /* Pointer to pd_idx_t n */
void            pd_free_cyclic(void **cyclic);

char           *pd_print_cyclic(void *cyclic);
void           *pd_copy_cyclic(void *cyclic); /* Make a new-memory copy */

void            pd_increment_cyclic(void *cyclic);
unsigned int    pd_ncyclics(void *cyclic);

bool            pd_cyclic_ok(void *cyclic);
int             pd_cyclic_cmp(const void *cyclicA,const void *cyclicB);

bool            pd_cyclics_unique(unsigned int ncyclics,pd_cyclic_t **cyclic_buf);

extern pd_iterops_t cyclic_ops;

#endif

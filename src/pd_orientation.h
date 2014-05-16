/*

  pd_orientation.h : This is a minimal group implementation
  which allows one to iterate over a collection of +/- orientations,
  for instance to generate all possible crossing signs for a 
  pd_code.

*/

#ifndef __PD_ORIENTATION_H_ 
#define __PD_ORIENTATION_H_ 1

typedef struct pd_orientation_struct {

  pd_idx_t   n;  /* Number of objects to keep track of */
  pd_or_t  *or;    /* PD_POS_ORIENTATION or PD_NEG_ORIENTATION for each */
 
} pd_orientation_t;

void           *pd_new_orientation(void *n); /* Pointer to pd_idx_t n */
void            pd_free_orientation(void **orientation);

char           *pd_print_orientation(void *orientation);
void           *pd_copy_orientation(void *orientation); /* Make a new-memory copy */

void            pd_increment_orientation(void *orientation);
unsigned int    pd_norientations(void *orientation);

bool            pd_orientation_ok(void *orientation);
int             pd_orientation_cmp(const void *orientationA,const void *orientationB);

bool            pd_orientations_unique(unsigned int norientations,pd_orientation_t **orientation_buf);

extern pd_iterops_t orientation_ops;

#endif

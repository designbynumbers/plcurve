/*

  pd_orientation.h : This is a minimal group implementation
  of Z/2Z, which feeds into the multidx data type to allow you 
  to iterate, say, over all the possible crossing signs in a diagram.

*/

#ifndef __PD_ORIENTATION_H_ 
#define __PD_ORIENTATION_H_ 1

typedef struct pd_orientation_struct {

  pd_or_t  or;    /* PD_POS_ORIENTATION or PD_NEG_ORIENTATION */
 
} pd_orientation_t;

void           *pd_new_orientation(void *nodata); /* Initial data ignored */
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

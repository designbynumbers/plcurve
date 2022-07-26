/*

  pd_dihedral.h : This is a minimal dihedral group
  implementation which allows one to select a dihedral
  group element.

*/

#ifndef __PD_DIHEDRAL_H_ 
#define __PD_DIHEDRAL_H_ 1

typedef struct pd_dihedral_struct {

  pd_idx_t n;     /* number of elts */
  pd_idx_t *map;  /* i -> map[i] */
  
  pd_or_t  orient;    /* rotation or reflection */
 
} pd_dihedral_t;

void           *pd_new_dihedral(void *n); /* Pointer to pd_idx_t n */
void            pd_free_dihedral(void **dihedral);

char           *pd_print_dihedral(void *dihedral);
void           *pd_copy_dihedral(void *dihedral); /* Make a new-memory copy */

void            pd_increment_dihedral(void *dihedral);
unsigned int    pd_ndihedrals(void *dihedral);

bool            pd_dihedral_ok(void *dihedral);
int             pd_dihedral_cmp(const void *dihedralA,const void *dihedralB);

bool            pd_dihedrals_unique(unsigned int ndihedrals,pd_dihedral_t **dihedral_buf);

extern pd_iterops_t dihedral_ops;

#endif

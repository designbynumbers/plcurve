/*

  pd_perm.h:  This is a minimal permutation group implementation
  which just allows one to loop over permutations. 

*/

#ifndef __PD_PERM_H_ 
#define __PD_PERM_H_ 1

typedef struct pd_perm_struct {

  pd_idx_t n;
  pd_idx_t *map;

  pd_idx_t pc_idx; /* Precomputed index */

} pd_perm_t;

void        *pd_new_perm(void *n);
void         pd_free_perm(void **perm);

char        *pd_print_perm(void *perm);
void        *pd_copy_perm(void *perm); /* Make a new-memory copy of perm */

void         pd_increment_perm(void *perm);
unsigned int pd_nperms(void *perm); /* Counts number of possible vals for this perm */

bool         pd_perm_ok(void *perm);

bool         pd_perms_eq(pd_perm_t *permA,pd_perm_t *permB);
int          pd_perm_cmp(const void *permA,const void *permB);

bool         pd_perm_is_e(pd_perm_t *perm); /* Check for identity perm (for testing) */
pd_perm_t   *pd_compose_perms(pd_perm_t *A,pd_perm_t *B); /* product permutation (A * B)(x) = A(B(x)). */
void         pd_stareq_perm(pd_perm_t *A,pd_perm_t *B);  /* A *= B (updates A in-place) */
unsigned int pd_perm_period(pd_perm_t *A); /* computes period of A in permutation group. */

bool         pd_perms_unique(unsigned int nperms,pd_perm_t **perm_buf);

void         pd_regenerate_pcidx(pd_perm_t *perm); 
             /* Recomputes the index of this perm's map in precomputed data.*/

bool         pd_perm_pcdata_ok(bool verbose); 
             /* Self-tests on precomputed data. Prints output if verbose == true. */

extern pd_iterops_t perm_ops;

#endif

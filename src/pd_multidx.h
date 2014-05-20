/* 

   pd_multidx.h : Part of the COLD project. Cantarella/2012.

*/

/* We start with a general iterator called a multi-index.
   The idea of a multi-index is basically object-oriented:
   it's an array of "objects" which have an "increment"
   function and an "nvals" function.

   We can increment these successively to count all
   n-tuples of the objects. Internally, we maintain a
   counter i[j] which records which state each of the objects
   is in 
   
        0 <= i[j] < nvals[j]                            

  The number of possible values for a multi-index can be
  counted with pd_multidx_nvals, and the values can be
  accessed in order with pd_increment_multidx(idx) which
  operates in place on idx.

  The objects and methods can be NULL. To use a multidx, we write 
  code such as:

  objdata = calloc(nobj,sizeof(object_data_t));
  int i;

  for(i=0;i<nobj;i++) { 

     (initialize objdata[i] in here)

  }

  pd_multidx_t *idx = pd_new_multidx(nobj,(void **)objdata,obj_ops);
  int nstates = pd_multidx_nvals(idx);

  for(i=0;i<nstates;i++,pd_increment_multidx(idx)) { 

     (do stuff depending on idx -- we can read the data from 
      the object in each position with something like 

     int j;
     for(j=0;j<nobj;j++) { 

         (object_t *)(idx->obj[j])->object_field

     }

  }

  pd_free_multidx(&idx);
  
*/

#ifndef __PD_MULTIDX_H_
#define __PD_MULTIDX_H_ 1

/* First, we define a set of methods for an iterator.
   Basically, an iterator is something which can be
   operated on according to these rules in a meaningful
   way. 

   All of our (group-valued) iterators elsewhere in the 
   program will obey this interface. */

typedef struct pd_iterops_struct {

  void *(*new)(void *);   /* Given a pointer to initial data, create a new object */
  void  (*free)(void **); /* Free everything allocated by new_obj, set pointer to NULL. */

  char *(*print)(void *); /* Allocates a string representing obj. */
  void *(*copy)(void *);  /* Allocates a new-memory copy of obj. */

  void  (*incr)(void *);  /* Increments (in place) an object */
  unsigned int (*nvals)(void *); /* After this many calls to incr, obj should return to initial state. */

  bool  (*ok)(void *);    /* Performs internal consistency checks on obj. */
  int   (*cmp)(const void *,const void *); /* Compares POINTERS to obj for searching and */
                                           /* sorting arrays of POINTERS to obj types. */
} pd_iterops_t;

bool pd_iterops_eq(pd_iterops_t opsA,pd_iterops_t opsB);
  
typedef struct pd_multidx_struct {

  pd_idx_t       nobj;
  void         **obj;               /* Array of nobj "object" pointers which can be freed */

  unsigned int  *i;                 /* Records state of each object. */
  unsigned int  *nvals;             /* Records number of states for each object. */

  pd_iterops_t   ops;               /* Pointers to functions for "object" operations */

} pd_multidx_t;

pd_multidx_t *pd_new_multidx(pd_idx_t nobj,void **initialdata,pd_iterops_t ops);

/* We are expected to provide an (allocated) buffer to initial data for objects
   which will be created by the function ops->new. */
 
void          pd_free_multidx(pd_multidx_t **idx);
void          pd_increment_multidx(pd_multidx_t *idx);
unsigned int  pd_multidx_nvals(pd_multidx_t *idx);
bool          pd_multidx_ok(pd_multidx_t *idx);

int           pd_multidx_cmp(const void *idxAp, const void *idxBp);
/* A comparison function for (pd_multidx_t *) types, so idxAp and idxBp are (pd_multidx_t **) */

bool          pd_multidxs_unique(unsigned int nmultidxs,pd_multidx_t **multidx_buf);
/* Uses the comparison function pd_multidx_cmp to check that all buffer elts are unique
   in a buffer of POINTERS to pd_multidx_t. */

pd_multidx_t *pd_copy_multidx(pd_multidx_t *idx);
/* Makes a new-memory copy of idx. */

#endif

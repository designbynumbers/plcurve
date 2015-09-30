void pd_compacting_copy(void *source, size_t obj_size, size_t nobj,
			pd_idx_t ndeletions, pd_idx_t *deletions,
			void **target_var,
			pd_idx_t **target_idx,
			pd_idx_t *ntarget) /* Returns the size of the target buffer. */

/* This primitive copies (in order) all the elements of source_buf EXCEPT those
   in the buffer "deletions" to target buf. We allocate "target_var" ourselves,
   and return in the target_idx a buffer of size nobj so that

   target_idx[index in source] = index of copy, if the object was copied
                               = PD_UNSET_IDX, if the object was deleted

   The basic idea is that we're going to sort the buffer of deletions, and
   then read through while copying to the target buffer, incrementing the
   deletions buffer when we find one of the elements to delete.

   This is linear time, and most importantly, we only need to write (and debug)
   it once. If we delete EVERYTHING (this can happen), then we set target_var to
   NULL, ntarget to zero, and all the entries of target_idx to PD_UNSET_IDX.

*/

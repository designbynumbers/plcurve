/* 

   pd_sortedbuf.h : Code to maintain and search sorted buffers of pd_idx_t.


*/

#ifndef __PD_SORTEDBUF_H__
#define __PD_SORTEDBUF_H__ 1

bool      pd_sortedbuf_contains(pd_idx_t *buf,pd_idx_t size,pd_idx_t i);
/* Returns true if the buffer contains the number "i". */

pd_idx_t *pd_sortedbuf_insert(pd_idx_t *buf,pd_idx_t *size,pd_idx_t i);
/* Inserts "i" into the buffer, maintaining sort order, and making the
   buffer bigger if needed (note that "size" is passed as a pointer so
   that it can be updated. pd_sortedbuf_insert will free the original
   buffer and copy it to new memory. */

bool      pd_buf_contains(pd_idx_t *buf,pd_idx_t size,pd_idx_t i,pd_idx_t *where);
/* Returns true if the buffer contains the number "i" and sets "where" to the index. */

#endif

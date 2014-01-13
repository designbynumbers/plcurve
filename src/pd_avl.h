/* Produced by texiweb from libavl.w. */

/* libavl - library for manipulation of binary trees.
   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2004 Free Software
   Foundation, Inc.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA.
*/

#ifndef PDAVL_H
#define PDAVL_H 1

#include <stddef.h>

/* Function types. */
typedef int pdavl_comparison_func (const void *pdavl_a, const void *pdavl_b,
                                 void *pdavl_param);
typedef void pdavl_item_func (void *pdavl_item, void *pdavl_param);
typedef void *pdavl_copy_func (void *pdavl_item, void *pdavl_param);

#ifndef LIBPDAVL_ALLOCATOR
#define LIBPDAVL_ALLOCATOR
/* Memory allocator. */
struct libpdavl_allocator
  {
    void *(*libpdavl_malloc) (struct libpdavl_allocator *, size_t libpdavl_size);
    void (*libpdavl_free) (struct libpdavl_allocator *, void *libpdavl_block);
  };
#endif

/* Default memory allocator. */
extern struct libpdavl_allocator pdavl_allocator_default;
void *pdavl_malloc (struct libpdavl_allocator *, size_t);
void pdavl_free (struct libpdavl_allocator *, void *);

/* Maximum PDAVL tree height. */
#ifndef PDAVL_MAX_HEIGHT
#define PDAVL_MAX_HEIGHT 92
#endif

/* Tree data structure. */
struct pdavl_table
  {
    struct pdavl_node *pdavl_root;          /* Tree's root. */
    pdavl_comparison_func *pdavl_compare;   /* Comparison function. */
    void *pdavl_param;                    /* Extra argument to |pdavl_compare|. */
    struct libpdavl_allocator *pdavl_alloc; /* Memory allocator. */
    size_t pdavl_count;                   /* Number of items in tree. */
    unsigned long pdavl_generation;       /* Generation number. */
  };

/* An PDAVL tree node. */
struct pdavl_node
  {
    struct pdavl_node *pdavl_link[2];  /* Subtrees. */
    void *pdavl_data;                /* Pointer to data. */
    signed char pdavl_balance;       /* Balance factor. */
  };

/* PDAVL traverser structure. */
struct pdavl_traverser
  {
    struct pdavl_table *pdavl_table;        /* Tree being traversed. */
    struct pdavl_node *pdavl_node;          /* Current node in tree. */
    struct pdavl_node *pdavl_stack[PDAVL_MAX_HEIGHT];
                                        /* All the nodes above |pdavl_node|. */
    size_t pdavl_height;                  /* Number of nodes in |pdavl_parent|. */
    unsigned long pdavl_generation;       /* Generation number. */
  };

/* Table functions. */
struct pdavl_table *pdavl_create (pdavl_comparison_func *, void *,
                              struct libpdavl_allocator *);
struct pdavl_table *pdavl_copy (const struct pdavl_table *, pdavl_copy_func *,
                            pdavl_item_func *, struct libpdavl_allocator *);
void pdavl_destroy (struct pdavl_table *, pdavl_item_func *);
void **pdavl_probe (struct pdavl_table *, void *);
void *pdavl_insert (struct pdavl_table *, void *);
void *pdavl_replace (struct pdavl_table *, void *);
void *pdavl_delete (struct pdavl_table *, const void *);
void *pdavl_find (const struct pdavl_table *, const void *);
void pdavl_assert_insert (struct pdavl_table *, void *);
void *pdavl_assert_delete (struct pdavl_table *, void *);

#define pdavl_count(table) ((size_t) (table)->pdavl_count)

/* Table traverser functions. */
void pdavl_t_init (struct pdavl_traverser *, struct pdavl_table *);
void *pdavl_t_first (struct pdavl_traverser *, struct pdavl_table *);
void *pdavl_t_last (struct pdavl_traverser *, struct pdavl_table *);
void *pdavl_t_find (struct pdavl_traverser *, struct pdavl_table *, void *);
void *pdavl_t_insert (struct pdavl_traverser *, struct pdavl_table *, void *);
void *pdavl_t_copy (struct pdavl_traverser *, const struct pdavl_traverser *);
void *pdavl_t_next (struct pdavl_traverser *);
void *pdavl_t_prev (struct pdavl_traverser *);
void *pdavl_t_cur (struct pdavl_traverser *);
void *pdavl_t_replace (struct pdavl_traverser *, void *);

#endif /* pdavl.h */

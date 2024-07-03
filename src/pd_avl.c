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

#ifdef HAVE_CONFIG_H
  #include<config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pd_avl.h"

/* Creates and returns a new table
   with comparison function |compare| using parameter |param|
   and memory allocator |allocator|.
   Returns |NULL| if memory allocation failed. */
struct pdavl_table *
pdavl_create (pdavl_comparison_func *compare, void *param,
            struct libpdavl_allocator *allocator)
{
  struct pdavl_table *tree;

  assert (compare != NULL);

  if (allocator == NULL)
    allocator = &pdavl_allocator_default;

  tree = allocator->libpdavl_malloc (allocator, sizeof *tree);
  if (tree == NULL)
    return NULL;

  tree->pdavl_root = NULL;
  tree->pdavl_compare = compare;
  tree->pdavl_param = param;
  tree->pdavl_alloc = allocator;
  tree->pdavl_count = 0;
  tree->pdavl_generation = 0;

  return tree;
}

/* Search |tree| for an item matching |item|, and return it if found.
   Otherwise return |NULL|. */
void *
pdavl_find (const struct pdavl_table *tree, const void *item)
{
  const struct pdavl_node *p;

  assert (tree != NULL && item != NULL);
  for (p = tree->pdavl_root; p != NULL; )
    {
      int cmp = tree->pdavl_compare (item, p->pdavl_data, tree->pdavl_param);

      if (cmp < 0)
        p = p->pdavl_link[0];
      else if (cmp > 0)
        p = p->pdavl_link[1];
      else /* |cmp == 0| */
        return p->pdavl_data;
    }

  return NULL;
}

/* Inserts |item| into |tree| and returns a pointer to |item|'s address.
   If a duplicate item is found in the tree,
   returns a pointer to the duplicate without inserting |item|.
   Returns |NULL| in case of memory allocation failure. */
void **
pdavl_probe (struct pdavl_table *tree, void *item)
{
  struct pdavl_node *y, *z; /* Top node to update balance factor, and parent. */
  struct pdavl_node *p, *q; /* Iterator, and parent. */
  struct pdavl_node *n;     /* Newly inserted node. */
  struct pdavl_node *w;     /* New root of rebalanced subtree. */
  int dir;                /* Direction to descend. */

  unsigned char da[PDAVL_MAX_HEIGHT]; /* Cached comparison results. */
  int k = 0;              /* Number of cached results. */

  assert (tree != NULL && item != NULL);

  z = (struct pdavl_node *) &tree->pdavl_root;
  y = tree->pdavl_root;
  dir = 0;
  for (q = z, p = y; p != NULL; q = p, p = p->pdavl_link[dir])
    {
      int cmp = tree->pdavl_compare (item, p->pdavl_data, tree->pdavl_param);
      if (cmp == 0)
        return &p->pdavl_data;

      if (p->pdavl_balance != 0)
        z = q, y = p, k = 0;
      da[k++] = dir = cmp > 0;
    }

  n = q->pdavl_link[dir] =
    tree->pdavl_alloc->libpdavl_malloc (tree->pdavl_alloc, sizeof *n);
  if (n == NULL)
    return NULL;

  tree->pdavl_count++;
  n->pdavl_data = item;
  n->pdavl_link[0] = n->pdavl_link[1] = NULL;
  n->pdavl_balance = 0;
  if (y == NULL)
    return &n->pdavl_data;

  for (p = y, k = 0; p != n; p = p->pdavl_link[da[k]], k++)
    if (da[k] == 0)
      p->pdavl_balance--;
    else
      p->pdavl_balance++;

  if (y->pdavl_balance == -2)
    {
      struct pdavl_node *x = y->pdavl_link[0];
      if (x->pdavl_balance == -1)
        {
          w = x;
          y->pdavl_link[0] = x->pdavl_link[1];
          x->pdavl_link[1] = y;
          x->pdavl_balance = y->pdavl_balance = 0;
        }
      else
        {
          assert (x->pdavl_balance == +1);
          w = x->pdavl_link[1];
          x->pdavl_link[1] = w->pdavl_link[0];
          w->pdavl_link[0] = x;
          y->pdavl_link[0] = w->pdavl_link[1];
          w->pdavl_link[1] = y;
          if (w->pdavl_balance == -1)
            x->pdavl_balance = 0, y->pdavl_balance = +1;
          else if (w->pdavl_balance == 0)
            x->pdavl_balance = y->pdavl_balance = 0;
          else /* |w->pdavl_balance == +1| */
            x->pdavl_balance = -1, y->pdavl_balance = 0;
          w->pdavl_balance = 0;
        }
    }
  else if (y->pdavl_balance == +2)
    {
      struct pdavl_node *x = y->pdavl_link[1];
      if (x->pdavl_balance == +1)
        {
          w = x;
          y->pdavl_link[1] = x->pdavl_link[0];
          x->pdavl_link[0] = y;
          x->pdavl_balance = y->pdavl_balance = 0;
        }
      else
        {
          assert (x->pdavl_balance == -1);
          w = x->pdavl_link[0];
          x->pdavl_link[0] = w->pdavl_link[1];
          w->pdavl_link[1] = x;
          y->pdavl_link[1] = w->pdavl_link[0];
          w->pdavl_link[0] = y;
          if (w->pdavl_balance == +1)
            x->pdavl_balance = 0, y->pdavl_balance = -1;
          else if (w->pdavl_balance == 0)
            x->pdavl_balance = y->pdavl_balance = 0;
          else /* |w->pdavl_balance == -1| */
            x->pdavl_balance = +1, y->pdavl_balance = 0;
          w->pdavl_balance = 0;
        }
    }
  else
    return &n->pdavl_data;
  z->pdavl_link[y != z->pdavl_link[0]] = w;

  tree->pdavl_generation++;
  return &n->pdavl_data;
}

/* Inserts |item| into |table|.
   Returns |NULL| if |item| was successfully inserted
   or if a memory allocation error occurred.
   Otherwise, returns the duplicate item. */
void *
pdavl_insert (struct pdavl_table *table, void *item)
{
  void **p = pdavl_probe (table, item);
  return p == NULL || *p == item ? NULL : *p;
}

/* Inserts |item| into |table|, replacing any duplicate item.
   Returns |NULL| if |item| was inserted without replacing a duplicate,
   or if a memory allocation error occurred.
   Otherwise, returns the item that was replaced. */
void *
pdavl_replace (struct pdavl_table *table, void *item)
{
  void **p = pdavl_probe (table, item);
  if (p == NULL || *p == item)
    return NULL;
  else
    {
      void *r = *p;
      *p = item;
      return r;
    }
}

/* Deletes from |tree| and returns an item matching |item|.
   Returns a null pointer if no matching item found. */
void *
pdavl_delete (struct pdavl_table *tree, const void *item)
{
  /* Stack of nodes. */
  struct pdavl_node *pa[PDAVL_MAX_HEIGHT]; /* Nodes. */
  unsigned char da[PDAVL_MAX_HEIGHT];    /* |pdavl_link[]| indexes. */
  int k;                               /* Stack pointer. */

  struct pdavl_node *p;   /* Traverses tree to find node to delete. */
  int cmp;              /* Result of comparison between |item| and |p|. */

  assert (tree != NULL && item != NULL);

  k = 0;
  p = (struct pdavl_node *) &tree->pdavl_root;
  for (cmp = -1; cmp != 0;
       cmp = tree->pdavl_compare (item, p->pdavl_data, tree->pdavl_param))
    {
      int dir = cmp > 0;

      pa[k] = p;
      da[k++] = dir;

      p = p->pdavl_link[dir];
      if (p == NULL)
        return NULL;
    }
  item = p->pdavl_data;

  if (p->pdavl_link[1] == NULL)
    pa[k - 1]->pdavl_link[da[k - 1]] = p->pdavl_link[0];
  else
    {
      struct pdavl_node *r = p->pdavl_link[1];
      if (r->pdavl_link[0] == NULL)
        {
          r->pdavl_link[0] = p->pdavl_link[0];
          r->pdavl_balance = p->pdavl_balance;
          pa[k - 1]->pdavl_link[da[k - 1]] = r;
          da[k] = 1;
          pa[k++] = r;
        }
      else
        {
          struct pdavl_node *s;
          int j = k++;

          for (;;)
            {
              da[k] = 0;
              pa[k++] = r;
              s = r->pdavl_link[0];
              if (s->pdavl_link[0] == NULL)
                break;

              r = s;
            }

          s->pdavl_link[0] = p->pdavl_link[0];
          r->pdavl_link[0] = s->pdavl_link[1];
          s->pdavl_link[1] = p->pdavl_link[1];
          s->pdavl_balance = p->pdavl_balance;

          pa[j - 1]->pdavl_link[da[j - 1]] = s;
          da[j] = 1;
          pa[j] = s;
        }
    }

  tree->pdavl_alloc->libpdavl_free (tree->pdavl_alloc, p);

  assert (k > 0);
  while (--k > 0)
    {
      struct pdavl_node *y = pa[k];

      if (da[k] == 0)
        {
          y->pdavl_balance++;
          if (y->pdavl_balance == +1)
            break;
          else if (y->pdavl_balance == +2)
            {
              struct pdavl_node *x = y->pdavl_link[1];
              if (x->pdavl_balance == -1)
                {
                  struct pdavl_node *w;
                  assert (x->pdavl_balance == -1);
                  w = x->pdavl_link[0];
                  x->pdavl_link[0] = w->pdavl_link[1];
                  w->pdavl_link[1] = x;
                  y->pdavl_link[1] = w->pdavl_link[0];
                  w->pdavl_link[0] = y;
                  if (w->pdavl_balance == +1)
                    x->pdavl_balance = 0, y->pdavl_balance = -1;
                  else if (w->pdavl_balance == 0)
                    x->pdavl_balance = y->pdavl_balance = 0;
                  else /* |w->pdavl_balance == -1| */
                    x->pdavl_balance = +1, y->pdavl_balance = 0;
                  w->pdavl_balance = 0;
                  pa[k - 1]->pdavl_link[da[k - 1]] = w;
                }
              else
                {
                  y->pdavl_link[1] = x->pdavl_link[0];
                  x->pdavl_link[0] = y;
                  pa[k - 1]->pdavl_link[da[k - 1]] = x;
                  if (x->pdavl_balance == 0)
                    {
                      x->pdavl_balance = -1;
                      y->pdavl_balance = +1;
                      break;
                    }
                  else
                    x->pdavl_balance = y->pdavl_balance = 0;
                }
            }
        }
      else
        {
          y->pdavl_balance--;
          if (y->pdavl_balance == -1)
            break;
          else if (y->pdavl_balance == -2)
            {
              struct pdavl_node *x = y->pdavl_link[0];
              if (x->pdavl_balance == +1)
                {
                  struct pdavl_node *w;
                  assert (x->pdavl_balance == +1);
                  w = x->pdavl_link[1];
                  x->pdavl_link[1] = w->pdavl_link[0];
                  w->pdavl_link[0] = x;
                  y->pdavl_link[0] = w->pdavl_link[1];
                  w->pdavl_link[1] = y;
                  if (w->pdavl_balance == -1)
                    x->pdavl_balance = 0, y->pdavl_balance = +1;
                  else if (w->pdavl_balance == 0)
                    x->pdavl_balance = y->pdavl_balance = 0;
                  else /* |w->pdavl_balance == +1| */
                    x->pdavl_balance = -1, y->pdavl_balance = 0;
                  w->pdavl_balance = 0;
                  pa[k - 1]->pdavl_link[da[k - 1]] = w;
                }
              else
                {
                  y->pdavl_link[0] = x->pdavl_link[1];
                  x->pdavl_link[1] = y;
                  pa[k - 1]->pdavl_link[da[k - 1]] = x;
                  if (x->pdavl_balance == 0)
                    {
                      x->pdavl_balance = +1;
                      y->pdavl_balance = -1;
                      break;
                    }
                  else
                    x->pdavl_balance = y->pdavl_balance = 0;
                }
            }
        }
    }

  tree->pdavl_count--;
  tree->pdavl_generation++;
  return (void *) item;
}

/* Refreshes the stack of parent pointers in |trav|
   and updates its generation number. */
static void
trav_refresh (struct pdavl_traverser *trav)
{
  assert (trav != NULL);

  trav->pdavl_generation = trav->pdavl_table->pdavl_generation;

  if (trav->pdavl_node != NULL)
    {
      pdavl_comparison_func *cmp = trav->pdavl_table->pdavl_compare;
      void *param = trav->pdavl_table->pdavl_param;
      struct pdavl_node *node = trav->pdavl_node;
      struct pdavl_node *i;

      trav->pdavl_height = 0;
      for (i = trav->pdavl_table->pdavl_root; i != node; )
        {
          assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
          assert (i != NULL);

          trav->pdavl_stack[trav->pdavl_height++] = i;
          i = i->pdavl_link[cmp (node->pdavl_data, i->pdavl_data, param) > 0];
        }
    }
}

/* Initializes |trav| for use with |tree|
   and selects the null node. */
void
pdavl_t_init (struct pdavl_traverser *trav, struct pdavl_table *tree)
{
  trav->pdavl_table = tree;
  trav->pdavl_node = NULL;
  trav->pdavl_height = 0;
  trav->pdavl_generation = tree->pdavl_generation;
}

/* Initializes |trav| for |tree|
   and selects and returns a pointer to its least-valued item.
   Returns |NULL| if |tree| contains no nodes. */
void *
pdavl_t_first (struct pdavl_traverser *trav, struct pdavl_table *tree)
{
  struct pdavl_node *x;

  assert (tree != NULL && trav != NULL);

  trav->pdavl_table = tree;
  trav->pdavl_height = 0;
  trav->pdavl_generation = tree->pdavl_generation;

  x = tree->pdavl_root;
  if (x != NULL)
    while (x->pdavl_link[0] != NULL)
      {
        assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
        trav->pdavl_stack[trav->pdavl_height++] = x;
        x = x->pdavl_link[0];
      }
  trav->pdavl_node = x;

  return x != NULL ? x->pdavl_data : NULL;
}

/* Initializes |trav| for |tree|
   and selects and returns a pointer to its greatest-valued item.
   Returns |NULL| if |tree| contains no nodes. */
void *
pdavl_t_last (struct pdavl_traverser *trav, struct pdavl_table *tree)
{
  struct pdavl_node *x;

  assert (tree != NULL && trav != NULL);

  trav->pdavl_table = tree;
  trav->pdavl_height = 0;
  trav->pdavl_generation = tree->pdavl_generation;

  x = tree->pdavl_root;
  if (x != NULL)
    while (x->pdavl_link[1] != NULL)
      {
        assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
        trav->pdavl_stack[trav->pdavl_height++] = x;
        x = x->pdavl_link[1];
      }
  trav->pdavl_node = x;

  return x != NULL ? x->pdavl_data : NULL;
}

/* Searches for |item| in |tree|.
   If found, initializes |trav| to the item found and returns the item
   as well.
   If there is no matching item, initializes |trav| to the null item
   and returns |NULL|. */
void *
pdavl_t_find (struct pdavl_traverser *trav, struct pdavl_table *tree, void *item)
{
  struct pdavl_node *p, *q;

  assert (trav != NULL && tree != NULL && item != NULL);
  trav->pdavl_table = tree;
  trav->pdavl_height = 0;
  trav->pdavl_generation = tree->pdavl_generation;
  for (p = tree->pdavl_root; p != NULL; p = q)
    {
      int cmp = tree->pdavl_compare (item, p->pdavl_data, tree->pdavl_param);

      if (cmp < 0)
        q = p->pdavl_link[0];
      else if (cmp > 0)
        q = p->pdavl_link[1];
      else /* |cmp == 0| */
        {
          trav->pdavl_node = p;
          return p->pdavl_data;
        }

      assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
      trav->pdavl_stack[trav->pdavl_height++] = p;
    }

  trav->pdavl_height = 0;
  trav->pdavl_node = NULL;
  return NULL;
}

/* Attempts to insert |item| into |tree|.
   If |item| is inserted successfully, it is returned and |trav| is
   initialized to its location.
   If a duplicate is found, it is returned and |trav| is initialized to
   its location.  No replacement of the item occurs.
   If a memory allocation failure occurs, |NULL| is returned and |trav|
   is initialized to the null item. */
void *
pdavl_t_insert (struct pdavl_traverser *trav, struct pdavl_table *tree, void *item)
{
  void **p;

  assert (trav != NULL && tree != NULL && item != NULL);

  p = pdavl_probe (tree, item);
  if (p != NULL)
    {
      trav->pdavl_table = tree;
      trav->pdavl_node =
        ((struct pdavl_node *)
         ((char *) p - offsetof (struct pdavl_node, pdavl_data)));
      trav->pdavl_generation = tree->pdavl_generation - 1;
      return *p;
    }
  else
    {
      pdavl_t_init (trav, tree);
      return NULL;
    }
}

/* Initializes |trav| to have the same current node as |src|. */
void *
pdavl_t_copy (struct pdavl_traverser *trav, const struct pdavl_traverser *src)
{
  assert (trav != NULL && src != NULL);

  if (trav != src)
    {
      trav->pdavl_table = src->pdavl_table;
      trav->pdavl_node = src->pdavl_node;
      trav->pdavl_generation = src->pdavl_generation;
      if (trav->pdavl_generation == trav->pdavl_table->pdavl_generation)
        {
          trav->pdavl_height = src->pdavl_height;
          memcpy (trav->pdavl_stack, (const void *) src->pdavl_stack,
                  sizeof *trav->pdavl_stack * trav->pdavl_height);
        }
    }

  return trav->pdavl_node != NULL ? trav->pdavl_node->pdavl_data : NULL;
}

/* Returns the next data item in inorder
   within the tree being traversed with |trav|,
   or if there are no more data items returns |NULL|. */
void *
pdavl_t_next (struct pdavl_traverser *trav)
{
  struct pdavl_node *x;

  assert (trav != NULL);

  if (trav->pdavl_generation != trav->pdavl_table->pdavl_generation)
    trav_refresh (trav);

  x = trav->pdavl_node;
  if (x == NULL)
    {
      return pdavl_t_first (trav, trav->pdavl_table);
    }
  else if (x->pdavl_link[1] != NULL)
    {
      assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
      trav->pdavl_stack[trav->pdavl_height++] = x;
      x = x->pdavl_link[1];

      while (x->pdavl_link[0] != NULL)
        {
          assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
          trav->pdavl_stack[trav->pdavl_height++] = x;
          x = x->pdavl_link[0];
        }
    }
  else
    {
      struct pdavl_node *y;

      do
        {
          if (trav->pdavl_height == 0)
            {
              trav->pdavl_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->pdavl_stack[--trav->pdavl_height];
        }
      while (y == x->pdavl_link[1]);
    }
  trav->pdavl_node = x;

  return x->pdavl_data;
}

/* Returns the previous data item in inorder
   within the tree being traversed with |trav|,
   or if there are no more data items returns |NULL|. */
void *
pdavl_t_prev (struct pdavl_traverser *trav)
{
  struct pdavl_node *x;

  assert (trav != NULL);

  if (trav->pdavl_generation != trav->pdavl_table->pdavl_generation)
    trav_refresh (trav);

  x = trav->pdavl_node;
  if (x == NULL)
    {
      return pdavl_t_last (trav, trav->pdavl_table);
    }
  else if (x->pdavl_link[0] != NULL)
    {
      assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
      trav->pdavl_stack[trav->pdavl_height++] = x;
      x = x->pdavl_link[0];

      while (x->pdavl_link[1] != NULL)
        {
          assert (trav->pdavl_height < PDAVL_MAX_HEIGHT);
          trav->pdavl_stack[trav->pdavl_height++] = x;
          x = x->pdavl_link[1];
        }
    }
  else
    {
      struct pdavl_node *y;

      do
        {
          if (trav->pdavl_height == 0)
            {
              trav->pdavl_node = NULL;
              return NULL;
            }

          y = x;
          x = trav->pdavl_stack[--trav->pdavl_height];
        }
      while (y == x->pdavl_link[0]);
    }
  trav->pdavl_node = x;

  return x->pdavl_data;
}

/* Returns |trav|'s current item. */
void *
pdavl_t_cur (struct pdavl_traverser *trav)
{
  assert (trav != NULL);

  return trav->pdavl_node != NULL ? trav->pdavl_node->pdavl_data : NULL;
}

/* Replaces the current item in |trav| by |new| and returns the item replaced.
   |trav| must not have the null item selected.
   The new item must not upset the ordering of the tree. */
void *
pdavl_t_replace (struct pdavl_traverser *trav, void *new)
{
  void *old;

  assert (trav != NULL && trav->pdavl_node != NULL && new != NULL);
  old = trav->pdavl_node->pdavl_data;
  trav->pdavl_node->pdavl_data = new;
  return old;
}

/* Destroys |new| with |pdavl_destroy (new, destroy)|,
   first setting right links of nodes in |stack| within |new|
   to null pointers to avoid touching uninitialized data. */
static void
copy_error_recovery (struct pdavl_node **stack, int height,
                     struct pdavl_table *new, pdavl_item_func *destroy)
{
  assert (stack != NULL && height >= 0 && new != NULL);

  for (; height > 2; height -= 2)
    stack[height - 1]->pdavl_link[1] = NULL;
  pdavl_destroy (new, destroy);
}

/* Copies |org| to a newly created tree, which is returned.
   If |copy != NULL|, each data item in |org| is first passed to |copy|,
   and the return values are inserted into the tree,
   with |NULL| return values taken as indications of failure.
   On failure, destroys the partially created new tree,
   applying |destroy|, if non-null, to each item in the new tree so far,
   and returns |NULL|.
   If |allocator != NULL|, it is used for allocation in the new tree.
   Otherwise, the same allocator used for |org| is used. */
struct pdavl_table *
pdavl_copy (const struct pdavl_table *org, pdavl_copy_func *copy,
          pdavl_item_func *destroy, struct libpdavl_allocator *allocator)
{
  struct pdavl_node *stack[2 * (PDAVL_MAX_HEIGHT + 1)];
  int height = 0;

  struct pdavl_table *new;
  const struct pdavl_node *x;
  struct pdavl_node *y;

  assert (org != NULL);
  new = pdavl_create (org->pdavl_compare, org->pdavl_param,
                    allocator != NULL ? allocator : org->pdavl_alloc);
  if (new == NULL)
    return NULL;
  new->pdavl_count = org->pdavl_count;
  if (new->pdavl_count == 0)
    return new;

  x = (const struct pdavl_node *) &org->pdavl_root;
  y = (struct pdavl_node *) &new->pdavl_root;
  for (;;)
    {
      while (x->pdavl_link[0] != NULL)
        {
          assert (height < 2 * (PDAVL_MAX_HEIGHT + 1));

          y->pdavl_link[0] =
            new->pdavl_alloc->libpdavl_malloc (new->pdavl_alloc,
                                           sizeof *y->pdavl_link[0]);
          if (y->pdavl_link[0] == NULL)
            {
              if (y != (struct pdavl_node *) &new->pdavl_root)
                {
                  y->pdavl_data = NULL;
                  y->pdavl_link[1] = NULL;
                }

              copy_error_recovery (stack, height, new, destroy);
              return NULL;
            }

          stack[height++] = (struct pdavl_node *) x;
          stack[height++] = y;
          x = x->pdavl_link[0];
          y = y->pdavl_link[0];
        }
      y->pdavl_link[0] = NULL;

      for (;;)
        {
          y->pdavl_balance = x->pdavl_balance;
          if (copy == NULL)
            y->pdavl_data = x->pdavl_data;
          else
            {
              y->pdavl_data = copy (x->pdavl_data, org->pdavl_param);
              if (y->pdavl_data == NULL)
                {
                  y->pdavl_link[1] = NULL;
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }
            }

          if (x->pdavl_link[1] != NULL)
            {
              y->pdavl_link[1] =
                new->pdavl_alloc->libpdavl_malloc (new->pdavl_alloc,
                                               sizeof *y->pdavl_link[1]);
              if (y->pdavl_link[1] == NULL)
                {
                  copy_error_recovery (stack, height, new, destroy);
                  return NULL;
                }

              x = x->pdavl_link[1];
              y = y->pdavl_link[1];
              break;
            }
          else
            y->pdavl_link[1] = NULL;

          if (height <= 2)
            return new;

          y = stack[--height];
          x = stack[--height];
        }
    }
}

/* Frees storage allocated for |tree|.
   If |destroy != NULL|, applies it to each data item in inorder. */
void
pdavl_destroy (struct pdavl_table *tree, pdavl_item_func *destroy)
{
  struct pdavl_node *p, *q;

  assert (tree != NULL);

  for (p = tree->pdavl_root; p != NULL; p = q)
    if (p->pdavl_link[0] == NULL)
      {
        q = p->pdavl_link[1];
        if (destroy != NULL && p->pdavl_data != NULL)
          destroy (p->pdavl_data, tree->pdavl_param);
        tree->pdavl_alloc->libpdavl_free (tree->pdavl_alloc, p);
      }
    else
      {
        q = p->pdavl_link[0];
        p->pdavl_link[0] = q->pdavl_link[1];
        q->pdavl_link[1] = p;
      }

  tree->pdavl_alloc->libpdavl_free (tree->pdavl_alloc, tree);
}

/* Allocates |size| bytes of space using |malloc()|.
   Returns a null pointer if allocation fails. */
void *
pdavl_malloc (struct libpdavl_allocator *allocator, size_t size)
{
  assert (allocator != NULL && size > 0);
  return malloc (size);
}

/* Frees |block|. */
void
pdavl_free (struct libpdavl_allocator *allocator, void *block)
{
  assert (allocator != NULL && block != NULL);
  free (block);
}

/* Default memory allocator that uses |malloc()| and |free()|. */
struct libpdavl_allocator pdavl_allocator_default =
  {
    pdavl_malloc,
    pdavl_free
  };

#undef NDEBUG
#include <assert.h>

/* Asserts that |pdavl_insert()| succeeds at inserting |item| into |table|. */
void
(pdavl_assert_insert) (struct pdavl_table *table, void *item)
{
  void **p = pdavl_probe (table, item);
  assert (p != NULL && *p == item);
}

/* Asserts that |pdavl_delete()| really removes |item| from |table|,
   and returns the removed item. */
void *
(pdavl_assert_delete) (struct pdavl_table *table, void *item)
{
  void *p = pdavl_delete (table, item);
  assert (p != NULL);
  return p;
}


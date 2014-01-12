/*
 *
 * Header information for the calls in mem.c
 *
 *  $Id: octmem.h,v 1.2 2004-10-05 19:10:47 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifndef __MEM_H
#define __MEM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

void reserve_octmem(int bytes_desired,void *mem,int memsize);
void reserve_aux_mem(int bytes_desired);
void free_octmem();
void *oct_alloc(int nelem,int size);
void *aux_alloc(int nelem,int size);

#endif

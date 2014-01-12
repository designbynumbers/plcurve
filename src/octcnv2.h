/*                                                                           *
 * octcnv2.h                                                                 *
 *                                                                           *
 * This file was automatically created by make_octcnv_h (a perl script) and  *
 * can either be recreated with that or edited by hand if it needs updating. *
 *                                                                           *
 * The number of bytes in an int when this file is used is:2                 */ 
/* And some other useful constants:                                          */

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

static int octree_max_bin = 32;
static int octree_max_bits = 5;

/* Here's a table to convert from binary to octal (by spreading the bits) */
static int bin_to_oct[32] = {
         0,          1,          8,          9,         64,         65,
        72,         73,        512,        513,        520,        521,
       576,        577,        584,        585,       4096,       4097,
      4104,       4105,       4160,       4161,       4168,       4169,
      4608,       4609,       4616,       4617,       4672,       4673,
 4680,  4681};

/*
 * @COPYRIGHT@
 *
 * Define TRUE and FALSE
 *
 *  $Id: truefalse.h,v 1.1.1.1 2005-07-01 00:44:41 cantarel Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

/* This file is part of vecttools.
   
vecttools is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

vecttools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vecttools; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifndef __LINKLIB_TRUEFALSE_H
#define __LINKLIB_TRUEFALSE_H

#ifndef FALSE
#define FALSE (1 == 0)
#endif /* FALSE */
#ifndef TRUE
#define TRUE (1 == 1)
#endif /* TRUE */

#endif

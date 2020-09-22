/* glplib1b.c (GNU C version) */

/*----------------------------------------------------------------------
-- Copyright (C) 2000, 2001, 2002, 2003, 2004 Andrew Makhorin,
-- Department for Applied Informatics, Moscow Aviation Institute,
-- Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcnet.ru>.
--
-- This file is part of GLPK (GNU Linear Programming Kit).
--
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
-- 02111-1307, USA.
----------------------------------------------------------------------*/

#include <sys/time.h>
#include "glplib.h"

double lib_get_time(void)
{     struct timeval tv;
      long int tv_sec0;
      double secs;
      /* tv_sec0 = (12:00:00 GMT 1/1/2000) - (00:00:00 GMT 1/1/1970) */
      tv_sec0 = 946728000;
      gettimeofday(&tv, NULL);
      secs = (double)(tv.tv_sec - tv_sec0) + 1e-6 * (double)tv.tv_usec;
      return secs;
}

/* eof */

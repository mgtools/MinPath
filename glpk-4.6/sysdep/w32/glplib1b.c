/* glplib1b.c (W32 version) */

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

#include <windows.h>
#include "glplib.h"

double lib_get_time(void)
{     SYSTEMTIME st;
      FILETIME ft0, ft;
      double secs;
      /* ft0 = 12:00:00 GMT January 1, 2000 */
      ft0.dwLowDateTime  = 0xBAA22000;
      ft0.dwHighDateTime = 0x01BF544F;
      GetSystemTime(&st);
      SystemTimeToFileTime(&st, &ft);
      /* ft = ft - ft0 */
      if (ft.dwLowDateTime >= ft0.dwLowDateTime)
      {  ft.dwLowDateTime  -= ft0.dwLowDateTime;
         ft.dwHighDateTime -= ft0.dwHighDateTime;
      }
      else
      {  ft.dwLowDateTime  += (0xFFFFFFFF - ft0.dwLowDateTime) + 1;
         ft.dwHighDateTime -= ft0.dwHighDateTime + 1;
      }
      secs = (4294967296.0 * (double)(LONG)ft.dwHighDateTime +
         (double)ft.dwLowDateTime) / 10000000.0;
      return secs;
}

/* eof */

/* glplib1a.c (W32 multi-thread DLL version) */

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

static DWORD dwTlsIndex;

#if defined(_MSC_VER)         /* Microsoft Visual C++ */
BOOL APIENTRY DllMain
#elif defined(__BORLANDC__)   /* Borland C++ */
BOOL WINAPI DllEntryPoint
#else
#error Unknown configuration
#endif
(     HINSTANCE hinstDLL,     /* DLL module handle */
      DWORD fdwReason,        /* reason called */
      LPVOID lpvReserved      /* reserved */
)
{     switch (fdwReason)
      {  /* The DLL is loading due to process initialization or a call
            to LoadLibrary. */
         case DLL_PROCESS_ATTACH:
            /* Allocate a TLS index. */
            dwTlsIndex = TlsAlloc();
            if (dwTlsIndex == 0xFFFFFFFF) return FALSE;
            /* Initialize the index for first thread. */
            TlsSetValue(dwTlsIndex, NULL);
#if 0
            /* Initialize GLPK library environment. */
            lib_init_env();
#endif
            break;
         /* The attached process creates a new thread. */
         case DLL_THREAD_ATTACH:
            /* Initialize the TLS index for this thread. */
            TlsSetValue(dwTlsIndex, NULL);
#if 0
            /* Initialize GLPK library environment. */
            lib_init_env();
#endif
            break;
         /* The thread of the attached process terminates. */
         case DLL_THREAD_DETACH:
            /* Free GLPK library environment. */
            lib_free_env();
            break;
         /* The DLL is unloading due to process termination or call to
            FreeLibrary. */
         case DLL_PROCESS_DETACH:
            /* Free GLPK library environment. */
            lib_free_env();
            /* Release the TLS index. */
            TlsFree(dwTlsIndex);
            break;
         default:
            break;
      }
      return TRUE;
}

void lib_set_ptr(void *ptr)
{     TlsSetValue(dwTlsIndex, ptr);
      return;
}

void *lib_get_ptr(void)
{     void *ptr;
      ptr = TlsGetValue(dwTlsIndex);
      return ptr;
}

/* eof */

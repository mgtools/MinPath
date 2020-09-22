/* mtsamp.c */

/*----------------------------------------------------------------------
-- This is an example program that illustrates how to use GLPK API in
-- multi-thread application on W32 platform (VC++ 6.0).
--
-- To run this program use the following command:
--
--    mtsamp.exe mps-file-1 mps-file-2 ... mps-file-n
--
-- Each mps file specified in the command line is processed by separate
-- thread.
--
-- To build multi-thread GLPK DLL the makefile 'w32vc6d.mak' should be
-- used.
----------------------------------------------------------------------*/

#include <stdlib.h>
#include <windows.h>
#include "glpk.h"

DWORD APIENTRY solve_mps(VOID *arg)
{     char *fname = arg;
      LPX *lp;
      lp = lpx_read_mps(fname);
      if (lp == NULL)
      {  print("Cannot read mps file `%s'", fname);
         return 1;
      }
      lpx_simplex(lp);
      lpx_delete_prob(lp);
      return 0;
}

#define MAX_THREADS 20

int main(int argc, char *argv[])
{     HANDLE h[1+MAX_THREADS];
      DWORD foo;
      int t;
      if (argc < 2)
         fault("At least one mps file must be specified");
      if (argc-1 > MAX_THREADS)
         fault("Too many mps files");
      for (t = 1; t < argc; t++)
      {  h[t] = CreateThread(NULL, 0, solve_mps, argv[t], 0, &foo);
         if (h[t] == NULL)
            fault("Unable to create thread for `%s'", argv[t]);
      }
      WaitForMultipleObjects(argc-1, &h[1], TRUE, INFINITE);
      print("GLPK multi-thread DLL seems to work");
      return 0;
}

/* eof */

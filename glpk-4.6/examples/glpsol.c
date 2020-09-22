/* glpsol.c */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "glplib.h"
#include "glplpx.h"
#include "glpmpl.h"

/*--------------------------------------------------------------------*/
/* This program is a stand-alone LP/MIP solver. For pure LP problems  */
/* either the simplex method or the primal-dual interior point method */
/* can be used. For MIP problems the branch-and-bound procedure based */
/* on the simplex method is used.                                     */
/*--------------------------------------------------------------------*/

static char *version = "GLPSOL -- GLPK LP/MIP Solver, Version 4.6";
/* version string */

static int format = 0;
/* type of input text file:
   0 - MPS
   1 - CPLEX LP
   2 - GNU MathProg
   3 - GNU LP */

static char *in_file = NULL;
/* name of input text file */

static char *in_data = NULL;
/* name of optional input text file, which contains data section; NULL
   means no separate data section is provided */

static char *display = NULL;
/* name of optional output text file, to which display output is sent;
   NULL means the output is sent to stdout */

static int dir = 0;
/* optimization direction flag:
   0       - not specified
   LPX_MIN - minimization
   LPX_MAX - maximization */

static int scale = 1;
/* if this flag is set, automatic scaling before solving the problem is
   performed; otherwise scaling is not used */

static int method = 0;
/* which method should be used for solving the problem:
   0 - simplex
   1 - interior point
   2 - branch-and-bound */

static int intopt = 0;
/* if this flag is set, use lpx_intopt rather than lpx_integer */

static char *out_sol = NULL;
/* name of output text file, to which the final solution should be sent
   in plain text format; NULL means no solution output */

static char *out_bnds = NULL;
/* name of output text file, to which sensitivity bounds should be sent
   in plain text format; NULL means no sensitivity output */

static int tmlim = -1;
/* solution time limit, in seconds */

static int check = 0;
/* if this flag is set, only input data checking is required */

static int orig = 0;
/* if this flag is set, try to use original names of rows and columns;
   otherwise use plain names */

static char *out_mps = NULL;
/* name of output text file, to which the problem should be written in
   MPS format; NULL means no MPS output */

static char *out_lpt = NULL;
/* name of output text file, to which the problem should be written in
   CPLEX LP format; NULL means no CPLEX LP output */

static char *out_txt = NULL;
/* name of output text file, to which the problem should be written in
   plain text format; NULL means no plain text output */

static char *out_glp = NULL;
/* name of output text file, to which the problem should be written in
   GNU LP format; NULL means no GNU LP output */

static char *newname = NULL;
/* new name which has to be assigned to the problem */

static int basis = 1;
/* which initial basis should be used:
   0 - standard initial basis
   1 - advanced initial basis */

static int price = 1;
/* which pricing technique should be used:
   0 - textbook pricing
   1 - steepest edge pricing */

static int relax = 1;
/* if this flag is set, the solver uses two-pass ratio test (for both
   primal and dual simplex) proposed by P.Harris; otherwise the solver
   uses the standard "textbook" ratio test */

static int presol = 1;
/* if this flag is set, the solver uses the LP presolver; otherwise the
   LP presolver is not used */

static int nomip = 0;
/* if this flag is set, the solver considers all integer variables as
   continuous (this allows solving MIP problem as pure LP) */

static int branch = 2;
/* which branching technique should be used:
   0 - on the first integer variable
   1 - on the last integer variable
   2 - using heuristic by Driebeek and Tomlin */

static int btrack = 2;
/* which backtracking technique should be used:
   0 - depth first search
   1 - breadth first search
   2 - the best projection heuristic */

/*----------------------------------------------------------------------
-- display_help - display help.
--
-- This routine displays help information about the program as required
-- by the GNU Coding Standards. */

static void display_help(char *my_name)
{     print("Usage: %s [options...] filename", my_name);
      print("");
      print("General options:");
      print("   --glp             read LP/MIP model in GNU LP format");
      print("   --mps             read LP/MIP problem in MPS format (de"
         "fault)");
      print("   --lpt             read LP/MIP problem in CPLEX LP forma"
         "t");
      print("   --math            read LP/MIP model written in GNU Math"
         "Prog modeling");
      print("                     language");
      print("   -m filename, --model filename");
      print("                     read model section and optional data "
         "section from");
      print("                     filename (the same as --math)");
      print("   -d filename, --data filename");
      print("                     read data section from filename (for "
         "--math only);");
      print("                     if model file also has data section, "
         "that section");
      print("                     is ignored");
      print("   -y filename, --display filename");
      print("                     send display output to filename (for "
         "--math only);");
      print("                     by default the output is sent to stdo"
         "ut");
      print("   --min             minimization");
      print("   --max             maximization");
      print("   --scale           scale problem (default)");
      print("   --noscale         do not scale problem");
      print("   --simplex         use simplex method (default)");
      print("   --interior        use interior point method (for pure L"
         "P only)");
      print("   -o filename, --output filename");
      print("                     write solution to filename in plain t"
         "ext format");
      print("   --bounds filename");
      print("                     write sensitivity bounds to filename "
         "in plain ");
      print("                     text format (LP only)");
      print("   --tmlim nnn       limit solution time to nnn seconds");
      print("   --check           do not solve problem, check input dat"
         "a only");
      print("   --name probname   change problem name to probname");
      print("   --plain           use plain names of rows and columns ("
         "default)");
      print("   --orig            try using original names of rows and "
         "columns");
      print("   --wglp filename   write problem to filename in GNU LP f"
         "ormat");
      print("   --wmps filename   write problem to filename in MPS form"
         "at");
      print("   --wlpt filename   write problem to filename in CPLEX LP"
         " format");
      print("   --wtxt filename   write problem to filename in plain te"
         "xt format");
      print("   -h, --help        display this help information and exi"
         "t");
      print("   -v, --version     display program version and exit");
      print("");
      print("Options specific to simplex method:");
      print("   --std             use standard initial basis of all sla"
         "cks");
      print("   --adv             use advanced initial basis (default)")
         ;
      print("   --steep           use steepest edge technique (default)"
         );
      print("   --nosteep         use standard \"textbook\" pricing");
      print("   --relax           use Harris' two-pass ratio test (defa"
         "ult)");
      print("   --norelax         use standard \"textbook\" ratio test")
         ;
      print("   --presol          use presolver (default; assumes --sca"
         "le and --adv)");
      print("   --nopresol        do not use presolver");
      print("");
      print("Options specific to MIP:");
      print("   --nomip           consider all integer variables as con"
         "tinuous");
      print("                     (allows solving MIP as pure LP)");
      print("   --first           branch on first integer variable");
      print("   --last            branch on last integer variable");
      print("   --drtom           branch using heuristic by Driebeck an"
         "d Tomlin");
      print("                     (default)");
      print("   --dfs             backtrack using depth first search");
      print("   --bfs             backtrack using breadth first search")
         ;
      print("   --bestp           backtrack using the best projection h"
         "euristic");
      print("                     (default)");
      print("");
      print("For description of the MPS and CPLEX LP formats see Refere"
         "nce Manual.");
      print("For description of the modeling language see \"GLPK: Model"
         "ing Language");
      print("GNU MathProg\". Both documents are included in the GLPK di"
         "stribution.");
      print("");
      print("See GLPK web page at <http://www.gnu.org/software/glpk/glp"
         "k.html>.");
      print("");
      print("Please report bugs to <bug-glpk@gnu.org>.");
      exit(EXIT_SUCCESS);
      /* no return */
}

/*----------------------------------------------------------------------
-- display_version - display version.
--
-- This routine displays version information for the program as required
-- by the GNU Coding Standards. */

static void display_version(void)
{     print("%s", version);
      print("Copyright (C) 2000, 01, 02, 03, 04 Andrew Makhorin <mao@ma"
         "i2.rcnet.ru>");
      print("This program is free software; you may redistribute it und"
         "er the terms of");
      print("the GNU General Public License. This program has absolutel"
         "y no warranty.");
      exit(EXIT_SUCCESS);
      /* no return */
}

/*----------------------------------------------------------------------
-- parse_cmdline - parse command-line parameters.
--
-- This routine parses parameters specified in the command line. */

#define p(str) (strcmp(argv[k], str) == 0)

static void parse_cmdline(int argc, char *argv[])
{     int k;
      for (k = 1; k < argc; k++)
      {  if (p("--mps"))
            format = 0;
         else if (p("--lpt"))
            format = 1;
         else if (p("--math") || p("-m") || p("--model"))
            format = 2;
         else if (p("--glp"))
            format = 3;
         else if (p("-d") || p("--data"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No data input file specifed");
               exit(EXIT_FAILURE);
            }
            if(in_data != NULL)
            {  print("Only one data input file allowed");
               exit(EXIT_FAILURE);
            }
            in_data = argv[k];
         }
         else if (p("-y") || p("--display"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No display output file specifed");
               exit(EXIT_FAILURE);
            }
            if (display != NULL)
            {  print("Only one display output file allowed");
               exit(EXIT_FAILURE);
            }
            display = argv[k];
         }
         else if (p("--min"))
            dir = LPX_MIN;
         else if (p("--max"))
            dir = LPX_MAX;
         else if (p("--scale"))
            scale = 1;
         else if (p("--noscale"))
            scale = 0;
         else if (p("--simplex"))
            method = 0;
         else if (p("--interior"))
            method = 1;
         else if (p("--intopt"))
            intopt = 1;
         else if (p("-o") || p("--output"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No solution output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_sol != NULL)
            {  print("Only one solution output file allowed");
               exit(EXIT_FAILURE);
            }
            out_sol = argv[k];
         }
         else if (p("--bounds"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No sensitivity bounds output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_bnds != NULL)
            {  print("Only one sensitivity bounds output file allowed");
               exit(EXIT_FAILURE);
            }
            out_bnds = argv[k];
         }
         else if (p("--tmlim"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No time limit specified");
               exit(EXIT_FAILURE);
            }
            if (str2int(argv[k], &tmlim) || tmlim < 0)
            {  print("Invalid time limit `%s'", argv[k]);
               exit(EXIT_FAILURE);
            }
         }
         else if (p("--check"))
            check = 1;
         else if (p("--name"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem name specified");
               exit(EXIT_FAILURE);
            }
            if (newname != NULL)
            {  print("Only one problem name allowed");
               exit(EXIT_FAILURE);
            }
            newname = argv[k];
         }
         else if (p("--plain"))
            orig = 0;
         else if (p("--orig"))
            orig = 1;
         else if (p("--wmps"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No MPS output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_mps != NULL)
            {  print("Only one MPS output file allowed");
               exit(EXIT_FAILURE);
            }
            out_mps = argv[k];
         }
         else if (p("--wlpt"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No CPLEX LP output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_lpt != NULL)
            {  print("Only one CPLEX LP output file allowed");
               exit(EXIT_FAILURE);
            }
            out_lpt = argv[k];
         }
         else if (p("--wtxt"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_txt != NULL)
            {  print("Only one problem output file allowed");
               exit(EXIT_FAILURE);
            }
            out_txt = argv[k];
         }
         else if (p("--wglp"))
         {  k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {  print("No problem output file specified");
               exit(EXIT_FAILURE);
            }
            if (out_glp != NULL)
            {  print("Only one problem output file allowed");
               exit(EXIT_FAILURE);
            }
            out_glp = argv[k];
         }
         else if (p("-h") || p("--help"))
            display_help(argv[0]);
         else if (p("-v") || p("--version"))
            display_version();
         else if (p("--std"))
            basis = 0;
         else if (p("--adv"))
            basis = 1;
         else if (p("--steep"))
            price = 1;
         else if (p("--nosteep"))
            price = 0;
         else if (p("--relax"))
            relax = 1;
         else if (p("--norelax"))
            relax = 0;
         else if (p("--presol"))
            presol = 1;
         else if (p("--nopresol"))
            presol = 0;
         else if (p("--nomip"))
            nomip = 1;
         else if (p("--first"))
            branch = 0;
         else if (p("--last"))
            branch = 1;
         else if (p("--drtom"))
            branch = 2;
         else if (p("--dfs"))
            btrack = 0;
         else if (p("--bfs"))
            btrack = 1;
         else if (p("--bestp"))
            btrack = 2;
         else if (argv[k][0] == '-' ||
                 (argv[k][0] == '-' && argv[k][1] == '-'))
         {  print("Invalid option `%s'; try %s --help",
               argv[k], argv[0]);
            exit(EXIT_FAILURE);
         }
         else
         {  if (in_file != NULL)
            {  print("Only one input file allowed");
               exit(EXIT_FAILURE);
            }
            in_file = argv[k];
         }
      }
      return;
}

/*----------------------------------------------------------------------
-- main - main program.
--
-- This main program is called by the control program and manages the
-- solving process. */

int main(int argc, char *argv[])
{     LPX *lp;
      MPL *mpl = NULL;
      int ret;
      double start;
      /* parse command line parameters */
      parse_cmdline(argc, argv);
      /* remove all output files specified in the command line */
      if (display != NULL) remove(display);
      if (out_sol != NULL) remove(out_sol);
      if (out_bnds != NULL) remove(out_bnds);
      if (out_mps != NULL) remove(out_mps);
      if (out_lpt != NULL) remove(out_lpt);
      if (out_txt != NULL) remove(out_txt);
      if (out_glp != NULL) remove(out_glp);
      /* read problem from the input file */
      if (in_file == NULL)
      {  print("No input file specified; try %s --help", argv[0]);
         exit(EXIT_FAILURE);
      }
      switch (format)
      {  case 0:
            lp = lpx_read_mps(in_file);
            if (lp == NULL)
            {  print("MPS file processing error");
               exit(EXIT_FAILURE);
            }
            break;
         case 1:
            lp = lpx_read_lpt(in_file);
            if (lp == NULL)
            {  print("CPLEX LP file processing error");
               exit(EXIT_FAILURE);
            }
            break;
         case 2:
#if 0 /* 01/VIII-2004 */
            lp = lpx_read_model(in_file, in_data, display);
            if (lp == NULL)
            {  print("Model processing error");
               exit(EXIT_FAILURE);
            }
#else
            /* initialize the translator database */
            mpl = mpl_initialize();
            /* read model section and optional data section */
            ret = mpl_read_model(mpl, in_file, in_data != NULL);
            if (ret == 4)
err:        {  print("Model processing error");
               exit(EXIT_FAILURE);
            }
            insist(ret == 1 || ret == 2);
            /* read data section, if necessary */
            if (in_data != NULL)
            {  insist(ret == 1);
               ret = mpl_read_data(mpl, in_data);
               if (ret == 4) goto err;
               insist(ret == 2);
            }
            /* generate model */
            ret = mpl_generate(mpl, display);
            if (ret == 4) goto err;
            /* extract problem instance */
            lp = lpx_extract_prob(mpl);
            insist(lp != NULL);
#endif
            if (lpx_get_num_rows(lp) == 0)
            {  print("Problem has no rows");
               exit(EXIT_FAILURE);
            }
            if (lpx_get_num_cols(lp) == 0)
            {  print("Problem has no columns");
               exit(EXIT_FAILURE);
            }
            break;
         case 3:
            lp = lpx_read_prob(in_file);
            if (lp == NULL)
            {  print("GNU LP file processing error");
               exit(EXIT_FAILURE);
            }
            break;
         default:
            insist(format != format);
      }
      /* change problem name (if required) */
      if (newname != NULL) lpx_set_prob_name(lp, newname);
      /* change optimization direction (if required) */
      if (dir != 0) lpx_set_obj_dir(lp, dir);
      /* write problem in MPS format (if required) */
      if (out_mps != NULL)
      {  lpx_set_int_parm(lp, LPX_K_MPSORIG, orig);
         ret = lpx_write_mps(lp, out_mps);
         if (ret != 0)
         {  print("Unable to write problem in MPS format");
            exit(EXIT_FAILURE);
         }
      }
      /* write problem in CPLEX LP format (if required) */
      if (out_lpt != NULL)
      {  lpx_set_int_parm(lp, LPX_K_LPTORIG, orig);
         ret = lpx_write_lpt(lp, out_lpt);
         if (ret != 0)
         {  print("Unable to write problem in CPLEX LP format");
            exit(EXIT_FAILURE);
         }
      }
      /* write problem in plain text format (if required) */
      if (out_txt != NULL)
      {  lpx_set_int_parm(lp, LPX_K_LPTORIG, orig);
         ret = lpx_print_prob(lp, out_txt);
         if (ret != 0)
         {  print("Unable to write problem in plain text format");
            exit(EXIT_FAILURE);
         }
      }
      /* write problem in GNU LP format (if required) */
      if (out_glp != NULL)
      {  ret = lpx_write_prob(lp, out_glp);
         if (ret != 0)
         {  print("Unable to write problem in GNU LP format");
            exit(EXIT_FAILURE);
         }
      }
      /* if only data check is required, skip computations */
      if (check) goto skip;
      /* scale the problem data (if required) */
      if (scale && (!presol || method == 1)) lpx_scale_prob(lp);
      /* build advanced initial basis (if required) */
      if (method == 0 && basis && !presol) lpx_adv_basis(lp);
      /* set some control parameters, which might be changed in the
         command line */
      lpx_set_int_parm(lp, LPX_K_PRICE, price);
      if (!relax) lpx_set_real_parm(lp, LPX_K_RELAX, 0.0);
      lpx_set_int_parm(lp, LPX_K_PRESOL, presol);
      lpx_set_int_parm(lp, LPX_K_BRANCH, branch);
      lpx_set_int_parm(lp, LPX_K_BTRACK, btrack);
      lpx_set_real_parm(lp, LPX_K_TMLIM, (double)tmlim);
      /* solve the problem */
      start = utime();
      switch (method)
      {  case 0:
            if (nomip || lpx_get_class(lp) == LPX_LP)
            {  ret = lpx_simplex(lp);
               if (presol && ret != LPX_E_OK && out_sol != NULL)
                  print("If you need actual output for non-optimal solu"
                     "tion, use --nopresol");
            }
            else
            {  method = 2;
               lpx_simplex(lp);
               if (!intopt)
                  lpx_integer(lp);
               else
                  lpx_intopt(lp);
            }
            break;
         case 1:
            if (nomip || lpx_get_class(lp) == LPX_LP)
               lpx_interior(lp);
            else
            {  print("Interior point method is not able to solve MIP pr"
                  "oblem; use --simplex");
               exit(EXIT_FAILURE);
            }
            break;
         default:
            insist(method != method);
      }
      /* display statistics */
      print("Time used:   %.1f secs", utime() - start);
      print("Memory used: %.1fM (%d bytes)",
         (double)lib_env_ptr()->mem_tpeak / (double)(1024 * 1024),
         lib_env_ptr()->mem_tpeak);
#if 1 /* 01/VIII-2004 */
      if (mpl != NULL && mpl_has_solve_stmt(mpl))
      {  int n, j, round;
         /* store the solution to the translator database */
         n = lpx_get_num_cols(lp);
         round = lpx_get_int_parm(lp, LPX_K_ROUND);
         lpx_set_int_parm(lp, LPX_K_ROUND, 1);
         switch (method)
         {  case 0:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_get_col_prim(lp, j));
               break;
            case 1:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_ipt_col_prim(lp, j));
               break;
            case 2:
               for (j = 1; j <= n; j++)
                  mpl_put_col_value(mpl, j, lpx_mip_col_val(lp, j));
               break;
            default:
               insist(method != method);
         }
         lpx_set_int_parm(lp, LPX_K_ROUND, round);
         /* perform postsolving */
         ret = mpl_postsolve(mpl, display);
         if (ret == 4)
         {  print("Model postsolving error");
            exit(EXIT_FAILURE);
         }
         insist(ret == 3);
      }
#endif
      /* write problem solution found by the solver (if required) */
      if (out_sol != NULL)
      {  switch (method)
         {  case 0:
               ret = lpx_print_sol(lp, out_sol);
               break;
            case 1:
               ret = lpx_print_ips(lp, out_sol);
               break;
            case 2:
               ret = lpx_print_mip(lp, out_sol);
               break;
            default:
               insist(method != method);
         }
         if (ret != 0)
         {  print("Unable to write problem solution");
            exit(EXIT_FAILURE);
         }
      }
      /* write sensitivity bounds information (if required) */
      if (out_bnds != NULL)
      {  if (method != 0)
         {  print("Cannot write sensitivity bounds information for inte"
               "rior-point or MIP solution");
            exit(EXIT_FAILURE);
         }
         ret = lpx_print_sens_bnds(lp, out_bnds);
         if (ret != 0)
         {  print("Unable to write sensitivity bounds information");
            exit(EXIT_FAILURE);
         }
      }
skip: /* delete the problem object */
      lpx_delete_prob(lp);
#if 1 /* 01/VIII-2004 */
      /* if the translator database exists, destroy it */
      if (mpl != NULL) mpl_terminate(mpl);
#endif
      /* check that no memory blocks are still allocated */
      insist(lib_env_ptr()->mem_total == 0);
      insist(lib_env_ptr()->mem_count == 0);
      /* return to the control program */
      return 0;
}

/* eof */

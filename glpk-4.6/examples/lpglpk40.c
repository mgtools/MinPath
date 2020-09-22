/* lpglpk40.c */

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

#include <float.h>
#include <string.h>
#include "glpk.h"
#include "lp.h"

/*----------------------------------------------------------------------
-- This is a non-trivial example of using GLPK as a base LP solver for
-- Concorde, a program for solving Traveling Salesman Problem (TSP).
--
-- If you wish to try GLPK and Concorde together, do the following:
--
-- 1. Install GLPK (see the file 'INSTALL' in GLPK distribution).
--
-- 2. Download Concorde from
--    http://www.math.princeton.edu/tsp/concorde.html
--
--    Note that the GLPK interface routines were tested with Concorde
--    version 991215.
--
-- 3. Unpack and unarchive Concorde tgz-file.
--
-- 4. Copy the file you are reading into the subdirectory 'concorde/LP/'
--    keeping its original name 'lpglpk40.c'.
--
-- 5. Configure Concorde in the usual way and then build Concorde using
--    the following command:
--
--       make LPSOLVER_INTERFACE=lpglpk40.c LPSOLVER_LIB=-lglpk
--
--    Note that all the command-line parameters are case sensitive.
--
-- 6. Run Concorde to solve your TSP instances. See also for TSPLIB 95,
--    a library of standard TSP instances; this library is available at
--    http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/
--    (some results given in the report tsplib95.ps are wrong; for more
--    correct results see the webpage .../TSPLIB95/stsp.html on the site
--    above).
--
-- Note that these interface routines do not include many important
-- features (for example, efficient strong branching), therefore hard
-- TSP instances cannot be solved in a reasonable time.
----------------------------------------------------------------------*/

static int MSGLEV = 0;
/* 0 - no output
   1 - output one line per one call to GLPK simplex solver
   2 - detailed output (for debugging) */

static char *SOLVER = "GLPK";

#define AT_LOWER 0
#define IS_BASIC 1
#define AT_UPPER 2

struct CClp { LPX *lp; };

struct CClp_warmstart { int nrows, ncols, *rstat, *cstat; };

struct CClp_info { int nrows, ncols, *rstat, *cstat; };

/*--------------------------------------------------------------------*/

int CClp_init(CClp **lp)
{     /* INITIALIZES the LP. */
      CClp_free(lp);
      (*lp) = umalloc(sizeof(CClp));
      (*lp)->lp = NULL;
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_force_perturb(CClp *lp)
{     /* FORCES a perturbation in the LP simplex solves. */
      insist(lp == lp);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_tune_small(CClp *lp)
{     /* SETS solver options for tiny problems. */
      insist(lp == lp);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_disable_presolve(CClp *lp)
{     /* DISABLES the solvers presolve. */
      insist(lp == lp);
      return 0;
}

/*--------------------------------------------------------------------*/

void CClp_free(CClp **lp)
{     /* FREES the LP, both the allocation in CClp_init () and the
         allocation in CClp_loadlp. */
      if ((*lp) != NULL)
      {  CClp_freelp(lp);
         ufree((*lp)), (*lp) = NULL;
      }
      return;
}

/*--------------------------------------------------------------------*/

void CClp_freelp(CClp **lp)
{     /* FREES the LP loaded in CClp_loadlp (or CClp_create), but does
         not free the data allocated by CClp_init. */
      if ((*lp) != NULL && (*lp)->lp != NULL)
         lpx_delete_prob((*lp)->lp), (*lp)->lp = NULL;
      return;
}

/*--------------------------------------------------------------------*/

int CClp_loadlp(CClp *lp, const char *name, int ncols, int nrows,
      int objsense, double *obj, double *rhs, char *sense, int *matbeg,
      int *matcnt, int *matind, double *matval, double *lb, double *ub)
{     /* LOADS the data into the LP.
         - name attaches a name to the LP (it can be used by the LP
           solver in io routines)
         - ncols and nrows give the number of columns and rows in the LP
         - objsense should be 1 for minimize and -1 for maximize
         - obj and rhs are arrays giving the objective function and rhs
         - sense is an array specifying 'L', 'E', or 'G' for each of the
           rows
         - matbeg, matcnt, matind, and matval give the coefficients of
           the constraint matrix in column by column order.
           matbeg gives the index of the start of each column;
           matcnt gives the number of coefficients in each column;
           matind gives the indices of the rows where the coefficients
           are located in the constraint matrix (so for column j, the
           indices are given in matcnt[j] locations starting at
           matind[matbeg[j]]; and matval gives the actual coefficients
           (organized like matind).
         - lb and ub are arrays giving the upper and lower bounds of the
           variables. */
      int i, j;
      /* create empty problem object */
      insist(lp->lp == NULL);
      lp->lp = lpx_create_prob();
      lpx_set_prob_name(lp->lp, (char *)name);
      /* set objective sense */
      switch (objsense)
      {  case +1:
            /* minimization */
            lpx_set_obj_dir(lp->lp, LPX_MIN); break;
         case -1:
            /* maximization */
            lpx_set_obj_dir(lp->lp, LPX_MAX); break;
         default:
            insist(objsense != objsense);
      }
      /* add rows */
      lpx_add_rows(lp->lp, nrows);
      for (i = 0; i < nrows; i++)
      {  int seqn, type;
         double lo, up;
         seqn = i+1;
         switch (sense[i])
         {  case 'L':
               type = LPX_UP, lo = 0.0, up = rhs[i];
               break;
            case 'E':
               type = LPX_FX, lo = up = rhs[i];
               break;
            case 'G':
               type = LPX_LO, lo = rhs[i], up = 0.0;
               break;
            default:
               insist(sense[i] != sense[i]);
         }
         lpx_set_row_bnds(lp->lp, seqn, type, lo, up);
      }
      /* add columns and constraint coefficients */
      lpx_add_cols(lp->lp, ncols);
      for (j = 0; j < ncols; j++)
      {  int seqn, type, k;
         double lo, up;
         seqn = j+1;
         lpx_set_col_coef(lp->lp, seqn, obj == NULL ? 0.0 : obj[j]);
         lo = lb[j], up = ub[j];
         /* check for finite bounds */
         insist(-1e12 <= lo && lo <= up && up <= +1e12);
         type = (lo == up ? LPX_FX : LPX_DB);
         lpx_set_col_bnds(lp->lp, seqn, type, lo, up);
         for (k = matbeg[j]; k < matbeg[j] + matcnt[j]; k++)
            matind[k]++;
         lpx_set_mat_col(lp->lp, seqn, matcnt[j],
            &matind[matbeg[j]] - 1, &matval[matbeg[j]] - 1);
         for (k = matbeg[j]; k < matbeg[j] + matcnt[j]; k++)
            matind[k]--;
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_create(CClp *lp, const char *name)
{     /* CREATES an empty lp. This supports an alternative to
         CClp_loadlp for loading a problem.
         - name attaches a name to the LP (it can be used by the LP
           solver in io routines). */
      insist(lp->lp == NULL);
      lp->lp = lpx_create_prob();
      lpx_set_prob_name(lp->lp, (char *)name);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_new_row(CClp *lp, char sense, double rhs)
{     /* ADDS a new empty row to the lp.
         - sense is 'L', 'E', or 'G' for a <=, =, or >= constraint;
         - rhs is the right-hand side of the row. */
      int seqn;
      lpx_add_rows(lp->lp, 1);
      seqn = lpx_get_num_rows(lp->lp);
      switch (sense)
      {  case 'L':
            lpx_set_row_bnds(lp->lp, seqn, LPX_UP, 0.0, rhs);
            break;
         case 'E':
            lpx_set_row_bnds(lp->lp, seqn, LPX_FX, rhs, rhs);
            break;
         case 'G':
            lpx_set_row_bnds(lp->lp, seqn, LPX_LO, rhs, 0.0);
            break;
         default:
            insist(sense != sense);
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_change_sense(CClp *lp, int row, char sense)
{     /* CHANGES the sense of a row.
         - row is the row number to change;
         - sense is 'L', 'E', or 'G' to change to <=, =, or >=. */
      int type;
      double lo, up, rhs;
      lpx_get_row_bnds(lp->lp, row+1, &type, &lo, &up);
      switch (type)
      {  case LPX_LO:
            rhs = lo; break;
         case LPX_UP:
            rhs = up; break;
         case LPX_FX:
            rhs = lo; break;
         default:
            insist(type != type);
      }
      switch (sense)
      {  case 'L':
            type = LPX_UP, lo = 0.0, up = rhs; break;
         case 'E':
            type = LPX_FX, lo = up = rhs; break;
         case 'G':
            type = LPX_LO, lo = rhs, up = 0.0; break;
         default:
            insist(sense != sense);
      }
      lpx_set_row_bnds(lp->lp, row+1, type, lo, up);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_opt(CClp *lp, int method)
{     /* CALLS designated LP solution method. */
      int stat, ret;
      if (MSGLEV >= 1)
      {  int m = lpx_get_num_rows(lp->lp);
         int n = lpx_get_num_cols(lp->lp);
         int nz = lpx_get_num_nz(lp->lp);
         print("CClp_opt: %-11s m = %d; n = %d; nz = %d",
            method == CClp_METHOD_DUAL ? "(dual)" : "(primal)",
            m, n, nz);
      }
      lpx_set_int_parm(lp->lp, LPX_K_DUAL, method == CClp_METHOD_DUAL);
      switch (MSGLEV)
      {  case 0:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 0);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 1000000);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 1e6);
            break;
         case 1:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 2);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 200);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 5.0);
            break;
         case 2:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 3);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 200);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 0.0);
            break;
         default:
            insist(MSGLEV != MSGLEV);
      }
      ret = lpx_simplex(lp->lp);
      if (ret == LPX_E_FAULT)
      {  if (MSGLEV >= 1) print("CClp_opt: restarting from advanced bas"
            "is...");
         lpx_adv_basis(lp->lp);
         ret = lpx_simplex(lp->lp);
      }
      if (ret != LPX_E_OK)
      {  print("CClp_opt: lpx_simplex failed; return code = %d", ret);
         ret = 1;
         goto done;
      }
      stat = lpx_get_status(lp->lp);
      if (stat == LPX_OPT)
         ret = 0;
      else if (stat == LPX_NOFEAS)
         ret = 2;
      else
      {  print("CClp_opt: optimization status = %d", stat);
         ret = 1;
      }
done: return ret;
}

/*--------------------------------------------------------------------*/

int CClp_limited_dualopt(CClp *lp, int iterationlim, int *status,
      double *objupperlim)
{     /* CALLS the dual simplex method with a limit on the number of
         pivots.
         - upperbound it is used to cutoff the dual simplex method (when
           the objective value reaches upperbound); it can be NULL;
         - status returns the status of the optimization (it can be
           NULL). */
      int stat, ret;
      insist(iterationlim == iterationlim);
      insist(objupperlim == objupperlim);
      if (MSGLEV >= 1)
      {  int m = lpx_get_num_rows(lp->lp);
         int n = lpx_get_num_cols(lp->lp);
         int nz = lpx_get_num_nz(lp->lp);
         print("CClp_limited_dualopt: m = %d; n = %d; nz = %d", m, n,
            nz);
      }
      lpx_set_int_parm(lp->lp, LPX_K_DUAL, 1);
      switch (MSGLEV)
      {  case 0:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 0);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 1000000);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 1e6);
            break;
         case 1:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 2);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 200);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 5.0);
            break;
         case 2:
            lpx_set_int_parm(lp->lp, LPX_K_MSGLEV, 3);
            lpx_set_int_parm(lp->lp, LPX_K_OUTFRQ, 200);
            lpx_set_real_parm(lp->lp, LPX_K_OUTDLY, 0.0);
            break;
         default:
            insist(MSGLEV != MSGLEV);
      }
      ret = lpx_simplex(lp->lp);
      if (ret == LPX_E_FAULT)
      {  if (MSGLEV >= 1) print("CClp_limited_dualopt: restarting from "
            "advanced basis...");
         lpx_adv_basis(lp->lp);
         ret = lpx_simplex(lp->lp);
      }
      if (ret != LPX_E_OK)
      {  print("CClp_limited_dualopt: lpx_simplex failed; return code ="
            " %d", ret);
         if (status) *status = CClp_FAILURE;
         ret = 1;
         goto done;
      }
      stat = lpx_get_status(lp->lp);
      if (stat == LPX_OPT)
      {  if (status) *status = CClp_SUCCESS;
         ret = 0;
      }
      else if (stat == LPX_NOFEAS)
      {  if (status) *status = CClp_INFEASIBLE;
         ret = 0;
      }
      else
      {  print("CClp_limited_dualopt: optimization status = %d", stat);
         if (status) *status = CClp_FAILURE;
         ret = 1;
      }
done: return ret;
}

/*--------------------------------------------------------------------*/

int CClp_addrows(CClp *lp, int newrows, int newnz, double *rhs,
      char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{     /* ADDS the rows to the LP.
         - newrows is the number of rows to be added;
         - newnz is the number of nonzero coefficients in the new rows;
         - rhs is an array of the rhs values for the new rows;
         - sense is 'L', 'E', or 'G' for each of the new rows;
         - rmatbeg, rmatind, and rmatval give the coefficients of the
           new rows in sparse format. The arrays can be freed after the
           call. */
      int i;
      lpx_add_rows(lp->lp, newrows);
      for (i = 0; i < newrows; i++)
      {  int seqn, type, k, t;
         double lo, up;
         seqn = lpx_get_num_rows(lp->lp) - newrows + i + 1;
         switch (sense[i])
         {  case 'L':
               type = LPX_UP, lo = 0.0, up = rhs[i];
               break;
            case 'E':
               type = LPX_FX, lo = up = rhs[i];
               break;
            case 'G':
               type = LPX_LO, lo = rhs[i], up = 0.0;
               break;
            default:
               insist(sense[i] != sense[i]);
         }
         lpx_set_row_bnds(lp->lp, seqn, type, lo, up);
         insist(rmatbeg != NULL);
         insist(rmatind != NULL);
         insist(rmatval != NULL);
         insist(rmatbeg[0] == 0);
         t = (i < newrows-1 ? rmatbeg[i+1] : newnz);
         for (k = rmatbeg[i]; k < t; k++) rmatind[k]++;
         lpx_set_mat_row(lp->lp, seqn, t - rmatbeg[i],
            &rmatind[rmatbeg[i]] - 1, &rmatval[rmatbeg[i]] - 1);
         for (k = rmatbeg[i]; k < t; k++) rmatind[k]--;
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_addcols(CClp *lp, int newcols, int newnz, double *obj,
      int *cmatbeg, int *cmatind, double *cmatval, double *lb,
      double *ub)
{     /* ADDS the columns to the LP. */
      int j;
      lpx_add_cols(lp->lp, newcols);
      for (j = 0; j < newcols; j++)
      {  int seqn, type, k, t;
         double lo, up;
         seqn = lpx_get_num_cols(lp->lp) - newcols + j + 1;
         lpx_set_col_coef(lp->lp, seqn, obj == NULL ? 0.0 : obj[j]);
         lo = lb[j], up = ub[j];
         /* check for finite bounds */
         insist(-1e10 <= lo && lo <= up && up <= +1e10);
         type = (lo == up ? LPX_FX : LPX_DB);
         lpx_set_col_bnds(lp->lp, seqn, type, lo, up);
         insist(cmatbeg != NULL);
         insist(cmatind != NULL);
         insist(cmatval != NULL);
         insist(cmatbeg[0] == 0);
         t = (j < newcols-1 ? cmatbeg[j+1] : newnz);
         for (k = cmatbeg[j]; k < t; k++) cmatind[k]++;
         lpx_set_mat_col(lp->lp, seqn, t - cmatbeg[j],
            &cmatind[cmatbeg[j]] - 1, &cmatval[cmatbeg[j]] - 1);
         for (k = cmatbeg[j]; k < t; k++) cmatind[k]--;
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_delete_row(CClp *lp, int i)
{     /* DELETES row i of the LP. */
      int num[1+1];
      num[1] = i+1;
      lpx_del_rows(lp->lp, 1, num);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_delete_set_of_rows(CClp *lp, int *delstat)
{     /* DELETES the rows corresponding to 1 entries in delstat.
         - delstat is a 0/1 array having an entry for each row. */
      int m = lpx_get_num_rows(lp->lp);
      int nrs = 0, i;
      int *num = ucalloc(1+m, sizeof(int));
      for (i = 0; i < m; i++)
         if (delstat[i]) num[++nrs] = i+1;
      if (nrs > 0) lpx_del_rows(lp->lp, nrs, num);
      ufree(num);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_delete_column(CClp *lp, int j)
{     /* DELETES column j from the LP. */
      int num[1+1];
      num[1] = j+1;
      lpx_del_cols(lp->lp, 1, num);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_delete_set_of_columns(CClp *lp, int *delstat)
{     /* DELETES the columns corresponding to the 1 entries in delstat.
         - delstat is a 0/1 array having an entry for each column. */
      insist(lp == lp);
      insist(delstat == delstat);
      fault("CClp_delete_set_of_columns: not implemented");
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_setbnd(CClp *lp, int j, char lower_or_upper, double bnd)
{     /* SETS the bound on the variable index by j.
         - lower_or_upper should be either 'L' or 'U'. */
      int type;
      double lo, up;
      lpx_get_col_bnds(lp->lp, j+1, &type, &lo, &up);
      if (type == LPX_FR || type == LPX_UP) lo = -DBL_MAX;
      if (type == LPX_FR || type == LPX_LO) up = +DBL_MAX;
      switch (lower_or_upper)
      {  case 'L':
            lo = bnd; break;
         case 'U':
            up = bnd; break;
         default:
            insist(lower_or_upper != lower_or_upper);
      }
      if (lo == -DBL_MAX && up == +DBL_MAX)
         insist(2 + 2 == 5);
      else if (up == +DBL_MAX)
         type = LPX_LO;
      else if (lo == -DBL_MAX)
         type = LPX_UP;
      else if (lo != up)
         type = LPX_DB;
      else
         type = LPX_FX;
      lpx_set_col_bnds(lp->lp, j+1, type, lo, up);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_get_warmstart(CClp *lp, CClp_warmstart **warm)
{     /* SAVES information for efficiently resolving the current lp in
         warm, for example, basis or norm information. */
      int i, j, tagx;
      CClp_free_warmstart(warm);
      (*warm) = umalloc(sizeof(CClp_warmstart));
      (*warm)->ncols = 0;
      (*warm)->nrows = 0;
      (*warm)->cstat = NULL;
      (*warm)->rstat = NULL;
      (*warm)->ncols = lpx_get_num_cols(lp->lp);
      if ((*warm)->ncols == 0)
      {  print("CClp_get_warmstart: no columns in LP");
         CClp_free_warmstart(warm);
         return 1;
      }
      (*warm)->nrows = lpx_get_num_rows(lp->lp);
      if ((*warm)->nrows == 0)
      {  print("CClp_get_warmstart: no rows in LP");
         CClp_free_warmstart(warm);
         return 1;
      }
      (*warm)->cstat = ucalloc((*warm)->ncols, sizeof(int));
      (*warm)->rstat = ucalloc((*warm)->nrows, sizeof(int));
      for (i = 1; i <= (*warm)->nrows; i++)
      {  lpx_get_row_info(lp->lp, i, &tagx, NULL, NULL);
         switch (tagx)
         {  case LPX_BS:
               (*warm)->rstat[i-1] = IS_BASIC; break;
            case LPX_NL:
               (*warm)->rstat[i-1] = AT_LOWER; break;
            case LPX_NU:
               (*warm)->rstat[i-1] = AT_UPPER; break;
            case LPX_NS:
               (*warm)->rstat[i-1] = AT_LOWER; break;
            default: insist(tagx != tagx);
         }
      }
      for (j = 1; j <= (*warm)->ncols; j++)
      {  lpx_get_col_info(lp->lp, j, &tagx, NULL, NULL);
         switch (tagx)
         {  case LPX_BS:
               (*warm)->cstat[j-1] = IS_BASIC; break;
            case LPX_NL:
               (*warm)->cstat[j-1] = AT_LOWER; break;
            case LPX_NU:
               (*warm)->cstat[j-1] = AT_UPPER; break;
            case LPX_NS:
               (*warm)->cstat[j-1] = AT_LOWER; break;
            default:
               insist(tagx != tagx);
         }
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_load_warmstart(CClp *lp, CClp_warmstart *warm)
{     /* RESTORES the warmstart information in warm. */
      int i, j, m, n, tagx, type, *rstat, *cstat;
      m = lpx_get_num_rows(lp->lp);
      n = lpx_get_num_cols(lp->lp);
      cstat = warm->cstat;
      rstat = warm->rstat;
      if (cstat == NULL || rstat == NULL)
      {  print("CClp_load_warmstart: no basis information");
         return 0;
      }
      for (j = 0; j < n; j++)
      {  if (cstat[j] == IS_BASIC)
            tagx = LPX_BS;
         else
         {  lpx_get_col_bnds(lp->lp, j+1, &type, NULL, NULL);
            switch (type)
            {  case LPX_FR:
                  tagx = LPX_NF; break;
               case LPX_LO:
                  tagx = LPX_NL; break;
               case LPX_UP:
                  tagx = LPX_NU; break;
               case LPX_DB:
                  tagx = (cstat[j] == AT_UPPER ? LPX_NU : LPX_NL);
                  break;
               case LPX_FX:
                  tagx = LPX_NS; break;
               default:
                  insist(type != type);
            }
         }
         lpx_set_col_stat(lp->lp, j+1, tagx);
      }
      for (i = 0; i < m; i++)
      {  if (rstat[i] == IS_BASIC)
            tagx = LPX_BS;
         else
         {  lpx_get_row_bnds(lp->lp, i+1, &type, NULL, NULL);
            switch (type)
            {  case LPX_FR:
                  tagx = LPX_NF; break;
               case LPX_LO:
                  tagx = LPX_NL; break;
               case LPX_UP:
                  tagx = LPX_NU; break;
               case LPX_DB:
                  tagx = (rstat[i] == AT_UPPER ? LPX_NU : LPX_NL);
                  break;
               case LPX_FX:
                  tagx = LPX_NS; break;
               default:
                  insist(type != type);
            }
         }
         lpx_set_row_stat(lp->lp, i+1, tagx);
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_build_warmstart(CClp_warmstart **warm, CClp_info *info)
{     /* BUILDS some warmstart information from the row/column
         information in info. */
      int i, j;
      CClp_free_warmstart(warm);
      (*warm) = umalloc(sizeof(CClp_warmstart));
      (*warm)->ncols = 0;
      (*warm)->nrows = 0;
      (*warm)->cstat = NULL;
      (*warm)->rstat = NULL;
      (*warm)->ncols = info->ncols;
      insist((*warm)->ncols > 0);
      (*warm)->nrows = info->nrows;
      insist((*warm)->nrows > 0);
      (*warm)->cstat = ucalloc((*warm)->ncols, sizeof(int));
      (*warm)->rstat = ucalloc((*warm)->nrows, sizeof(int));
      for (j = 0; j < (*warm)->ncols; j++)
         (*warm)->cstat[j] = info->cstat[j];
      for (i = 0; i < (*warm)->nrows; i++)
         (*warm)->rstat[i] = info->rstat[i];
      return 0;
}

/*--------------------------------------------------------------------*/

void CClp_free_warmstart(CClp_warmstart **warm)
{     /* FREES the memory used by warm. */
      if (*warm != NULL)
      {  if ((*warm)->cstat != NULL)
            ufree((*warm)->cstat), (*warm)->cstat = NULL;
         if ((*warm)->rstat != NULL)
            ufree((*warm)->rstat), (*warm)->rstat = NULL;
         ufree((*warm)), (*warm) = NULL;
      }
      return;
}

/*--------------------------------------------------------------------*/

int CClp_sread_warmstart(CC_SFILE *file, CClp_warmstart **warm)
{     /* READS warmstart information from the file. */
      int i, j, ccount, rcount;
      char name[5];
      CClp_free_warmstart(warm);
      for (i = 0; i < 4; i++)
         if (CCutil_sread_char(file, &name[i])) goto fail;
      name[4] = '\0';
      if (strncmp(name, SOLVER, 4))
      {  print("CClp_sread_warmstart: warmstart for another solver (%s)"
            " ignored", name);
         return 0;
      }
      if (CCutil_sread_int(file, &ccount)) goto fail;
      if (CCutil_sread_int(file, &rcount)) goto fail;
      (*warm) = umalloc(sizeof(CClp_warmstart));
      (*warm)->ncols = 0;
      (*warm)->nrows = 0;
      (*warm)->cstat = NULL;
      (*warm)->rstat = NULL;
      (*warm)->cstat = ucalloc(ccount, sizeof(int));
      (*warm)->rstat = ucalloc(rcount, sizeof(int));
      for (j = 0; j < ccount; j++)
         if (CCutil_sread_bits(file, &(((*warm)->cstat)[j]), 2))
            goto fail;
      for (i = 0; i < rcount; i++)
         if (CCutil_sread_bits(file, &(((*warm)->rstat)[i]), 1))
            goto fail;
      (*warm)->ncols = ccount;
      (*warm)->nrows = rcount;
      return 0;
fail: CClp_free_warmstart(warm);
      return 1;
}

/*--------------------------------------------------------------------*/

int CClp_swrite_warmstart(CC_SFILE *file, CClp_warmstart *warm)
{     /* WRITES warmstart information from the f. */
      int i, j;
      const char *name = SOLVER;
      for (i = 0; i < 4; i++)
         if (CCutil_swrite_char(file, name[i])) return 1;
      if (CCutil_swrite_int(file, warm->ncols)) return 1;
      if (CCutil_swrite_int(file, warm->nrows)) return 1;
      for (j = 0; j < warm->ncols; j++)
         if (CCutil_swrite_bits (file, warm->cstat[j], 2)) return 1;
      for (i = 0; i < warm->nrows; i++)
         if (CCutil_swrite_bits (file, warm->rstat[i], 1)) return 1;
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_get_info(CClp *lp, CClp_info **info)
{     /* BUILDS information useful for efficiently answering questions
         about the status of rows and columns. */
      int i, j, tagx;
      CClp_free_info(info);
      (*info) = umalloc(sizeof(CClp_info));
      (*info)->ncols = 0;
      (*info)->nrows = 0;
      (*info)->cstat = NULL;
      (*info)->rstat = NULL;
      (*info)->ncols = lpx_get_num_cols(lp->lp);
      insist((*info)->ncols > 0);
      (*info)->nrows = lpx_get_num_rows(lp->lp);
      insist((*info)->nrows > 0);
      (*info)->cstat = ucalloc((*info)->ncols, sizeof(int));
      (*info)->rstat = ucalloc((*info)->nrows, sizeof(int));
      for (i = 1; i <= (*info)->nrows; i++)
      {  lpx_get_row_info(lp->lp, i, &tagx, NULL, NULL);
         (*info)->rstat[i-1] = tagx == LPX_BS ? 1 : 0;
      }
      for (j = 1; j <= (*info)->ncols; j++)
      {  lpx_get_col_info(lp->lp, j, &tagx, NULL, NULL);
         switch (tagx)
         {  case LPX_BS:
               (*info)->cstat[j-1] = 0; break;
            case LPX_NL:
               (*info)->cstat[j-1] = 1; break;
            case LPX_NU:
               (*info)->cstat[j-1] = 2; break;
            case LPX_NS:
               (*info)->cstat[j-1] = 1; break;
            default:
               insist(tagx != tagx);
         }
      }
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_create_info(CClp_info **info, int nrows, int ncols)
{     /* CREATES a structure for storing information about the status
         of rows and columns. */
      int i, j;
      CClp_free_info(info);
      (*info) = umalloc(sizeof(CClp_info));
      (*info)->ncols = 0;
      (*info)->nrows = 0;
      (*info)->cstat = NULL;
      (*info)->rstat = NULL;
      (*info)->ncols = ncols;
      insist(ncols > 0);
      (*info)->nrows = nrows;
      insist(nrows > 0);
      (*info)->cstat = ucalloc((*info)->ncols, sizeof(int));
      (*info)->rstat = ucalloc((*info)->nrows, sizeof(int));
      for (j = 0; j < ncols; j++) (*info)->cstat[j] = 0;
      for (i = 0; i < nrows; i++) (*info)->rstat[i] = 0;
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_is_col_active(CClp_info *info, int j)
{     /* returns 1 if column j is active, 0 otherwise.
         "active" means participating in the current solution (for
         example, it could mean basic or nonbasic at upper bound) */
      if (j < 0 || j >= info->ncols) return 0;
      return info->cstat[j] == 1 || info->cstat[j] == 2;
}

/*--------------------------------------------------------------------*/

int CClp_is_row_active(CClp_info *info, int i)
{     /* returns 1 if row i is active, 0 otherwise.
         "active" means participating in the current solution (for
         example, it could mean the row's slack is non-basic) */
      if (i < 0 || i >= info->nrows) return 0;
      return info->rstat[i] == 0;
}

/*--------------------------------------------------------------------*/

void CClp_set_col_active(CClp_info *info, int j)
{     /* marks column j as active */
      if (j >= 0 && j < info->ncols) info->cstat[j] = 1;
      return;
}

/*--------------------------------------------------------------------*/

void CClp_set_col_inactive(CClp_info *info, int j)
{     /*  marks column j as inactive */
      if (j >= 0 && j < info->ncols) info->cstat[j] = 0;
      return;
}

/*--------------------------------------------------------------------*/

void CClp_set_col_upper(CClp_info *info, int j)
{     /* marks column j as active at upper bound */
      if (j >= 0 && j < info->ncols) info->cstat[j] = 2;
      return;
}

/*--------------------------------------------------------------------*/

void CClp_set_row_active(CClp_info *info, int i)
{     /* marks row i as active */
      if (i >= 0 && i < info->nrows) info->rstat[i] = 0;
      return;
}

/*--------------------------------------------------------------------*/

void CClp_set_row_inactive(CClp_info *info, int i)
{     /* marks row i as inactive */
      if (i >= 0 && i < info->nrows) info->rstat[i] = 1;
}

/*--------------------------------------------------------------------*/

void CClp_free_info(CClp_info **info)
{     /* FREES the memory used by info. */
      if ((*info) != NULL)
      {  if ((*info)->cstat != NULL)
            ufree((*info)->cstat), (*info)->cstat = NULL;
         if ((*info)->rstat != NULL)
            ufree((*info)->rstat), (*info)->rstat = NULL;
         ufree((*info)), (*info) = NULL;
      }
      return;
}

/*--------------------------------------------------------------------*/

int CClp_x(CClp *lp, double *x)
{     /* RETURNS the current LP solution.
         - x should be an array of length at least ncols. */
      int ncols, j;
      ncols = lpx_get_num_cols(lp->lp);
      insist(ncols > 0);
      for (j = 0; j < ncols; j++)
         lpx_get_col_info(lp->lp, j+1, NULL, &x[j], NULL);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_rc(CClp *lp, double *rc)
{     /* RETURNS the current reduced costs.
         - rc should be an array of length at least ncols. */
      int ncols, j;
      ncols = lpx_get_num_cols(lp->lp);
      insist(ncols > 0);
      for (j = 0; j < ncols; j++)
         lpx_get_col_info(lp->lp, j+1, NULL, NULL, &rc[j]);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_pi(CClp *lp, double *pi)
{     /* RETURNS the dual values on the constraints.
         - pi should be an array of length at least nrows. */
      int nrows, i;
      nrows = lpx_get_num_rows(lp->lp);
      insist(nrows > 0);
      for (i = 0; i < nrows; i++)
         lpx_get_row_info(lp->lp, i+1, NULL, NULL, &pi[i]);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_objval(CClp *lp, double *obj)
{     /* RETURNS the objective value of the lp. */
      if (obj != NULL) *obj = lpx_get_obj_val(lp->lp);
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_nrows(CClp *lp)
{     /* RETURNS the number of rows in the LP. */
      return lpx_get_num_rows(lp->lp);
}

/*--------------------------------------------------------------------*/

int CClp_ncols(CClp *lp)
{     /* RETURNS the number of columns in the LP. */
      return lpx_get_num_cols(lp->lp);
}

/*--------------------------------------------------------------------*/

int CClp_nnonzeros(CClp *lp)
{     /* RETURNS the number of nonzeros in the LP. */
      return lpx_get_num_nz(lp->lp);
}

/*--------------------------------------------------------------------*/

int CClp_status(CClp *lp, int *status)
{     /* CHECKS whether the current lp is infeasible or whether an
         optimal solution has been found. It returns an error if the LP
         has not not been optimized.
         - lp is the lp;
         - status returns 0 if the lp has an optimal solution and 1 if
           it is infeasible. */
      insist(lp == lp);
      insist(status == status);
      fault("CClp_status: not implemented");
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_getweight(CClp *lp, int nrows, int *rmatbeg, int *rmatind,
      double *rmatval, double *weight)
{     /* COMPUTES the duals of the steepest edge norms for the n rows
         specified in rmatbeg, rmatind, and rmatval.
         - weight returns the array of weights; the array should be at
           least nrows long. */
      insist(lp == lp);
      insist(nrows == nrows);
      insist(rmatbeg == rmatbeg);
      insist(rmatind == rmatind);
      insist(rmatval == rmatval);
      insist(weight == weight);
      fault("CClp_getweight: not implemented");
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_dump_lp(CClp *lp, const char *fname)
{     /* WRITES the LP to file fname. */
      insist(lp == lp);
      insist(fname == fname);
      fault("CClp_dump_lp: not implemented");
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_getgoodlist(CClp *lp, int *goodlist, int *goodlen_p,
      double *downpen, double *uppen)
{     /* RETURNS an array of the column indices corresponding to
         variables that move in both the up and down directions. This
         is a useful list of candidates for strong branching.
         - goodlist, downpen and uppen should be arrays of length at
           least ncols. */
      insist(lp == lp);
      insist(goodlist == goodlist);
      insist(downpen == downpen);
      insist(uppen == uppen);
      /* not implemented */
      *goodlen_p = 0;
      return 0;
}

/*--------------------------------------------------------------------*/

int CClp_strongbranch(CClp *lp, int *candidatelist, int ncand,
      double *downpen, double *uppen, int iterations, double upperbound)
{     /* RETURNS estimates of the lp values obtained by setting each of
         the ncand variables listed in candidatelist to 0 and 1. The
         estimates are obtained by performing iterations pivots of dual
         simplex method. upperbound is used to cutoff the dual simplex
         method. downpen and uppen should never be > upperbound.
         -downpen and uppen should be arrays of length at least ncand */
      int i;
      insist(lp == lp);
      insist(candidatelist == candidatelist);
      insist(iterations == iterations);
      /* not implemented */
      for (i = 0; i < ncand; i++)
      {  downpen[i] = upperbound;
         uppen[i] = upperbound;
      }
      return 0;
}

/* eof */

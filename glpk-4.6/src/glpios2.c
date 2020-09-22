/* glpios2.c */

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

#include <stddef.h>
#include "glpios.h"
#include "glplib.h"
#include "glplpx.h"

/**********************************************************************/
/* * *                LP SOLVER INTERFACE ROUTINES                * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_extract_lp - extract LP relaxation of current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- LPX *ios_extract_lp(IOS *ios);
--
-- *Description*
--
-- The routine ios_extract_lp constructs an LP problem object which is
-- LP relaxation of the current subproblem.
--
-- Ordinal numbers of rows and columns in the LP problem object are the
-- same as in the current subproblem.
--
-- The LP problem object created by this routine is a separate program
-- object not connected with the integer optimization suite. Therefore
-- any changes applied to it do not affect the current subproblem.
--
-- *Returns*
--
-- The routine ios_extract_lp returns a pointer to the constructed LPX
-- object (for technical reasons the pointer is casted to void *). */

struct load_info
{     /* transitional information passed to the routine next_aij */
      IOS *ios;
      /* pointer to the main data structure */
      int m;
      /* number of rows */
      int n;
      /* number of columns */
      int j;
      /* number of the current column being read */
      int len;
      /* number of constraint coefficients in j-th column which are not
         read yet */
      int *ind; /* int ind[1+m]; */
      /* row indices of constraint coefficients */
      double *val; /* double val[1+m]; */
      /* numeric values of constraint coefficients */
};

static double next_aij(void *_info, int *i, int *j)
{     /* "read" a next element of the constraint matrix */
      struct load_info *info = _info;
      double val;
      while (info->len == 0)
      {  /* either it is the first call or the current column has been
            completely read; choose the first/next column */
         if (info->j < info->n)
            info->j++;
         else
         {  /* the entire matrix has been completely read */
            info->j = 0;
            *i = *j = 0;
            return 0.0;
         }
         /* obtain j-th column of the constraint matrix */
         info->len = ios_get_mat_col(info->ios, info->j, info->ind,
            info->val);
      }
      /* "read" the next element from the current column */
      insist(1 <= info->len && info->len <= info->m);
      *i = info->ind[info->len];
      *j = info->j;
      val = info->val[info->len];
      info->len--;
      return val;
}

void *ios_extract_lp(IOS *ios)
{     LPX *lp;
      struct load_info info;
      int m, n, i, j, type, stat, dir;
      double lb, ub, coef;
      char *name;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_extract_lp: current subproblem does not exist");
      /* create LP problem object */
      lp = lpx_create_prob();
      /* set optimization direction */
      dir = ios->dir;
      switch (dir)
      {  case IOS_MIN: dir = LPX_MIN; break;
         case IOS_MAX: dir = LPX_MAX; break;
         default: insist(dir != dir);
      }
      lpx_set_obj_dir(lp, dir);
      /* create rows */
      m = ios_get_num_rows(ios);
      if (m > 0) lpx_add_rows(lp, m);
      for (i = 1; i <= m; i++)
      {  /* set row name */
         name = ios_get_row_name(ios, i);
         lpx_set_row_name(lp, i, name);
         /* set row type and bounds */
         type = ios_get_row_bnds(ios, i, &lb, &ub);
         switch (type)
         {  case IOS_FR: type = LPX_FR; break;
            case IOS_LO: type = LPX_LO; break;
            case IOS_UP: type = LPX_UP; break;
            case IOS_DB: type = LPX_DB; break;
            case IOS_FX: type = LPX_FX; break;
            default: insist(type != type);
         }
         lpx_set_row_bnds(lp, i, type, lb, ub);
         /* set row status */
         stat = ios_get_row_soln(ios, i, NULL, NULL);
         switch (stat)
         {  case IOS_BS: stat = LPX_BS; break;
            case IOS_NL: stat = LPX_NL; break;
            case IOS_NU: stat = LPX_NU; break;
            case IOS_NF: stat = LPX_NF; break;
            case IOS_NS: stat = LPX_NS; break;
            default: insist(stat != stat);
         }
         lpx_set_row_stat(lp, i, stat);
      }
      /* create columns */
      n = ios_get_num_cols(ios);
      if (n > 0) lpx_add_cols(lp, n);
      for (j = 1; j <= n; j++)
      {  /* set column name */
         name = ios_get_col_name(ios, j);
         lpx_set_col_name(lp, j, name);
         /* set column type and bounds */
         type = ios_get_col_bnds(ios, j, &lb, &ub);
         switch (type)
         {  case IOS_FR: type = LPX_FR; break;
            case IOS_LO: type = LPX_LO; break;
            case IOS_UP: type = LPX_UP; break;
            case IOS_DB: type = LPX_DB; break;
            case IOS_FX: type = LPX_FX; break;
            default: insist(type != type);
         }
         lpx_set_col_bnds(lp, j, type, lb, ub);
         /* set column objective coefficient */
         coef = ios_get_obj_coef(ios, j);
         lpx_set_col_coef(lp, j, coef);
         /* set column status */
         stat = ios_get_col_soln(ios, j, NULL, NULL);
         switch (stat)
         {  case IOS_BS: stat = LPX_BS; break;
            case IOS_NL: stat = LPX_NL; break;
            case IOS_NU: stat = LPX_NU; break;
            case IOS_NF: stat = LPX_NF; break;
            case IOS_NS: stat = LPX_NS; break;
            default: insist(stat != stat);
         }
         lpx_set_col_stat(lp, j, stat);
      }
      /* set constant term of the objective function */
      coef = ios_get_obj_coef(ios, 0);
      lpx_set_obj_c0(lp, coef);
      /* load the constraint matrix */
      info.ios = ios;
      info.m = m;
      info.n = n;
      info.j = 0;
      info.len = 0;
      info.ind = ucalloc(1+m, sizeof(int));
      info.val = ucalloc(1+m, sizeof(double));
      lpx_load_mat(lp, &info, next_aij);
      insist(ios_get_num_nz(ios) == lpx_get_num_nz(lp));
      ufree(info.ind);
      ufree(info.val);
      /* LP relaxation has been extracted */
      return lp;
}

/*----------------------------------------------------------------------
-- eval_lp_obj - compute objective value for LP relaxation.
--
-- This routine computes the value of the objective function for basic
-- solution of LP relaxaton using the formula:
--
--    Z = c[1]*x[1] + ... + c[n]*x[n] + c[0],
--
-- where x[1], ..., x[n] are structural variables, c[1], ..., c[n] are
-- objective coefficients, c[0] is the constant term. */

static double eval_lp_obj(IOS *ios)
{     int n, j;
      double x, obj;
      n = ios_get_num_cols(ios);
      obj = ios_get_obj_coef(ios, 0);
      for (j = 1; j <= n; j++)
      {  ios_get_col_soln(ios, j, &x, NULL);
         obj += ios_get_obj_coef(ios, j) * x;
      }
      return obj;
}

/*----------------------------------------------------------------------
-- eval_lp_sum - compute sum of infeasibilities for LP relaxation.
--
-- This routine computes the sum of primal infeasibilities for basic
-- solution of LP relaxation using the formula:
--
--        m+n ( 0,            if l[k] <= x[k] <= u[k]
--    S = Sum ( l[k] - x[k],  if x[k] < l[k]
--        k=1 ( x[k] - u[k],  if x[k] > u[k]
--
-- where x[1], ..., x[m] are auxiliary variables, x[1], ..., x[n] are
-- structural variables, l[1], ..., l[m+n] are lower bounds, u[1], ...,
-- u[m+n] are upper bounds. */

static double eval_lp_sum(IOS *ios)
{     int m, n, i, j, k, type;
      double lb, ub, x, sum;
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      sum = 0.0;
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  i = k;
            type = ios_get_row_bnds(ios, i, &lb, &ub);
            ios_get_row_soln(ios, i, &x, NULL);
         }
         else
         {  j = k - m;
            type = ios_get_col_bnds(ios, j, &lb, &ub);
            ios_get_col_soln(ios, j, &x, NULL);
         }
         switch (type)
         {  case IOS_FR:
               break;
            case IOS_LO:
               if (x < lb) sum += lb - x;
               break;
            case IOS_UP:
               if (x > ub) sum += x - ub;
               break;
            case IOS_DB:
            case IOS_FX:
               if (x < lb) sum += lb - x;
               if (x > ub) sum += x - ub;
               break;
            default:
               insist(type != type);
         }
      }
      return sum;
}

/*----------------------------------------------------------------------
-- ios_put_lp_soln - store basic solution of LP relaxation.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_put_lp_soln(IOS *ios, LPX *lp);
--
-- *Description*
--
-- The routine ios_put_lp_soln stores basic solution of LP relaxation
-- of the current subproblem obtained by the LP solver and contained in
-- the LP problem object, which the parameter lp points to, back to the
-- current subproblem (for some technical reasons a pointer to the LP
-- problem object is casted to void *). */

void ios_put_lp_soln(IOS *ios, void *_lp)
{     LPX *lp = _lp;
      int m, n, i, j, stat;
      double prim, dual;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_put_lp_soln: current subproblem does not exist");
      /* determine number of rows and columns */
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      if (!(lpx_get_num_rows(lp) == m && lpx_get_num_cols(lp) == n))
         fault("ios_put_lp_soln: lp problem object seems to have been c"
            "hanged");
      /* store primal status */
      stat = lpx_get_prim_stat(lp);
      switch (stat)
      {  case LPX_P_UNDEF:  stat = IOS_UNDEF;  break;
         case LPX_P_FEAS:   stat = IOS_FEAS;   break;
         case LPX_P_INFEAS: stat = IOS_INFEAS; break;
         case LPX_P_NOFEAS: stat = IOS_NOFEAS; break;
         default: insist(stat != stat);
      }
      ios->p_stat = stat;
      /* store dual status */
      stat = lpx_get_dual_stat(lp);
      switch (stat)
      {  case LPX_D_UNDEF:  stat = IOS_UNDEF;  break;
         case LPX_D_FEAS:   stat = IOS_FEAS;   break;
         case LPX_D_INFEAS: stat = IOS_INFEAS; break;
         case LPX_D_NOFEAS: stat = IOS_NOFEAS; break;
         default: insist(stat != stat);
      }
      ios->d_stat = stat;
      /* store basic solution components for rows */
      for (i = 1; i <= m; i++)
      {  IOSROW *row;
         lpx_get_row_info(lp, i, &stat, &prim, &dual);
         switch (stat)
         {  case LPX_BS: stat = IOS_BS; break;
            case LPX_NL: stat = IOS_NL; break;
            case LPX_NU: stat = IOS_NU; break;
            case LPX_NF: stat = IOS_NF; break;
            case LPX_NS: stat = IOS_NS; break;
            default: insist(stat != stat);
         }
         ios_set_row_stat(ios, i, stat);
         row = ios_get_row_ptr(ios, i);
         row->prim = prim;
         row->dual = dual;
      }
      /* store basic solution components for columns */
      for (j = 1; j <= n; j++)
      {  IOSCOL *col;
         lpx_get_col_info(lp, j, &stat, &prim, &dual);
         switch (stat)
         {  case LPX_BS: stat = IOS_BS; break;
            case LPX_NL: stat = IOS_NL; break;
            case LPX_NU: stat = IOS_NU; break;
            case LPX_NF: stat = IOS_NF; break;
            case LPX_NS: stat = IOS_NS; break;
            default: insist(stat != stat);
         }
         ios_set_col_stat(ios, j, stat);
         col = ios_get_col_ptr(ios, j);
         col->prim = prim;
         col->dual = dual;
      }
      /* store the objective function value */
      ios->lp_obj = eval_lp_obj(ios);
      /* store the sum of primal infeasibilities */
      ios->lp_sum = eval_lp_sum(ios);
      return;
}

/*----------------------------------------------------------------------
-- ios_solve_root - solve initial LP relaxation.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_solve_root(IOS *ios);
--
-- *Description*
--
-- The routine ios_solve_root solves initial LP relaxation of the root
-- subproblem using the primal simplex method.
--
-- *Returns*
--
-- If LP relaxation has been solved successfully, the routine returns
-- zero. Otherwise the routine returns non-zero. */

int ios_solve_root(IOS *ios)
{     LPX *lp;
      int ret;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_solve_root: current subproblem does not exist");
      if (ios_get_curr_node(ios) != 1)
         fault("ios_solve_root: current subproblem is not the root subp"
            "roblem");
      /* extract LP relaxation of the root subproblem */
      lp = ios_extract_lp(ios);
      /* set some control parameters */
      lpx_set_int_parm(lp, LPX_K_DUAL, 0);
      lpx_set_int_parm(lp, LPX_K_ITCNT, ios->it_cnt);
      /* scale the LP relaxation, if required */
      if (ios->scale) lpx_scale_prob(lp);
      /* construct an initial basis to start from */
      lpx_set_int_parm(lp,
         LPX_K_MSGLEV, ios->msg_lev <= 1 ? ios->msg_lev : 3);
      switch (ios->init_lp)
      {  case 0:
            lpx_std_basis(lp); break;
         case 1:
            lpx_adv_basis(lp); break;
         case 2:
            break;
         default:
            insist(ios != ios);
      }
      /* try to solve the LP relaxation */
      lpx_set_int_parm(lp,
         LPX_K_MSGLEV, ios->msg_lev <= 1 ? ios->msg_lev : 2);
      ret = lpx_simplex(lp);
      if (ret == LPX_E_OK) ret = 0; else ret = 1;
      /* store basic solution back to the current subproblem */
      ios_put_lp_soln(ios, lp);
      /* and also copy back the simplex iteration count */
      ios->it_cnt = lpx_get_int_parm(lp, LPX_K_ITCNT);
      /* delete the LP relaxation */
      lpx_delete_prob(lp);
      return ret;
}

/*----------------------------------------------------------------------
-- eval_pi_opt - compute multipliers (optimal solution).
--
-- This routine computes Lagrange multipliers which are needed for the
-- application procedure for generating columns. It is assumed that the
-- current LP relaxation has optimal solution.
--
-- The computations are based on the main formula:
--
--    d[k] = c[k] - pi' * A~[k],                                     (1)
--
-- where d[k] is the reduced cost of (structural or auxiliary) variable
-- x[k], c[k] is the objective coefficient at x[k], pi is the vector of
-- Lagrange multipliers, A~[k] is the column of the expanded constraint
-- matrix A~ = (I | -A) for x[k].
--
-- If x[i] is an auxiliary variable, c[i] is zero, and A~[i] is a unity
-- column. Thus:
--
--    pi[i] = - d[i],                                                (2)
--
-- where d[i] is the reduced cost (dual value) of x[i] already computed
-- for the original objective function. */

static void eval_pi_opt(IOS *ios, LPX *lp)
{     int m, i;
      double d;
      /* the basic solution must be optimal */
      insist(ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS);
      /* compute Lagrange multipliers for all rows */
      m = ios_get_num_rows(ios);
      for (i = 1; i <= m; i++)
      {  /* obtain dual value for i-th row */
         lpx_get_row_info(lp, i, NULL, NULL, &d);
         /* compute pi[i] using formula (2) */
         ios_get_row_ptr(ios, i)->pi = 0.0 - d;
      }
      return;
}

/*----------------------------------------------------------------------
-- eval_pi_nfs - compute multipliers (no feasible solution).
--
-- This routine computes Lagrange multipliers which are needed for the
-- application procedure for generating columns. It is assumed that the
-- current LP relaxation has no primal feasible solution.
--
-- The computations are based on the main formula:
--
--    d[k] = c[k] - pi' * A~[k],                                     (1)
--
-- where d[k] is the reduced cost of (structural or auxiliary) variable
-- x[k], c[k] is the objective coefficient at x[k], pi is the vector of
-- Lagrange multipliers, A~[k] is the column of the expanded constraint
-- matrix A~ = (I | -A) for x[k].
--
-- Objective coefficients c[k] used in this routine correspond to the
-- auxiliary objective which is the sum of primal infeasibilities to be
-- minimized:
--
--              0, if lb[k] <= x[k] <= ub[k]
--    c[k] =   -1, if x[k] < lb[k]                                   (2)
--             +1, if x[k] > ub[k]
--
-- where lb[k] are ub[k] are, respectively, the lower and upper bounds
-- of (structural or auxiliary) variable x[k].
--
-- NOTE: If the original objective sense is maximization, the auxiliary
-- objective is taken with the minus sign in order to keep the original
-- sense, i.e. all c[k] in (2) have opposite signs.
--
-- If x[i] is an auxiliary variable, A~[i] is a unity column. Thus:
--
--    pi[i] = c[i] - d[i],                                           (3)
--
-- where c[i] and d[i] are, respectively, the objective coefficient (2)
-- and the reduced cost (dual value) of x[i] computed for the auxiliary
-- objective function. */

static void eval_pi_nfs(IOS *ios, LPX *lp)
{     int m, n, i, j, k, type, stat;
      double one, lb, ub, x, *c, d;
      /* there must be no primal feasible solution */
      insist(ios->p_stat == IOS_NOFEAS);
      /* if the original objective sense is maximization, the auxiliary
         objective function should be taken with the minus sign */
      switch (ios->dir)
      {  case IOS_MIN: one = +1.0; break;
         case IOS_MAX: one = -1.0; break;
         default: insist(lp != lp);
      }
      /* compute coefficients of the auxiliary objective function using
         formula (2) */
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      c = ucalloc(1+m+n, sizeof(double));
      for (k = 1; k <= m+n; k++)
      {  if (k <= m)
         {  i = k;
            /* obtain type and bounds of i-th row */
            type = ios_get_row_bnds(ios, i, &lb, &ub);
            /* obtain status and primal value of i-th row */
            stat = ios_get_row_soln(ios, i, &x, NULL);
         }
         else
         {  j = k - m;
            /* obtain type and bounds of j-th column */
            type = ios_get_col_bnds(ios, j, &lb, &ub);
            /* obtain status and primal value of j-th column */
            stat = ios_get_col_soln(ios, j, &x, NULL);
         }
         /* only bounds of basic variables can be violated */
         c[k] = 0.0;
         if (stat == IOS_BS)
         {  switch (type)
            {  case IOS_FR:
                  break;
               case IOS_LO:
                  if (x < lb) c[k] = - one;
                  break;
               case IOS_UP:
                  if (x > ub) c[k] = + one;
                  break;
               case IOS_DB:
               case IOS_FX:
                  if (x < lb) c[k] = - one;
                  if (x > ub) c[k] = + one;
                  break;
               default:
                  insist(type != type);
            }
         }
      }
      /* substitute auxiliary variables into the auxiliary objective
         function to express the latter through structural variables;
         objective coefficients at auxiliary variables are kept in
         c[1], ..., c[m] while new objective coefficients at structural
         variables are accumulated in c[m+1], ..., c[m+n] */
      {  int *ind, len;
         double *val;
         ind = ucalloc(1+n, sizeof(int));
         val = ucalloc(1+n, sizeof(double));
         for (i = 1; i <= m; i++)
         {  if (c[i] == 0.0) continue;
            len = lpx_get_mat_row(lp, i, ind, val);
            for (k = 1; k <= len; k++)
            {  j = ind[k];
               insist(1 <= j && j <= n);
               c[m+j] += c[i] * val[k];
            }
         }
         ufree(ind);
         ufree(val);
      }
      /* set coefficients of the auxiliary objective function */
      lpx_set_obj_coef(lp, 0, 0.0);
      for (j = 1; j <= n; j++) lpx_set_obj_coef(lp, j, c[m+j]);
      /* compute components of basic solution to the LP relaxation for
         the auxiliary objective function */
      insist(lpx_warm_up(lp) == LPX_E_OK);
      /* compute Lagrange multipliers for all rows */
      for (i = 1; i <= m; i++)
      {  /* obtain dual value for i-th row */
         lpx_get_row_info(lp, i, NULL, NULL, &d);
         /* compute pi[i] using formula (3) */
         ios_get_row_ptr(ios, i)->pi = c[i] - d;
      }
      ufree(c);
      return;
}

/*----------------------------------------------------------------------
-- ios_solve_node - solve LP relaxation of current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_solve_node(IOS *ios);
--
-- *Description*
--
-- The routine ios_solve_node solves or re-optimizes LP relaxation of
-- the current subproblem using the primal or dual simplex method.
--
-- If the column generation is enabled, the routine also computes the
-- Lagrange multipliers for all rows.
--
-- *Returns*
--
-- If LP relaxation has been solved successfully, the routine returns
-- zero. Otherwise the routine returns non-zero. */

int ios_solve_node(IOS *ios)
{     LPX *lp;
      int ret;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_solve_node: current subproblem does not exist");
      /* extract LP relaxation of the current subproblem */
      lp = ios_extract_lp(ios);
      /* set some control parameters */
      lpx_set_int_parm(lp, LPX_K_DUAL, 1);
      lpx_set_int_parm(lp, LPX_K_ITCNT, ios->it_cnt);
      /* scale the LP relaxation, if required */
      if (ios->scale) lpx_scale_prob(lp);
      /* if the incumbent objective value is already known, use it to
         prematurely terminate the dual simplex search (this technique
         MUST NOT be used if the column generation is ENABLED) */
      if (!ios->col_gen && ios->found)
      {  switch (ios->dir)
         {  case IOS_MIN:
               lpx_set_real_parm(lp, LPX_K_OBJUL, ios->best);
               break;
            case IOS_MAX:
               lpx_set_real_parm(lp, LPX_K_OBJLL, ios->best);
               break;
            default:
               insist(ios != ios);
         }
      }
      /* try to solve/re-optimize the LP relaxation */
      lpx_set_int_parm(lp,
         LPX_K_MSGLEV, ios->msg_lev <= 1 ? ios->msg_lev : 2);
      lpx_set_real_parm(lp,
         LPX_K_OUTDLY, ios->msg_lev <= 2 ? ios->out_dly : 0.0);
      ret = lpx_simplex(lp);
      if (ret == LPX_E_OK || ret == LPX_E_OBJLL || ret == LPX_E_OBJUL)
         ret = 0;
      else
         ret = 1;
      /* store basic solution back to the current subproblem */
      ios_put_lp_soln(ios, lp);
      /* and also copy back the simplex iteration count */
      ios->it_cnt = lpx_get_int_parm(lp, LPX_K_ITCNT);
      /* if the column generation is enabled and the basic solution is
         either optimal or proven to be primal infeasible, compute the
         Lagrange multipliers for all rows */
      if (ios->col_gen && ret == 0)
      {  switch (lpx_get_status(lp))
         {  case LPX_OPT:
               eval_pi_opt(ios, lp);
               break;
            case LPX_NOFEAS:
               eval_pi_nfs(ios, lp);
               break;
            default:
               break;
         }
      }
      /* delete the LP relaxation */
      lpx_delete_prob(lp);
      return ret;
}

/* eof */

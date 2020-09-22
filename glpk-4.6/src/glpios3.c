/* glpios3.c */

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
#include <math.h>
#include <stdio.h>
#include "glpios.h"
#include "glplib.h"

/**********************************************************************/
/* * *                  IOS FUNCTIONARY ROUTINES                  * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- show_progress - display progress of the search.
--
-- This routine displays some information about progress of the search.
--
-- This information includes:
--
-- the current number of iterations performed by the simplex solver;
--
-- the incumbent objective value, i.e. the objective value for the best
-- known integer feasible solution, which is upper (minimization) or
-- lower (maximization) global bound for optimal solution of the entire
-- problem;
--
-- the objective bound for subproblem corresponding to the pseudo-root,
-- which is lower (minimization) or upper (maximization) global bound
-- for optimal solution of the entire problem;
--
-- the number of open (active) subproblems;
--
-- the number of completely explored subproblems, i.e. whose nodes have
-- been already removed from the tree. */

static void show_progress(IOS *ios)
{     int p, a_cnt, n_cnt, t_cnt;
      char best[50], glob[50], *rho;
      /* format best known integer feasible solution */
      if (ios->found)
         sprintf(best, "%17.9e", ios->best);
      else
         sprintf(best, "%17s", "not found yet");
      /* determine reference number of the pseudo-root subproblem */
      p = ios_pseudo_root(ios);
      /* format objective bound for the pseudo-root subproblem */
      if (p == 0)
         sprintf(glob, "%17s", "tree is empty");
      else
      {  double bound;
         bound = ios_get_npd_ptr(ios, p)->bound;
         if (bound == -DBL_MAX)
            sprintf(glob, "%17s", "-inf");
         else if (bound == +DBL_MAX)
            sprintf(glob, "%17s", "+inf");
         else
            sprintf(glob, "%17.9e", bound);
      }
      /* choose the relation sign between global bounds */
      switch (ios->dir)
      {  case IOS_MIN: rho = ">="; break;
         case IOS_MAX: rho = "<="; break;
         default: insist(ios != ios);
      }
      /* determine current size of the tree */
      iet_get_tree_size(ios->iet, &a_cnt, &n_cnt, &t_cnt);
      /* display progress of the search */
      print("+%6d:   ip_obj = %s %s %17s (%d; %d)", ios->it_cnt, best,
         rho, glob, a_cnt, t_cnt - n_cnt);
      ios->tm_lag = utime();
      return;
}

/*----------------------------------------------------------------------
-- clear_lp_soln - invalidate basic solution of LP relaxation.
--
-- This routine invalidates basic solution of the current LP relaxation
-- in the beginning of the minor loop. */

static void clear_lp_soln(IOS *ios)
{     int m, n, i, j;
      ios->p_stat = IOS_UNDEF;
      ios->d_stat = IOS_UNDEF;
      ios->lp_obj = 0.0;
      ios->lp_sum = 0.0;
      ios->ii_cnt = 0;
      ios->ii_sum = 0.0;
      m = ios_get_num_rows(ios);
      for (i = 1; i <= m; i++)
      {  IOSROW *row;
         row = ios_get_row_ptr(ios, i);
         row->prim = 0.0;
         row->dual = 0.0;
         row->pi = 0.0;
      }
      n = ios_get_num_cols(ios);
      for (j = 1; j <= n; j++)
      {  IOSCOL *col;
         col = ios_get_col_ptr(ios, j);
         col->prim = 0.0;
         col->dual = 0.0;
         col->frac = 0;
      }
      return;
}

/*----------------------------------------------------------------------
-- check_integrality - check integrality of basic solution.
--
-- This routine checks if the basic solution of LP relaxation of the
-- current subproblem satisfies to integrality conditions, i.e. that all
-- variables of integer kind have integral primal values. (The solution
-- is *not* assumed to be primal feasible.)
--
-- For each variable of integer kind the routine computes the following
-- quantity:
--
--    ii(x[j]) = min(x[j] - floor(x[j]), ceil(x[j]) - x[j]),         (1)
--
-- which is a measure of the integer infeasibility (non-integrality) of
-- x[j] (for example, ii(2.1) = 0.1, ii(3.7) = 0.3, ii(5.0) = 0). It is
-- understood that 0 <= ii(x[j]) <= 0.5, and variable x[j] is integer
-- feasible if ii(x[j]) = 0. However, due to floating-point arithmetic
-- the routine checks less restrictive condition:
--
--    ii(x[j]) <= tol_int,                                           (2)
--
-- where tol_int is a given tolerance (small positive number) and marks
-- each variable which does not satisfy to (2) as integer infeasible by
-- setting its fractionality flag.
--
-- In order to characterize integer infeasibility of the basic solution
-- in the whole the routine computes two parameters: ii_cnt, which is
-- the number of variables with the fractionality flag set, and ii_sum,
-- which is the sum of integer infeasibilities (1). */

static void check_integrality(IOS *ios)
{     int j, kind, type, stat;
      double lb, ub, x, temp1, temp2;
      /* walk through the set of columns (structural variables) */
      for (j = 1; j <= ios->iet->n; j++)
      {  /* determine kind of j-th column */
         kind = ios_get_col_kind(ios, j);
         /* if the column is continuous, skip it */
         if (kind == IOS_NUM) continue;
         /* otherwise the column can be only of integer kind */
         insist(kind == IOS_INT);
         /* obtain the status and primal value of the column */
         stat = ios_get_col_soln(ios, j, &x, NULL);
         /* if the column is non-basic, it is integer feasible */
         if (stat != IOS_BS) continue;
         /* obtain the type and bounds of the column */
         type = ios_get_col_bnds(ios, j, &lb, &ub);
         /* if the column's primal value is close to the lower bound,
            the column is integer feasible within given tolerance */
         if (type == IOS_LO || type == IOS_DB || type == IOS_FX)
         {  temp1 = lb - ios->tol_int;
            temp2 = lb + ios->tol_int;
            if (temp1 <= x && x <= temp2) continue;
            /* if the solution is primal feasible, the lower bound must
               not be violated */
            if (ios->p_stat == IOS_FEAS) insist(x >= lb);
         }
         /* if the column's primal value is close to the upper bound,
            the column is integer feasible within given tolerance */
         if (type == IOS_UP || type == IOS_DB || type == IOS_FX)
         {  temp1 = ub - ios->tol_int;
            temp2 = ub + ios->tol_int;
            if (temp1 <= x && x <= temp2) continue;
            /* if the solution is primal feasible, the upper bound must
               not be violated */
            if (ios->p_stat == IOS_FEAS) insist(x <= ub);
         }
         /* if the column's primal value is close to nearest integer,
            the column is integer feasible within given tolerance */
         temp1 = floor(x + 0.5) - ios->tol_int;
         temp2 = floor(x + 0.5) + ios->tol_int;
         if (temp1 <= x && x <= temp2) continue;
         /* otherwise the column is integer infeasible */
         ios_get_col_ptr(ios, j)->frac = 1;
         /* increase the number of fractional-valued columns */
         ios->ii_cnt++;
         /* compute the sum of integer infeasibilities */
         temp1 = x - floor(x);
         temp2 = ceil(x) - x;
         insist(temp1 > 0.0 && temp2 > 0.0);
         ios->ii_sum += (temp1 <= temp2 ? temp1 : temp2);
      }
      if (ios->msg_lev >= 3)
      {  if (ios->ii_cnt == 0)
            print("There are no fractional columns");
         else if (ios->ii_cnt == 1)
            print("There is one fractional column, integer infeasibilit"
               "y is %.3e", ios->ii_sum);
         else
            print("There are %d fractional columns, integer infeasibili"
               "ty is %.3e", ios->ii_cnt, ios->ii_sum);
      }
      return;
}

/*----------------------------------------------------------------------
-- generate_rows - generate rows.
--
-- This routine calls the application procedure to generate additional
-- rows to be included in the current subproblem. */

static void generate_rows(IOS *ios)
{     int old_m, nrs;
      /* save the current number of rows */
      old_m = ios_get_num_rows(ios);
      /* call the application procedure to generate rows */
      ios->event = IOS_V_GENROW;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      /* determine how many rows have been generated */
      nrs = ios_get_num_rows(ios) - old_m;
      insist(nrs >= 0);
      if (ios->msg_lev >= 3)
      {  if (nrs == 1)
            print("One row has been generated");
         else if (nrs > 1)
            print("%d rows have been generated", nrs);
      }
      /* if at least one row has been generated, LP relaxation needs to
         be re-optimized */
      if (nrs > 0) ios->r_flag = 1;
      return;
}

/*----------------------------------------------------------------------
-- generate_cols - generate columns.
--
-- This routine calls the application procedure to generate additional
-- columns to be included in the current subproblem. */

static void generate_cols(IOS *ios)
{     int old_n, ncs;
      /* save the current number of columns */
      old_n = ios_get_num_cols(ios);
      /* call the application procedure to generate columns */
      ios->event = IOS_V_GENCOL;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      /* determine how many columns have been generated */
      ncs = ios_get_num_cols(ios) - old_n;
      insist(ncs >= 0);
      if (ios->msg_lev >= 3)
      {  if (ncs == 1)
            print("One column has been generated");
         else if (ncs > 1)
            print("%d columns have been generated", ncs);
      }
      /* if at least one column has been generated, LP relaxation needs
         to be re-optimized */
      if (ncs > 0) ios->r_flag = 1;
      return;
}

/*----------------------------------------------------------------------
-- set_local_bound - set local bound for current subproblem.
--
-- This routine sets the local bound for integer optimal solution of
-- the current subproblem, which is optimal solution of LP relaxation.
-- If the latter is fractional while the objective function is known to
-- be integral for any integer feasible point, the bound is strengthen
-- by rounding up (minimization) or down (maximization). */

static void set_local_bound(IOS *ios)
{     double bound;
      bound = ios->lp_obj;
      if (ios->int_obj)
      {  /* the objective function is known to be integral */
         double temp;
         temp = floor(bound + 0.5);
         if (temp - 1e-5 <= bound && bound <= temp + 1e-5)
            bound = temp;
         else
         {  switch (ios->dir)
            {  case IOS_MIN:
                  bound = ceil(bound); break;
               case IOS_MAX:
                  bound = floor(bound); break;
               default:
                  insist(ios != ios);
            }
         }
      }
      ios_get_npd_ptr(ios, ios_get_curr_node(ios))->bound = bound;
      if (ios->msg_lev >= 3)
         print("Local bound is %.9e", bound);
      return;
}

/*----------------------------------------------------------------------
-- is_branch_hopeful - check if specified branch is hopeful.
--
-- This routine checks if the specified subproblem can have an integer
-- optimal solution which is better than the best known one.
--
-- The check is based on comparison of the local objective bound stored
-- in the subproblem descriptor and the incumbent objective value which
-- is the global objective bound.
--
-- If there is a chance that the specified subproblem can have a better
-- integer optimal solution, the routine returns non-zero. Otherwise, if
-- the corresponding branch can pruned, zero is returned. */

static int is_branch_hopeful(IOS *ios, int p)
{     int ret = 1;
      if (ios->found)
      {  double bound, eps;
         bound = ios_get_npd_ptr(ios, p)->bound;
         eps = ios->tol_obj * (1.0 + fabs(ios->best));
         switch (ios->dir)
         {  case IOS_MIN:
               if (bound >= ios->best - eps) ret = 0;
               break;
            case IOS_MAX:
               if (bound <= ios->best + eps) ret = 0;
               break;
            default:
               insist(ios != ios);
         }
      }
      return ret;
}

/*----------------------------------------------------------------------
-- generate_cuts - generate cutting planes.
--
-- This routine calls the application procedure to generate additional
-- cutting planes to be included in the current subproblem. */

static void generate_cuts(IOS *ios)
{     int old_m, nrs;
      /* save the current number of rows */
      old_m = ios_get_num_rows(ios);
      /* call the application procedure to generate cuts */
      ios->event = IOS_V_GENCUT;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      /* determine how many cuts have been generated */
      nrs = ios_get_num_rows(ios) - old_m;
      insist(nrs >= 0);
      if (ios->msg_lev >= 3)
      {  if (nrs == 1)
            print("One cut has been generated");
         else if (nrs > 1)
            print("%d cuts have been generated", nrs);
      }
      /* if at least one cut has been generated, LP relaxation needs to
         be re-optimized */
      if (nrs > 0) ios->r_flag = 1;
      return;
}

/*----------------------------------------------------------------------
-- ios_branch_first - choose first column to branch on.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_branch_first(IOS *ios, int *next);
--
-- *Description*
--
-- The routine ios_branch_first chooses the first column in the column
-- list, which is of integer kind and whose primal value is fractional.
--
-- The routine also selects the branch to be solved next, where integer
-- infeasibility of the chosen column is less than in other branch. The
-- corresponding flag is stored to location next: -1 means down branch,
-- +1 means up branch.
--
-- *Returns*
--
-- The routine ios_branch_first returns the column number. */

int ios_branch_first(IOS *ios, int *next)
{     int n, j;
      double x;
      if (ios->event != IOS_V_BRANCH)
         fault("ios_branch_first: event != IOS_V_BRANCH; improper call "
            "sequence");
      /* choose column to branch on */
      n = ios_get_num_cols(ios);
      for (j = 1; j <= n; j++)
         if (ios_is_col_frac(ios, j)) break;
      insist(1 <= j && j <= n);
      /* obtain primal value of the column */
      ios_get_col_soln(ios, j, &x, NULL);
      /* select branch to be solved next */
      if (next != NULL)
         *next = x - floor(x) < ceil(x) - x ? -1 : +1;
      return j;
}

/*----------------------------------------------------------------------
-- ios_branch_last - choose last column to branch on.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_branch_last(IOS *ios, int *next);
--
-- *Description*
--
-- The routine ios_branch_last chooses the last column in the column
-- list, which is of integer kind and whose primal value is fractional.
--
-- The routine also selects the branch to be solved next, where integer
-- infeasibility of the chosen column is less than in other branch. The
-- corresponding flag is stored to location next: -1 means down branch,
-- +1 means up branch.
--
-- *Returns*
--
-- The routine ios_branch_last returns the column number. */

int ios_branch_last(IOS *ios, int *next)
{     int n, j;
      double x;
      if (ios->event != IOS_V_BRANCH)
         fault("ios_branch_last: event != IOS_V_BRANCH; improper call s"
            "equence");
      /* choose column to branch on */
      n = ios_get_num_cols(ios);
      for (j = n; j >= 1; j--)
         if (ios_is_col_frac(ios, j)) break;
      insist(1 <= j && j <= n);
      /* obtain primal value of the column */
      ios_get_col_soln(ios, j, &x, NULL);
      /* select branch to be solved next */
      if (next != NULL)
         *next = x - floor(x) < ceil(x) - x ? -1 : +1;
      return j;
}

/*----------------------------------------------------------------------
-- ios_branch_drtom - choose column using Driebeck-Tomlin heuristic.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_branch_drtom(IOS *ios, int *next);
--
-- *Description*
--
-- The routine ios_branch_drtom chooses a column to branch on using the
-- heuristic proposed by Driebeck and Tomlin.
--
-- The routine also selects the branch to be solved next and stores the
-- corresponding flag to location next: -1 means down branch, +1 means
-- up branch.
--
-- *References*
--
-- This routine is based on the heuristic proposed in:
--
-- Driebeck N.J. An algorithm for the solution of mixed-integer
-- programming problems, Management Science, 12: 576-87 (1966)
--
-- and improved in:
--
-- Tomlin J.A. Branch and bound methods for integer and non-convex
-- programming, in J.Abadie (ed.), Integer and Nonlinear Programming,
-- North-Holland, Amsterdam, pp. 437-50 (1970).
--
-- Note that this heuristic is time-expensive, because computing
-- one-step degradation (see the routine below) requires one BTRAN for
-- every fractional-valued column. */

#include "glplpx.h"

int ios_branch_drtom(IOS *ios, int *_next)
{     LPX *lp;
      int m, n, j, jj, k, t, next, kase, len, stat, *ind;
      double x, dk, alfa, delta_j, delta_k, delta_z, dz_dn, dz_up,
         dd_dn, dd_up, degrad, *val;
      if (ios->event != IOS_V_BRANCH)
         fault("ios_branch_drtom: event != IOS_V_BRANCH; improper call "
            "sequence");
      /* basic solution of LP relaxation is assumed to be optimal */
      insist(ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS);
      /* extract LP relaxation of the current subproblem */
      lp = ios_extract_lp(ios);
      /* factorize the basis to access the simplex table */
      if (ios->scale) lpx_scale_prob(lp);
      insist(lpx_warm_up(lp) == LPX_E_OK);
      if (lpx_get_status(lp) != LPX_OPT)
      {  /* this may only happen due to excessive round-off errors */
         if (ios->msg_lev >= 2)
            print("ios_branch_drtom: numeric instability");
         /* we must not change the basis to attain optimality, so use
            dirty engineering (don't worry, this is ok :+) */
#if 0
         lp->p_stat = LPX_P_FEAS;
         lp->d_stat = LPX_D_FEAS;
#else
         lpx_put_solution(lp, LPX_P_FEAS, LPX_D_FEAS, NULL, NULL, NULL,
            NULL, NULL, NULL);
         insist(lpx_is_b_avail(lp));
#endif
         insist(lpx_get_status(lp) == LPX_OPT);
      }
      /* determine number of rows and columns */
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      /* allocate working arrays */
      ind = ucalloc(1+n, sizeof(int));
      val = ucalloc(1+n, sizeof(double));
      /* nothing has been chosen so far */
      jj = 0;
      degrad = -1.0;
      /* walk through the list of columns (structural variables) */
      for (j = 1; j <= n; j++)
      {  /* if j-th column is not marked as fractional, skip it */
         if (!ios_is_col_frac(ios, j)) continue;
         /* obtain (fractional) value of j-th column in basic solution
            of LP relaxation */
         ios_get_col_soln(ios, j, &x, NULL);
         /* since the value of j-th column is fractional, the column is
            basic; compute corresponding row of the simplex table */
         len = lpx_eval_tab_row(lp, m+j, ind, val);
         /* the following fragment computes a change in the objective
            function: delta Z = new Z - old Z, where old Z is the
            objective value in the current optimal basis, and new Z is
            the objective value in the adjacent basis, for two cases:
            1) if new upper bound ub' = floor(x[j]) is introduced for
               j-th column (down branch);
            2) if new lower bound lb' = ceil(x[j]) is introduced for
               j-th column (up branch);
            since in both cases the solution remaining dual feasible
            becomes primal infeasible, one implicit simplex iteration
            is performed to determine the change delta Z;
            it is obvious that new Z, which is never better than old Z,
            is a lower (minimization) or upper (maximization) bound of
            the objective function for down- and up-branches. */
         for (kase = -1; kase <= +1; kase += 2)
         {  /* if kase < 0, the new upper bound of x[j] is introduced;
               in this case x[j] should decrease in order to leave the
               basis and go to its new upper bound */
            /* if kase > 0, the new lower bound of x[j] is introduced;
               in this case x[j] should increase in order to leave the
               basis and go to its new lower bound */
            /* apply the dual ratio test in order to determine which
               auxiliary or structural variable should enter the basis
               to keep dual feasibility */
            k = lpx_dual_ratio_test(lp, len, ind, val, kase, 1e-8);
            /* if no non-basic variable has been chosen, LP relaxation
               of corresponding branch being primal infeasible and dual
               unbounded has no primal feasible solution; in this case
               the change delta Z is formally set to infinity */
            if (k == 0)
            {  delta_z = (ios->dir == IOS_MIN ? +DBL_MAX : -DBL_MAX);
               goto skip;
            }
            /* row of the simplex table that corresponds to non-basic
               variable x[k] choosen by the dual ratio test is:
                  x[j] = ... + alfa * x[k] + ...
               where alfa is the influence coefficient (an element of
               the simplex table row) */
            /* determine the coefficient alfa */
            for (t = 1; t <= len; t++) if (ind[t] == k) break;
            insist(1 <= t && t <= len);
            alfa = val[t];
            /* since in the adjacent basis the variable x[j] becomes
               non-basic, knowing its value in the current basis we can
               determine its change delta x[j] = new x[j] - old x[j] */
            delta_j = (kase < 0 ? floor(x) : ceil(x)) - x;
            /* and knowing the coefficient alfa we can determine the
               corresponding change delta x[k] = new x[k] - old x[k],
               where old x[k] is a value of x[k] in the current basis,
               and new x[k] is a value of x[k] in the adjacent basis */
            delta_k = delta_j / alfa;
            /* Tomlin notices that if the variable x[k] is of integer
               kind, its change cannot be less (eventually) than one in
               the magnitude */
            if (k > m && ios_get_col_kind(ios, k-m) == IOS_INT)
            {  /* x[k] is structural integer variable */
               if (fabs(delta_k - floor(delta_k + 0.5)) > 1e-3)
               {  if (delta_k > 0.0)
                     delta_k = ceil(delta_k);  /* +3.14 -> +4 */
                  else
                     delta_k = floor(delta_k); /* -3.14 -> -4 */
               }
            }
            /* now determine the status and reduced cost of x[k] in the
               current basis */
            if (k <= m)
               stat = ios_get_row_soln(ios, k, NULL, &dk);
            else
               stat = ios_get_col_soln(ios, k-m, NULL, &dk);
            /* if the current basis is dual degenerative, some reduced
               costs which are close to zero may have wrong sign due to
               round-off errors, so correct the sign of d[k] */
            switch (ios->dir)
            {  case IOS_MIN:
                  if (stat == IOS_NL && dk < 0.0 ||
                      stat == IOS_NU && dk > 0.0 ||
                      stat == IOS_NF) dk = 0.0;
                  break;
               case IOS_MAX:
                  if (stat == IOS_NL && dk > 0.0 ||
                      stat == IOS_NU && dk < 0.0 ||
                      stat == IOS_NF) dk = 0.0;
                  break;
               default:
                  insist(ios != ios);
            }
            /* now knowing the change of x[k] and its reduced cost d[k]
               we can compute the corresponding change in the objective
               function delta Z = new Z - old Z = d[k] * delta x[k];
               note that due to Tomlin's modification new Z can be even
               worse than in the adjacent basis */
            delta_z = dk * delta_k;
skip:       /* new Z is never better than old Z, therefore the change
               delta Z is always non-negative (in case of minimization)
               or non-positive (in case of maximization) */
            switch (ios->dir)
            {  case IOS_MIN: insist(delta_z >= 0.0); break;
               case IOS_MAX: insist(delta_z <= 0.0); break;
               default: insist(ios != ios);
            }
            /* save the change in the objective fnction for down- and
               up-branches, respectively */
            if (kase < 0) dz_dn = delta_z; else dz_up = delta_z;
         }
         /* thus, in down-branch no integer feasible solution can be
            better than Z + dz_dn, and in up-branch no integer feasible
            solution can be better than Z + dz_up, where Z is value of
            the objective function in the current basis */
         /* following the heuristic by Driebeck and Tomlin we choose a
            column (i.e. structural variable) which provides largest
            degradation of the objective function in some of branches;
            besides, we select the branch with smaller degradation to
            be solved next and keep other branch with larger degradation
            in the active list hoping to minimize the number of further
            backtrackings */
         if (degrad < fabs(dz_dn) || degrad < fabs(dz_up))
         {  jj = j;
            if (fabs(dz_dn) < fabs(dz_up))
            {  /* select down branch to be solved next */
               next = -1;
               degrad = fabs(dz_up);
            }
            else
            {  /* select up branch to be solved next */
               next = +1;
               degrad = fabs(dz_dn);
            }
            /* save the objective changes for printing */
            dd_dn = dz_dn, dd_up = dz_up;
            /* if down- or up-branch has no feasible solution, we does
               not need to consider other candidates (in principle, the
               corresponding branch could be pruned right now) */
            if (degrad == DBL_MAX) break;
         }
      }
      /* free working arrays */
      ufree(ind);
      ufree(val);
      /* delete LP relaxation */
      lpx_delete_prob(lp);
      /* something must be chosen */
      insist(1 <= jj && jj <= n);
      if (ios->msg_lev >= 3)
      {  print("ios_branch_drtom: column %d chosen to branch on", jj);
         if (fabs(dd_dn) == DBL_MAX)
            print("ios_branch_drtom: down-branch is infeasible");
         else
            print("ios_branch_drtom: down-branch bound is %.9e",
               ios->lp_obj + dd_dn);
         if (fabs(dd_up) == DBL_MAX)
            print("ios_branch_drtom: up-branch   is infeasible");
         else
            print("ios_branch_drtom: up-branch   bound is %.9e",
               ios->lp_obj + dd_up);
      }
      if (_next != NULL) *_next = next;
      return jj;
}

#define LPX ???

/*----------------------------------------------------------------------
-- ios_branch_on - perform branching on specified column.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_branch_on(IOS *ios, int j, int next);
--
-- *Description*
--
-- The routine ios_branch_on performs branching on j-th column of the
-- current subproblem. The specified column (structural variable) must
-- be of integer kind and have a fractional value in basic solution of
-- LP relaxation of the current subproblem (i.e. only columns for which
-- ios_is_col_frac returns non-zero are valid candidates to branch on).
--
-- Let x be j-th structural variable, and beta be its primal fractional
-- value in the current basic solution. Branching on j-th variable is
-- dividing the current subproblem into two new subproblems, which are
-- identical to the current subproblem with the following exception: in
-- the first subproblem that begins the down-branch x has a new upper
-- bound x <= floor(beta), and in the second subproblem that begins the
-- up-branch x has a new lower bound x >= ceil(beta).
--
-- The parameter next specifies which subproblem should be solved next
-- to continue the search:
--
-- -1 means that one which begins the down-branch;
-- +1 means that one which begins the up-branch.
--
-- Note that the application program can call this routine only once on
-- processing the event IOS_V_BRANCH. */

void ios_branch_on(IOS *ios, int j, int next)
{     int p, type, clone[1+2];
      double beta, lb, ub, new_lb, new_ub;
      if (ios->event != IOS_V_BRANCH)
         fault("ios_branch_on: event != IOS_V_BRANCH; improper call seq"
            "uence");
      if (ios->b_flag)
         fault("ios_branch_on: branching already done");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_branch_on: j = %d; column number out of range", j);
      if (!ios_get_col_ptr(ios, j)->frac)
         fault("ios_branch_on: j = %d; non-fractional column not allowe"
            "d to branch on", j);
      if (!(next == -1 || next == +1))
         fault("ios_branch_on: next = %d; invalid parameter", next);
      /* obtain primal value of j-th column in basic solution */
      ios_get_col_soln(ios, j, &beta, NULL);
      if (ios->msg_lev >= 3)
         print("Branching on column %d, primal value is %.9e", j, beta);
      /* determine the reference number of the current subproblem */
      p = ios_get_curr_node(ios);
      /* freeze the current subproblem */
      ios_freeze_node(ios);
      /* create two clones of the current subproblem; the first clone
         begins the down-branch, the second one begins the up-branch */
      ios_clone_node(ios, p, 2, clone);
      if (ios->msg_lev >= 3)
         print("Node %d begins down branch, node %d begins up branch",
            clone[1], clone[2]);
      /* set new upper bound of j-th column in the first subproblem */
      ios_revive_node(ios, clone[1]);
      new_ub = floor(beta);
      type = ios_get_col_bnds(ios, j, &lb, &ub);
      switch (type)
      {  case IOS_FR:
            type = IOS_UP;
            break;
         case IOS_LO:
            insist(lb <= new_ub);
            type = (lb < new_ub ? IOS_DB : IOS_FX);
            break;
         case IOS_UP:
            insist(new_ub <= ub - 1.0);
            break;
         case IOS_DB:
            insist(lb <= new_ub && new_ub <= ub - 1.0);
            type = (lb < new_ub ? IOS_DB : IOS_FX);
            break;
         default:
            insist(type != type);
      }
      ios_set_col_bnds(ios, j, type, lb, new_ub);
      ios_freeze_node(ios);
      /* set new lower bound of j-th column in the second subproblem */
      ios_revive_node(ios, clone[2]);
      new_lb = ceil(beta);
      type = ios_get_col_bnds(ios, j, &lb, &ub);
      switch (type)
      {  case IOS_FR:
            type = IOS_LO;
            break;
         case IOS_LO:
            insist(lb + 1.0 <= new_lb);
            break;
         case IOS_UP:
            insist(new_lb <= ub);
            type = (new_lb < ub ? IOS_DB : IOS_FX);
            break;
         case IOS_DB:
            insist(lb + 1.0 <= new_lb && new_lb <= ub);
            type = (new_lb < ub ? IOS_DB : IOS_FX);
            break;
         default:
            insist(type != type);
      }
      ios_set_col_bnds(ios, j, type, new_lb, ub);
      ios_freeze_node(ios);
      /* set the branching flag */
      ios->b_flag = 1;
      /* revive the specified subproblem to be solved next */
      ios_revive_node(ios, clone[next < 0 ? 1 : 2]);
      return;
}

/*----------------------------------------------------------------------
-- cleanup_the_tree - prune hopeless branches of the tree.
--
-- This routine walks through the active list and checks the local
-- bound for every active subproblem. If the local bound indicates that
-- the subproblem cannot have integer optimal solution better than the
-- incumbent objective value, the routine deletes such subproblem that,
-- in turn, involves pruning the corresponding branch of the tree. */

static void cleanup_the_tree(IOS *ios)
{     int p, next, count;
      /* the global bound must exist */
      insist(ios->found);
      /* count is the number of active subproblems deleted */
      count = 0;
      /* walk through the list of active subproblems */
      for (p = ios_get_next_node(ios, 0); p != 0; p = next)
      {  /* save the reference number of the next active subproblem */
         next = ios_get_next_node(ios, p);
         /* if the branch is hopeless, prune it */
         if (!is_branch_hopeful(ios, p))
            ios_delete_node(ios, p), count++;
      }
      if (ios->msg_lev >= 3)
      {  if (count == 1)
            print("One hopeless branch has been pruned");
         else if (count > 1)
            print("%d hopeless branches have been pruned", count);
      }
      return;
}

/*----------------------------------------------------------------------
-- ios_select_fifo - select subproblem using FIFO heuristic.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_select_fifo(IOS *ios);
--
-- *Description*
--
-- The routine ios_select_fifo selects active subproblem to be solved
-- next using FIFO heuristic. This is the subproblem, which was added
-- to the active list most early.
--
-- If this heuristic is used, nodes of the search tree are visited in
-- the breadth-first-search manner.
--
-- NOTE: This heuristic is absolutely impractical. It is included here
-- for completeness only.
--
-- *Returns*
--
-- The routine ios_select_fifo returns the reference number (positive
-- integer) of the subproblem selected. */

int ios_select_fifo(IOS *ios)
{     if (ios->event != IOS_V_SELECT)
         fault("ios_select_fifo: event != IOS_V_SELECT; improper call s"
            "equence");
      return ios_get_next_node(ios, 0);
}

/*----------------------------------------------------------------------
-- ios_select_lifo - select subproblem using LIFO heuristic.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_select_lifo(IOS *ios);
--
-- *Description*
--
-- The routine ios_select_lifo selects active subproblem to be solved
-- next using LIFO heuristic. This is the subproblem, which was added
-- to the active list most recently.
--
-- If this heuristic is used, nodes of the search tree are visited in
-- the depth-first-search manner.
--
-- *Returns*
--
-- The routine ios_select_lifo returns the reference number (positive
-- integer) of the subproblem selected. */

int ios_select_lifo(IOS *ios)
{     if (ios->event != IOS_V_SELECT)
         fault("ios_select_lifo: event != IOS_V_SELECT; improper call s"
            "equence");
      return ios_get_prev_node(ios, 0);
}

/*----------------------------------------------------------------------
-- ios_select_node - select subproblem to continue the search.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_select_node(IOS *ios, int p);
--
-- *Description*
--
-- The routine ios_select_node selects the specified subproblem whose
-- reference number is p to continue the search.
--
-- Note that the application program can call this routine only once on
-- processing the event IOS_V_SELECT. */

void ios_select_node(IOS *ios, int p)
{     if (ios->event != IOS_V_SELECT)
         fault("ios_select_node: event != IOS_V_SELECT; improper call s"
            "equence");
      if (ios->t_flag)
         fault("ios_select_node: subproblem already selected");
      if (ios_get_node_lev(ios, p) < 0)
         fault("ios_select_node: p = %d; invalid subproblem reference n"
            "umber", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_select_node: p = %d; selecting inactive subproblem "
            "not allowed", p);
      /* set the backtracking flag */
      ios->t_flag = 1;
      /* revive the specified subproblem to be solved next */
      ios_revive_node(ios, p);
      return;
}

/*----------------------------------------------------------------------
-- ios_driver - integer optimization driver routine.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_driver(void (*appl)(IOS *ios, void *info), void *info);
--
-- *Description*
--
-- The routine ios_driver is a driver routine which manages the process
-- of solving the integer optimization problem.
--
-- The parameter appl specifies an entry point to the event-driven
-- application procedure supplied by the user. This procedure is called
-- by the driver at certain points of the optimization process in order
-- to perform some application specific actions. Its first parameter is
-- a pointer to a data structure (which is the main data structure of
-- the integer optimization suite) created by the driver. This pointer
-- should be passed to every suite routine called from the application
-- procedure. The second parameter is the transit informational pointer
-- passed to the routine ios_driver.
--
-- *Returns*
--
-- If the search has been successfully completed, the routine returns
-- zero. Otherwise the routine returns non-zero. */

int ios_driver(void (*appl)(IOS *ios, void *info), void *info)
{     IOS *ios;
      int p, ret;
      /* create the main data structure */
      ios = ios_create_tree(appl, info);
      /* revive the root subproblem (initially it is empty) */
      ios_revive_node(ios, 1);
      /* call the application procedure to perform initialization */
      ios->event = IOS_V_INIT;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      if (ios->msg_lev >= 3)
      {  int m, n, nz;
         m = ios_get_num_rows(ios);
         n = ios_get_num_cols(ios);
         nz = ios_get_num_nz(ios);
         print("Root subproblem has %d row%s, %d column%s, and %d non-z"
            "ero%s", m, m == 1 ? "" : "s", n, n == 1 ? "" : "s", nz,
            nz == 1 ? "" : "s");
      }
      /* check that the initial root subproblem is not empty, i.e. has
         at least one row and one column */
      if (!(ios_get_num_rows(ios) > 0 && ios_get_num_cols(ios) > 0))
      {  if (ios->msg_lev >= 1)
            print("ios_driver: root subproblem has no rows/columns");
         ret = 1;
         goto done;
      }
      /* now the optimization direction is known; initialize the lower
         (minimization) or upper (maximization) objective bound for the
         root subproblem */
      switch (ios->dir)
      {  case IOS_MIN:
            ios_get_npd_ptr(ios, 1)->bound = -DBL_MAX;
            break;
         case IOS_MAX:
            ios_get_npd_ptr(ios, 1)->bound = +DBL_MAX;
            break;
         default:
            insist(ios != ios);
      }
      /* solve LP relaxation of the initial root subproblem */
      if (ios->msg_lev >= 2)
         print("Solving initial LP relaxation...");
      if (ios_solve_root(ios) != 0)
      {  if (ios->msg_lev >= 1)
            print("ios_driver: cannot solve initial LP relaxation");
         ret = 1;
         goto done;
      }
      /* analyze status of the initial basic solution */
      if (ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS)
      {  /* LP relaxation has optimal solution */
         if (ios->msg_lev >= 2)
            print("FOUND OPTIMAL SOLUTION TO LP RELAXATION");
      }
      else if (ios->p_stat == IOS_FEAS && ios->d_stat == IOS_NOFEAS)
      {  /* LP relaxation has unbounded solution */
         if (ios->msg_lev >= 2)
         {  print("LP RELAXATION HAS UNBOUNDED SOLUTION");
            print("INTEGER OPTIMIZATION IMPOSSIBLE");
         }
         ret = 1;
         goto done;
      }
      else if (ios->p_stat == IOS_NOFEAS)
      {  /* LP relaxation has no primal feasible solution */
         if (ios->msg_lev >= 2)
            print("LP RELAXATION HAS NO FEASIBLE SOLUTION");
         /* if the column generation is disabled, there are no chances
            to attain primal feasibility */
         if (!ios->col_gen)
         {  if (ios->msg_lev >= 2)
               print("INTEGER OPTIMIZATION IMPOSSIBLE");
            ret = 1;
            goto done;
         }
      }
      else
      {  /* other cases are impossible */
         insist(ios != ios);
      }
      /* if the column generation is disabled, optimal solution of LP
         relaxation is the lower (minimization) or upper (maximization)
         objective bound for the root subproblem */
      if (!ios->col_gen)
      {  insist(ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS);
         ios_get_npd_ptr(ios, 1)->bound = ios->lp_obj;
      }
      /* initialization has been finished */
      if (ios->msg_lev >= 2)
      {  print("Integer optimization begins...");
         show_progress(ios);
      }
LOOP: /* major loop starts here; a branch to this label is performed
         whenever some subproblem has been just selected from the active
         list and made current; in the end of this loop the subproblem
         will be either fathomed (and then its branch will be pruned) or
         divided into other subproblems */
      /* determine the reference number of the current subproblem */
      p = ios_get_curr_node(ios);
      insist(p != 0);
      if (ios->msg_lev >= 3)
      {  int level = ios_get_node_lev(ios, p);
         print("-------------------------------------------------------"
            "-----------------");
         print("Processing node %d at level %d", p, level);
      }
loop: /* minor loop starts here; a branch to this label is performed
         either in the beginning of the major loop when LP relaxation of
         the current subproblem needs to be solved for the first time or
         whenever the current subproblem has been changed and therefore
         its LP relaxation needs to be re-optimized */
      /* invalidate basic solution of the current LP relaxation */
      clear_lp_soln(ios);
      /* clear the re-optimization flag */
      ios->r_flag = 0;
      /* clear the branching flag */
      ios->b_flag = 0;
      /* clear the backtracking flag */
      ios->t_flag = 0;
      /* display the current progress of the search */
      if (ios->msg_lev >= 3 || ios->msg_lev >= 2 &&
          utime() - ios->tm_lag >= ios->out_frq - 0.001)
         show_progress(ios);
      /* solve LP relaxation of the current subproblem */
      if (ios->msg_lev >= 3)
         print("Solving LP relaxation: %d row(s), %d column(s), %d non-"
            "zero(s)", ios_get_num_rows(ios), ios_get_num_cols(ios),
            ios_get_num_nz(ios));
      if (ios_solve_node(ios) != 0)
      {  if (ios->msg_lev >= 1)
         {  if (ios->msg_lev >= 2) show_progress(ios);
            print("ios_driver: cannot solve current LP relaxation");
         }
         ret = 1;
         goto done;
      }
      /* check integrality of the basic solution */
      check_integrality(ios);
      /* analyze status of the basic solution */
      if (ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS)
      {  /* LP relaxation has optimal solution */
         if (ios->msg_lev >= 3)
            print("Found optimal solution to LP relaxation");
         /* if the column generation is enabled, make sure that there
            are no variables which can improve the optimal solution and
            are still not included in the current subproblem */
         if (ios->col_gen)
         {  generate_cols(ios);
            if (ios->r_flag) goto loop;
         }
      }
      else if (ios->p_stat == IOS_FEAS && ios->d_stat == IOS_NOFEAS)
      {  /* LP relaxation has unbounded solution */
         /* since the current subproblem cannot have a larger feasible
            region than its parent, there is something wrong */
         if (ios->msg_lev >= 1)
            print("ios_driver: current LP relaxation has unbounded solu"
               "tion");
         ret = 1;
         goto done;
      }
      else if (ios->p_stat == IOS_INFEAS && ios->d_stat == IOS_FEAS)
      {  /* LP relaxation has no primal solution which is better than
            the incumbent objective value */
         /* this is possible only if the column generation is disabled
            and the dual simplex prematurely terminated the search */
         insist(!ios->col_gen);
         insist(ios->found);
         if (ios->msg_lev >= 3)
            print("LP relaxation has no solution better than incumbent "
               "objective value");
         /* prune the branch */
         goto fath;
      }
      else if (ios->p_stat == IOS_NOFEAS)
      {  /* LP relaxation has no primal feasible solution */
         if (ios->msg_lev >= 3)
            print("LP relaxation has no feasible solution");
         /* if the column generation is enabled, make sure that there
            are no variables which can reduce primal infeasibility and
            are still not included in the current subproblem */
         if (ios->col_gen)
         {  generate_cols(ios);
            if (ios->r_flag) goto loop;
         }
         /* prune the branch */
         goto fath;
      }
      else
      {  /* other cases cannot appear */
         insist(ios != ios);
      }
      /* at this point basic solution of LP relaxation of the current
         subproblem is optimal for the complete set of columns */
      insist(ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS);
      /* thus, it defines a local bound for integer optimal solution of
         the current subproblem */
      set_local_bound(ios);
      /* if the local bound indicates that integer optimal solution of
         the current subproblem cannot be better than the global bound,
         prune the branch */
      if (!is_branch_hopeful(ios, p))
      {  if (ios->msg_lev >= 3)
            print("Current branch is hopeless and can be pruned");
         goto fath;
      }
      /* if the row generation is enabled, make sure that there are no
         constraints which are violated but not included in the current
         subproblem */
      if (ios->row_gen)
      {  generate_rows(ios);
         if (ios->r_flag) goto loop;
      }
      /* at this point basic solution of LP relaxation of the current
         subproblem is optimal for the complete set of rows and columns
         and is better than the global bound */
      insist(ios->p_stat == IOS_FEAS && ios->d_stat == IOS_FEAS);
      /* if the basic solution satifies to all integrality conditions,
         it is a new, better integer feasible solution */
      if (ios->ii_cnt == 0)
      {  if (ios->msg_lev >= 3)
            print("New integer feasible solution found");
         ios->found = 1;
         ios->best = ios->lp_obj;
         if (ios->msg_lev >= 2) show_progress(ios);
         /* call the application procedure to make it happy */
         ios->event = IOS_V_BINGO;
         ios->appl(ios, ios->info);
         ios->event = IOS_V_NONE;
         /* the current subproblem is fathomed; prune its branch */
         goto fath;
      }
      /* basic solution of LP relaxation of the current subproblem is
         integer infeasible */
      /* if the cutting plane generation is enabled, try to add one or
         more cutting planes in order to cut off the current fractional
         point from the integer polytope */
      if (ios->cut_gen)
      {  generate_cuts(ios);
         if (ios->r_flag) goto loop;
      }
      /* the current subproblem is near to become inactive; store some
         additional information to its descriptor */
      ios_get_npd_ptr(ios, p)->ii_cnt = ios->ii_cnt;
      ios_get_npd_ptr(ios, p)->ii_sum = ios->ii_sum;
      /* call the application procedure to perform branching */
      ios->event = IOS_V_BRANCH;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      if (!ios->b_flag)
         fault("ios_driver: branching required but not performed");
      /* continue the search from the subproblem specified */
      goto LOOP;
fath: /* the current subproblem has been fathomed */
      if (ios->msg_lev >= 3)
         print("Node %d fathomed", p);
      /* freeze the current subproblem */
      ios_freeze_node(ios);
      /* and prune the corresponding branch of the tree */
      ios_delete_node(ios, p);
      /* if a new integer feasible solution has just been found, other
         branches may become hopeless and therefore should be pruned */
      if (ios->found) cleanup_the_tree(ios);
      /* if the active list is empty, the search is finished */
      if (ios_get_next_node(ios, 0) == 0) goto fini;
      /* call the application procedure to perform backtracking */
      ios->event = IOS_V_SELECT;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      if (!ios->t_flag)
         fault("ios_driver: backtracking required but not performed");
      /* contiunue the search from the subproblem selected */
      goto LOOP;
fini: /* the search has been successfully completed */
      if (ios->msg_lev >= 3)
         print("Active list is empty!");
      if (ios->msg_lev >= 2)
      {  show_progress(ios);
         if (ios->found)
            print("INTEGER OPTIMAL SOLUTION FOUND");
         else
            print("PROBLEM HAS NO INTEGER FEASIBLE SOLUTION");
      }
      ret = 0;
done: /* if the current subproblem still exists, freeze it */
      if (iet_get_curr_node(ios->iet) != 0) ios_freeze_node(ios);
      /* delete all active subproblems which are still in the tree */
      for (;;)
      {  p = ios_get_next_node(ios, 0);
         if (p == 0) break;
         ios_delete_node(ios, p);
      }
      /* call the application procedure to perform termination */
      ios->event = IOS_V_TERM;
      ios->appl(ios, ios->info);
      ios->event = IOS_V_NONE;
      /* delete all the data structures to free the memory */
      ios_delete_tree(ios);
      /* return to the calling program */
      return ret;
}

/* eof */

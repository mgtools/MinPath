/* glpios1.c */

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
#include <stddef.h>
#include "glpios.h"
#include "glplib.h"

/**********************************************************************/
/* * *               LOW-LEVEL MAINTENANCE ROUTINES               * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_attach_npd - attach extension to subproblem descriptor.
--
-- This routine creates and attaches extension to the descriptor of the
-- subproblem whose reference number is p.
--
-- Detaching is performed automatically in the hook routine on deleting
-- the corresponding subproblem. */

void ios_attach_npd(IOS *ios, int p)
{     IOSNPD *node;
      node = ios_get_npd_ptr(ios, p);
      insist(node == NULL);
      node = dmp_get_atom(ios->npd_pool);
      node->bound = 0.0;
      node->ii_cnt = 0;
      node->ii_sum = 0.0;
      ios_set_npd_ptr(ios, p, node);
      return;
}

/*----------------------------------------------------------------------
-- ios_attach_rgd - attach extension to row global descriptor.
--
-- This routine creates and attaches extension to the global descriptor
-- of i-th row of the current subproblem.
--
-- Detaching is performed automatically in the hook routine on deleting
-- the corresponding row. */

void ios_attach_rgd(IOS *ios, int i)
{     IOSRGD *rgd;
      rgd = ios_get_rgd_ptr(ios, i);
      insist(rgd == NULL);
      rgd = dmp_get_atom(ios->rgd_pool);
      rgd->mark = 0;
      rgd->link = NULL;
      ios_set_rgd_ptr(ios, i, rgd);
      return;
}

/*----------------------------------------------------------------------
-- ios_attach_cgd - attach extension to column global descriptor.
--
-- This routine creates and attaches extension to the global descriptor
-- of j-th column of the current subproblem.
--
-- Detaching is performed automatically in the hook routine on deleting
-- the corresponding column. */

void ios_attach_cgd(IOS *ios, int j)
{     IOSCGD *cgd;
      cgd = ios_get_cgd_ptr(ios, j);
      insist(cgd == NULL);
      cgd = dmp_get_atom(ios->cgd_pool);
      cgd->kind = IOS_NUM;
      cgd->mark = 0;
      cgd->link = NULL;
      ios_set_cgd_ptr(ios, j, cgd);
      return;
}

/*----------------------------------------------------------------------
-- ios_attach_row - attach extension to row local descriptor.
--
-- This routine creates and attaches extension to the local descriptor
-- of i-th row of the current subproblem. */

void ios_attach_row(IOS *ios, int i)
{     IOSROW *row;
      row = ios_get_row_ptr(ios, i);
      insist(row == NULL);
      row = dmp_get_atom(ios->row_pool);
      row->prim = 0.0;
      row->dual = 0.0;
      row->pi = 0.0;
      ios_set_row_ptr(ios, i, row);
      return;
}

/*----------------------------------------------------------------------
-- ios_attach_col - attach extension to column local descriptor.
--
-- This routine creates and attaches extension to the local descriptor
-- of j-th column of the current subproblem. */

void ios_attach_col(IOS *ios, int j)
{     IOSCOL *col;
      col = ios_get_col_ptr(ios, j);
      insist(col == NULL);
      col = dmp_get_atom(ios->col_pool);
      col->prim = 0.0;
      col->dual = 0.0;
      col->frac = 0;
      ios_set_col_ptr(ios, j, col);
      return;
}

/*----------------------------------------------------------------------
-- ios_detach_row - detach extension from row local descriptor.
--
-- This routine detaches and deletes extension of the local descriptor
-- of i-th row of the current subproblem. */

void ios_detach_row(IOS *ios, int i)
{     IOSROW *row;
      row = ios_get_row_ptr(ios, i);
      insist(row != NULL);
      dmp_free_atom(ios->row_pool, row);
      row = NULL;
      ios_set_row_ptr(ios, i, row);
      return;
}

/*----------------------------------------------------------------------
-- ios_detach_col - detach extension from column local descriptor.
--
-- This routine detaches and deletes extension of the local descriptor
-- of j-th column of the current subproblem. */

void ios_detach_col(IOS *ios, int j)
{     IOSCOL *col;
      col = ios_get_col_ptr(ios, j);
      insist(col != NULL);
      dmp_free_atom(ios->col_pool, col);
      col = NULL;
      ios_set_col_ptr(ios, j, col);
      return;
}

/*----------------------------------------------------------------------
-- ios_hook_routine - callback interface to enumeration tree.
--
-- This hook routine is a callback interface to the enumeration tree.
-- It is called whenever a subproblem, row, or column is being deleted.
-- Being called it, in turn, calls the application procedure in order
-- to notify the latter about deletion and then detaches and deletes
-- extension of the global descriptor of the corresponding object. */

void ios_hook_routine(void *info, int what, char *name, void *link)
{     IOS *ios = info;
      switch (what)
      {  case IET_ND:
            /* some subproblem is being deleted */
            ios->hook_name = name;
            ios->hook_link.npd = link;
            /* notify the application procedure */
            ios->event = IOS_V_DELSUB;
            ios->appl(ios, ios->info);
            ios->event = IOS_V_NONE;
            /* delete extension of the subproblem descriptor */
            dmp_free_atom(ios->npd_pool, ios->hook_link.npd);
            ios->hook_name = NULL;
            ios->hook_link.npd = NULL;
            break;
         case IET_RD:
            /* some row is being deleted */
            ios->hook_name = name;
            ios->hook_link.rgd = link;
            /* notify the application procedure */
            ios->event = IOS_V_DELROW;
            ios->appl(ios, ios->info);
            ios->event = IOS_V_NONE;
            /* delete extension of the row global descriptor */
            dmp_free_atom(ios->rgd_pool, ios->hook_link.rgd);
            ios->hook_name = NULL;
            ios->hook_link.rgd = NULL;
            break;
         case IET_CD:
            /* some column is being deleted */
            ios->hook_name = name;
            ios->hook_link.cgd = link;
            /* notify the application procedure */
            ios->event = IOS_V_DELCOL;
            ios->appl(ios, ios->info);
            ios->event = IOS_V_NONE;
            /* delete extension of the column global descriptor */
            dmp_free_atom(ios->cgd_pool, ios->hook_link.cgd);
            ios->hook_name = NULL;
            ios->hook_link.cgd = NULL;
            break;
         default:
            insist(what != what);
      }
      return;
}

/**********************************************************************/
/* * *                  TREE MANAGEMENT ROUTINES                  * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_create_tree - create integer optimization suite.
--
-- This routines creates and initializes the main data structure of the
-- integer optimization suite.
--
-- At the beginning the implicit enumeration tree consists of the only
-- root subproblem which is empty, i.e. has no rows and columns, and is
-- assigned the reference number 1. Note that the root subproblem is in
-- the frozen state and therefore needs to be revived. */

IOS *ios_create_tree(void (*appl)(IOS *ios, void *info), void *info)
{     IOS *ios;
      /* create the main data structure */
      ios = umalloc(sizeof(IOS));
      ios->npd_pool = dmp_create_pool(sizeof(IOSNPD));
      ios->rgd_pool = dmp_create_pool(sizeof(IOSRGD));
      ios->cgd_pool = dmp_create_pool(sizeof(IOSCGD));
      ios->row_pool = dmp_create_pool(sizeof(IOSROW));
      ios->col_pool = dmp_create_pool(sizeof(IOSCOL));
      ios->iet = NULL;
      ios->hook_name = NULL;
      ios->hook_link.npd = NULL;
      ios->hook_link.rgd = NULL;
      ios->hook_link.cgd = NULL;
      ios->dir = IOS_MIN;
      ios->int_obj = 0;
      ios->row_gen = 0;
      ios->col_gen = 0;
      ios->cut_gen = 0;
      ios->found = 0;
      ios->best = 0.0;
      ios->p_stat = IOS_UNDEF;
      ios->d_stat = IOS_UNDEF;
      ios->lp_obj = 0.0;
      ios->lp_sum = 0.0;
      ios->ii_cnt = 0;
      ios->ii_sum = 0.0;
      ios->msg_lev = 2;
      ios->init_lp = 1;
      ios->scale = 0;
      ios->tol_int = 1e-5;
      ios->tol_obj = 1e-7;
      ios->out_frq = 5.0;
      ios->out_dly = 10.0;
      ios->it_cnt = 0;
      ios->tm_beg = utime();
      ios->tm_lag = 0.0;
      ios->appl = appl;
      ios->info = info;
      ios->event = IOS_V_NONE;
      ios->r_flag = 0;
      ios->b_flag = 0;
      ios->t_flag = 0;
      /* create the implicit enumeration tree */
      ios->iet = iet_create_tree();
      iet_install_hook(ios->iet, ios_hook_routine, ios);
      /* attach extension to the descriptor of the root subproblem */
      ios_attach_npd(ios, 1);
      return ios;
}

/*----------------------------------------------------------------------
-- ios_revive_node - revive specified subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_revive_node(IOS *ios, int p);
--
-- *Description*
--
-- The routine ios_revive_node revives the specified subproblem, whose
-- reference number is p, and thereby makes it the current subproblem.
-- Note that the specified subproblem must be active. Besides, if the
-- current subproblem already exists, it must be frozen before reviving
-- another subproblem. */

void ios_revive_node(IOS *ios, int p)
{     int m, n, i, j;
      if (ios_get_node_lev(ios, p) < 0)
         fault("ios_revive_node: p = %d; invalid subproblem reference n"
            "umber", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_revive_node: p = %d; reviving inactive subproblem n"
            "ot allowed", p);
      if (ios_get_curr_node(ios) != 0)
         fault("ios_revive_node: current subproblem already exists");
      /* revive the specified subproblem */
      iet_revive_node(ios->iet, p);
      /* attach extensions to local descriptors of rows and columns of
         the subproblem just revived */
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      for (i = 1; i <= m; i++) ios_attach_row(ios, i);
      for (j = 1; j <= n; j++) ios_attach_col(ios, j);
      /* basic solution components are not computed yet */
      ios->p_stat = IOS_UNDEF;
      ios->d_stat = IOS_UNDEF;
      ios->lp_obj = 0.0;
      ios->ii_cnt = 0;
      ios->ii_sum = 0.0;
      return;
}

/*----------------------------------------------------------------------
-- ios_freeze_node - freeze current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_freeze_node(IOS *ios);
--
-- *Description*
--
-- The routine ios_freeze_node freezes the current subproblem. */

void ios_freeze_node(IOS *ios)
{     int m, n, i, j;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_freeze_node: current subproblem does not exist");
      /* detach extensions from local descriptors of rows and columns
         of the current subproblem */
      m = ios_get_num_rows(ios);
      n = ios_get_num_cols(ios);
      for (i = 1; i <= m; i++) ios_detach_row(ios, i);
      for (j = 1; j <= n; j++) ios_detach_col(ios, j);
      /* freeze the current subproblem */
      iet_freeze_node(ios->iet);
      return;
}

/*----------------------------------------------------------------------
-- ios_clone_node - clone specified subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_clone_node(IOS *ios, int p, int nnn, int ref[]);
--
-- *Description*
--
-- The routine ios_clone_node clones the specified subproblem, whose
-- reference number is p, creating its nnn exact copies. Note that the
-- specified subproblem must be active and must be in the frozen state
-- (i.e. it must not be the current subproblem).
--
-- Each clone, an exact copy of the specified subproblem, becomes a new
-- active subproblem added to the end of the active list. After cloning
-- the specified subproblem becomes inactive.
--
-- The reference numbers of clone subproblems are stored to locations
-- ref[1], ..., ref[nnn]. */

void ios_clone_node(IOS *ios, int p, int nnn, int ref[])
{     IOSNPD *orig, *node;
      if (ios_get_node_lev(ios, p) < 0)
         fault("ios_clone_node: p = %d; invalid subproblem reference nu"
            "mber", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_clone_node: p = %d; cloning inactive subproblem not"
            " allowed", p);
      if (ios_get_curr_node(ios) == p)
         fault("ios_clone_node: p = %d; cloning current subproblem not "
            "allowed", p);
      if (nnn < 1)
         fault("ios_clone_node: nnn = %d; invalid number of clone subpr"
            "oblems", nnn);
      /* obtain pointer to extension of the descriptor of the specified
         subproblem */
      orig = ios_get_npd_ptr(ios, p);
      /* create nnn clones of the specified subproblem */
      iet_clone_node(ios->iet, p, nnn);
      /* now the clone subproblems just created are added in the end of
         the active list */
      p = ios_get_prev_node(ios, 0);
      for (p = p; nnn > 0; p = ios_get_prev_node(ios, p))
      {  /* store the clone reference number */
         ref[nnn--] = p;
         /* attach extension of the descriptor to the clone subproblem
            and inherit some information from the parent subproblem */
         ios_attach_npd(ios, p);
         node = ios_get_npd_ptr(ios, p);
         node->bound = orig->bound;
      }
      insist(nnn == 0);
      return;
}

/*----------------------------------------------------------------------
-- ios_delete_node - delete specified subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_delete_node(IOS *ios, int p);
--
-- *Description*
--
-- The routine ios_delete_node deletes the specified subproblem, whose
-- reference number is p. The subproblem must be active and must be in
-- the frozen state (i.e. it must not be the current subproblem).
--
-- Note that deletion is performed recursively, i.e. if a subproblem to
-- be deleted is the only child of its parent, the parent subproblem is
-- also deleted, etc. */

void ios_delete_node(IOS *ios, int p)
{     if (ios_get_node_lev(ios, p) < 0)
         fault("ios_delete_node: p = %d; invalid subproblem reference n"
            "umber", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_delete_node: p = %d; deleting inactive subproblem n"
            "ot allowed", p);
      if (ios_get_curr_node(ios) == p)
         fault("ios_delete_node: p = %d; deleting current subproblem no"
            "t allowed", p);
      iet_delete_node(ios->iet, p);
      return;
}

/*----------------------------------------------------------------------
-- ios_delete_tree - delete integer optimization suite.
--
-- This routine deletes all the data structures and thereby frees all
-- the memory allocated to the integer optimization suite. Note that if
-- the current subproblem exists, it must be frozen before call to this
-- routine. */

void ios_delete_tree(IOS *ios)
{     if (ios_get_curr_node(ios) != 0)
         fault("ios_delete_tree: current subproblem still exists");
      /* delete implicit enumeration tree */
      iet_delete_tree(ios->iet);
      /* now all atoms must return to their memory pools */
      insist(ios->npd_pool->count == 0);
      insist(ios->rgd_pool->count == 0);
      insist(ios->cgd_pool->count == 0);
      insist(ios->row_pool->count == 0);
      insist(ios->col_pool->count == 0);
      /* free all the memory */
      dmp_delete_pool(ios->npd_pool);
      dmp_delete_pool(ios->rgd_pool);
      dmp_delete_pool(ios->cgd_pool);
      dmp_delete_pool(ios->row_pool);
      dmp_delete_pool(ios->col_pool);
      ufree(ios);
      return;
}

/**********************************************************************/
/* * *                  TREE EXPLORING ROUTINES                   * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_get_curr_node - determine current active subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_curr_node(IOS *ios);
--
-- *Returns*
--
-- The routine ios_get_curr_node returns the reference number of the
-- current active subproblem. However, if the current subproblem does
-- not exist, the routine returns zero. */

int ios_get_curr_node(IOS *ios)
{     int p;
      p = iet_get_curr_node(ios->iet);
      return p;
}

/*----------------------------------------------------------------------
-- ios_get_next_node - determine next active subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_next_node(IOS *ios, int p);
--
-- *Returns*
--
-- If the parameter p is zero, the routine ios_get_next_node returns
-- the reference number of the first active subproblem. However, if the
-- tree is empty, zero is returned.
--
-- If the parameter p is not zero, it must specify the reference number
-- of some active subproblem, in which case the routine returns the
-- reference number of the next active subproblem. However, if there is
-- no next active subproblem in the list, zero is returned.
--
-- All subproblems in the active list are ordered chronologically, i.e.
-- subproblem A precedes subproblem B if A was created before B. */

int ios_get_next_node(IOS *ios, int p)
{     if (p == 0) goto next;
      if (ios_get_node_lev(ios, p) < 0)
         fault("ios_get_next_node: p = %d; invalid subproblem reference"
            " number", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_get_next_node: p = %d; subproblem not in the active"
            " list", p);
next: p = iet_get_next_node(ios->iet, p);
      return p;
}

/*----------------------------------------------------------------------
-- ios_get_prev_node - determine previous active subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_prev_node(IOS *ios, int p);
--
-- *Returns*
--
-- If the parameter p is zero, the routine ios_get_prev_node returns
-- the reference number of the last active subproblem. However, if the
-- tree is empty, zero is returned.
--
-- If the parameter p is not zero, it must specify the reference number
-- of some active subproblem, in which case the routine returns the
-- reference number of the previous active subproblem. However, if there
-- is no previous active subproblem in the list, zero is returned.
--
-- All subproblems in the active list are ordered chronologically, i.e.
-- subproblem A precedes subproblem B if A was created before B. */

int ios_get_prev_node(IOS *ios, int p)
{     if (p == 0) goto prev;
      if (ios_get_node_lev(ios, p) < 0)
         fault("ios_get_prev_node: p = %d; invalid subproblem reference"
            " number", p);
      if (ios_get_node_cnt(ios, p) > 0)
         fault("ios_get_prev_node: p = %d; subproblem not in the active"
            " list", p);
prev: p = iet_get_prev_node(ios->iet, p);
      return p;
}

/*----------------------------------------------------------------------
-- ios_get_up_node - determine parent subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_up_node(IOS *ios, int p);
--
-- *Returns*
--
-- The parameter p must specify the reference number of some (active or
-- inactive) subproblem, in which case the routine ios_get_up_node
-- returns the reference number of its parent subproblem. However, if
-- the specified subproblem is the root of the tree and therefore has no
-- parent, the routine returns zero. */

int ios_get_up_node(IOS *ios, int p)
{     if (ios_get_node_lev(ios, p) < 0)
         fault("ios_get_up_node: p = %d; invalid subproblem reference n"
            "umber", p);
      p = iet_get_up_node(ios->iet, p);
      return p;
}

/*----------------------------------------------------------------------
-- ios_get_node_lev - determine subproblem level.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_node_lev(IOS *ios, int p);
--
-- *Returns*
--
-- The routine ios_get_node_lev returns the level which the subproblem
-- whose reference number is p has in the tree. (The root subproblem
-- has the level 0, and the level of any other subproblem is the level
-- of its parent plus one.) However, if the parameter p is not a valid
-- subproblem reference number, the routine returns negative value. */

int ios_get_node_lev(IOS *ios, int p)
{     int level;
      level = iet_get_node_lev(ios->iet, p);
      return level;
}

/*----------------------------------------------------------------------
-- ios_get_node_cnt - determine number of child subproblems.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_node_cnt(IOS *ios, int p);
--
-- *Returns*
--
-- The routine ios_get_node_cnt returns the number of child subproblems
-- for the subproblem whose reference number is p.
--
-- Zero means the subproblem p is active and therefore has no childs.
--
-- Positive value means the subproblem p is inactive where the value is
-- the number of its childs.
--
-- Negative value means the reference number p is invalid. */

int ios_get_node_cnt(IOS *ios, int p)
{     int count;
      count = iet_get_node_cnt(ios->iet, p);
      return count;
}

/*----------------------------------------------------------------------
-- ios_pseudo_root - find pseudo-root of the tree.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_pseudo_root(IOS *ios);
--
-- *Description*
--
-- The routine ios_pseudo_root finds so-called pseudo-root of the tree,
-- a node which is the highest-leveled common ancestor of all active
-- nodes including the current one . (Juenger calls such node "the root
-- of the remaining tree".)
--
-- If walk from the root of the tree, the pseudo-root is the first node
-- which has more than one child:
--
--                 root -->    A              Level 0
--                             |
--                             B              Level 1
--                             |
--          pseudo-root -->    C              Level 2
--                            / \
--                          /     \
--                        /         \
--                       D           G        Level 3
--                     /   \       /   \
--                    E     F     H     I     Level 4
--                   ...   ...   ...   ...
--
-- (However, if the tree has the only active node, the pseudo-root is
-- that active node by definition.)
--
-- The subproblem associated with the pseudo-root has the property that
-- the local lower (minimization) or upper (maximization) bound of its
-- integer optimal solution is the global bound for the entire problem.
--
-- *Returns*
--
-- The routine ios_pseudo_root returns the subproblem reference number
-- which corresponds to the pseudo-root. However, if the tree is empty,
-- zero is returned. */

int ios_pseudo_root(IOS *ios)
{     int p;
      p = iet_pseudo_root(ios->iet);
      return p;
}

/**********************************************************************/
/* * *               SUBPROBLEM MODIFYING ROUTINES                * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_set_obj_dir - set optimization direction flag.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_obj_dir(IOS *ios, int dir);
--
-- *Description*
--
-- The routine ios_set_obj_dir sets the optimization direction flag as
-- specified by the parameter dir:
--
-- IOS_MIN means minimization;
-- IOS_MAX means maximization.
--
-- NOTE: Changing the optimization direction has the global effect. */

void ios_set_obj_dir(IOS *ios, int dir)
{     if (!(dir == IOS_MIN || dir == IOS_MAX))
         fault("ios_set_obj_dir: dir = %d; invalid parameter");
      ios->dir = dir;
      return;
}

/*----------------------------------------------------------------------
-- ios_add_rows - add new rows to current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_add_rows(IOS *ios, int nrs);
--
-- *Description*
--
-- The routine ios_add_rows adds nrs rows (constraints) to the current
-- subproblem. New rows are always added to the end of the row list, so
-- the ordinal numbers assigned to existing rows are not changed.
--
-- Each new row is initially free (unbounded) and has empty list of the
-- constraint coefficients. */

void ios_add_rows(IOS *ios, int nrs)
{     int m, i;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_add_rows: current subproblem does not exist");
      if (nrs < 1)
         fault("ios_add_rows: nrs = %d; invalid number of rows", nrs);
      /* determine the number of rows before addition */
      m = iet_get_num_rows(ios->iet);
      /* add nrs rows to the current subproblem */
      iet_add_rows(ios->iet, nrs);
      /* attach extensions to descriptors of new rows */
      for (i = m+1; i <= m+nrs; i++)
      {  ios_attach_rgd(ios, i);
         ios_attach_row(ios, i);
      }
      return;
}

/*----------------------------------------------------------------------
-- ios_add_cols - add new columns to current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_add_cols(IOS *ios, int ncs);
--
-- *Description*
--
-- The routine ios_add_cols adds ncs columns (structural variables) to
-- the current subproblem. New columns are always added to the end of
-- the column list, so the ordinal numbers assigned to existng columns
-- are not changed.
--
-- Each new column is initially fixed at zero and has empty list of the
-- constraint coefficients. */

void ios_add_cols(IOS *ios, int ncs)
{     int n, j;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_add_cols: current subproblem does not exist");
      if (ncs < 1)
         fault("ios_add_cols: ncs = %d; invalid number of columns",
            ncs);
      /* determine the number of columns before addition */
      n = iet_get_num_cols(ios->iet);
      /* add ncs columns to the current subproblem */
      iet_add_cols(ios->iet, ncs);
      /* attach extensions to descriptors of new columns */
      for (j = n+1; j <= n+ncs; j++)
      {  ios_attach_cgd(ios, j);
         ios_attach_col(ios, j);
      }
      return;
}

/*----------------------------------------------------------------------
-- ios_check_name - check correctness of symbolic name.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_check_name(IOS *ios, char *name);
--
-- *Description*
--
-- The routine ios_check_name checks if given symbolic name is correct
-- (a name is correct if it contains 1 up to 255 graphic characters).
--
-- *Returns*
--
-- If the symbolic name is correct, the routine returns zero. Otherwise
-- non-zero is returned. */

int ios_check_name(IOS *ios, char *name)
{     int ret;
      ret = iet_check_name(ios->iet, name);
      return ret;
}

/*----------------------------------------------------------------------
-- ios_set_row_name - assign symbolic name to row.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_row_name(IOS *ios, int i, char *name);
--
-- *Description*
--
-- The routine ios_set_row_name assigns a given symbolic name to i-th
-- row of the current subproblem.
--
-- If the parameter name is NULL, the routine just erases existing name
-- of the i-th row.
--
-- NOTE: Changing the row name has the global effect. */

void ios_set_row_name(IOS *ios, int i, char *name)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_row_name: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_set_row_name: i = %d; row number out of range", i);
      if (name != NULL && ios_check_name(ios, name))
         fault("ios_set_row_name: i = %d; invalid name", i);
      iet_set_row_name(ios->iet, i, name);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_col_name - assign symbolic name to column.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_col_name(IOS *ios, int j, char *name);
--
-- *Description*
--
-- The routine ios_set_col_name assigns a given symbolic name to j-th
-- column of the current subproblem.
--
-- If the parameter name is NULL, the routine just erases existing name
-- of the j-th column.
--
-- NOTE: Changing the column name has the global effect. */

void ios_set_col_name(IOS *ios, int j, char *name)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_col_name: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_col_name: j = %d; column number out of range",
            j);
      if (name != NULL && ios_check_name(ios, name))
         fault("ios_set_col_name: j = %d; invalid name", j);
      iet_set_col_name(ios->iet, j, name);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_row_attr - assign attributes to row.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_row_attr(IOS *ios, int i, int mark, void *link);
--
-- *Description*
--
-- The routine ios_set_row_attr assigns two specific attributes to i-th
-- row of the current subproblem.
--
-- The first attribute is the mark, an integer, which can be used by the
-- application to distinguish between rows of different variety.
--
-- The second attribute is the link, a pointer, which can be used by the
-- application to associate additional information with the row.
--
-- NOTE: Changing the row attributes has the global effect. */

void ios_set_row_attr(IOS *ios, int i, int mark, void *link)
{     IOSRGD *rgx;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_set_row_attr: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_set_row_attr: i = %d; row number out of range", i);
      rgx = ios_get_rgd_ptr(ios, i);
      rgx->mark = mark;
      rgx->link = link;
      return;
}

/*----------------------------------------------------------------------
-- ios_set_col_attr - assign attributes to column.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_col_attr(IOS *ios, int j, int mark, void *link);
--
-- *Description*
--
-- The routine ios_set_col_attr assigns two specific attributes to j-th
-- column of the current subproblem.
--
-- The first attribute is the mark, an integer, which can be used by the
-- application to distinguish between columns of different variety.
--
-- The second attribute is the link, a pointer, which can be used by the
-- application to associate additional information with the column.
--
-- NOTE: Changing the column attributes has the global effect. */

void ios_set_col_attr(IOS *ios, int j, int mark, void *link)
{     IOSCGD *cgx;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_set_col_attr: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_col_attr: j = %d; column number out of range",
            j);
      cgx = ios_get_cgd_ptr(ios, j);
      cgx->mark = mark;
      cgx->link = link;
      return;
}

/*----------------------------------------------------------------------
-- ios_set_col_kind - set column kind.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_col_kind(IOS *ios, int j, int kind);
--
-- *Description*
--
-- The routine ios_set_col_kind sets the kind of structural variable
-- corresponding to j-th column of the current subproblem as specified
-- by the parameter kind:
--
-- LPX_NUM - continuous variable;
-- LPX_INT - integer variable (needs integer bounds).
--
-- NOTE: Changing the column kind has the global effect. */

void ios_set_col_kind(IOS *ios, int j, int kind)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_col_kind: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_col_kind: j = %d; column number out of range",
            j);
      if (!(kind == IOS_NUM || kind == IOS_INT))
         fault("ios_set_col_kind: j = %d; kind = %d; invalid column kin"
            "d", j, kind);
      if (kind == IOS_INT)
      {  int type;
         double lb, ub;
         type = ios_get_col_bnds(ios, j, &lb, &ub);
         if ((type == IOS_LO || type == IOS_DB) && lb != floor(lb))
            fault("ios_set_col_kind: j = %d; lb = %.*g; integer column "
               "needs integer lower bound", j, DBL_DIG, lb);
         if ((type == IOS_UP || type == IOS_DB) && ub != floor(ub))
            fault("ios_set_col_kind: j = %d; ub = %.*g; integer column "
               "needs integer upper bound", j, DBL_DIG, ub);
         if (type == IOS_FX && lb != floor(lb))
            fault("ios_set_col_kind: j = %d; fx = %.*g; integer column "
               "needs integer fixed value", j, DBL_DIG, lb);
      }
      ios_get_cgd_ptr(ios, j)->kind = kind;
      return;
}

/*----------------------------------------------------------------------
-- ios_set_row_bnds - set row type and bounds.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_row_bnds(IOS *ios, int i, int type, double lb,
--    double ub);
--
-- *Description*
--
-- The routine ios_set_row_bnds sets the type and bounds of i-th row of
-- the current subproblem.
--
-- The parameters type, lb, and ub specify the type, lower bound, and
-- upper bound, respectively, as shown below:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IOS_FR   -inf <  x <  +inf   free variable
--    IOS_LO     lb <= x <  +inf   lower bound
--    IOS_UP   -inf <  x <=  ub    upper bound
--    IOS_DB     lb <= x <=  ub    double bound
--    IOS_FX           x  =  lb    fixed variable
--
-- where x is auxiliary variable associated with i-th row.
--
-- If the row has no lower bound, the parameter lb is ignored. If the
-- row has no upper bound, the parameter ub is ignored. If the row is
-- of fixed type, the parameter lb is used, and the parameter ub is
-- ignored. */

void ios_set_row_bnds(IOS *ios, int i, int type, double lb, double ub)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_row_bnds: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_set_row_bnds: i = %d; row number out of range", i);
      if (!(type == IOS_FR || type == IOS_LO || type == IOS_UP ||
            type == IOS_DB || type == IOS_FX))
         fault("ios_set_row_bnds: i = %d; type = %d; invalid row type",
            i, type);
      iet_set_row_bnds(ios->iet, i, type, lb, ub);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_col_bnds - set column type and bounds.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_col_bnds(IOS *ios, int j, int type, double lb,
--    double ub);
--
-- *Description*
--
-- The routine ios_set_col_bnds sets the type and bounds of j-th column
-- of the current subproblem.
--
-- The parameters type, lb, and ub specify the type, lower bound, and
-- upper bound, respectively, as shown below:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IOS_FR   -inf <  x <  +inf   free variable
--    IOS_LO     lb <= x <  +inf   lower bound
--    IOS_UP   -inf <  x <=  ub    upper bound
--    IOS_DB     lb <= x <=  ub    double bound
--    IOS_FX           x  =  lb    fixed variable
--
-- where x is structural variable associated with j-th column.
--
-- If the column has no lower bound, the parameter lb is ignored. If
-- the column has no upper bound, the parameter ub is ignored. If the
-- column is of fixed type, the parameter lb is used, and the parameter
-- ub is ignored. */

void ios_set_col_bnds(IOS *ios, int j, int type, double lb, double ub)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_col_bnds: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_col_bnds: j = %d; column number out of range",
            j);
      if (!(type == IOS_FR || type == IOS_LO || type == IOS_UP ||
            type == IOS_DB || type == IOS_FX))
         fault("ios_set_col_bnds: j = %d; type = %d; invalid column typ"
            "e", j, type);
      if (ios_get_col_kind(ios, j) == IOS_INT)
      {  if ((type == IOS_LO || type == IOS_DB) && lb != floor(lb))
            fault("ios_set_col_bnds: j = %d; lb = %.*g; integer column "
               "needs integer lower bound", j, DBL_DIG, lb);
         if ((type == IOS_UP || type == IOS_DB) && ub != floor(ub))
            fault("ios_set_col_bnds: j = %d; ub = %.*g; integer column "
               "needs integer upper bound", j, DBL_DIG, ub);
         if (type == IOS_FX && lb != floor(lb))
            fault("ios_set_col_bnds: j = %d; fx = %.*g; integer column "
               "needs integer fixed value", j, DBL_DIG, lb);
      }
      iet_set_col_bnds(ios->iet, j, type, lb, ub);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_obj_coef - set objective coefficient or constant term.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_obj_coef(IOS *ios, int j, double coef);
--
-- *Description*
--
-- The routine ios_set_obj_coef sets the objective coefficient at j-th
-- column of the current subproblem.
--
-- If the parameter j is 0, the routine sets the constant term (shift)
-- of the objective function of the current subproblem. */

void ios_set_obj_coef(IOS *ios, int j, double coef)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_obj_coef: current subproblem does not exist");
      if (!(0 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_obj_coef: j = %d; column number out of range",
            j);
      iet_set_obj_coef(ios->iet, j, coef);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_mat_row - replace row of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_mat_row(IOS *ios, int i, int len, int ind[],
--    double val[]);
--
-- *Description*
--
-- The routine ios_set_mat_row replaces the contents of i-th row of the
-- constraint matrix of the current subproblem.
--
-- Column indices and numeric values of new row elements must be placed
-- in locations ind[1], ..., ind[len] and val[1], ..., val[len], resp.,
-- where 0 <= len <= n is new length of i-th row, and n is the number of
-- columns in the current subproblem. Note that zero elements as well as
-- elements with identical column indices are not allowed.
--
-- If the parameter len is zero, the both parameters ind and val can be
-- specified as NULL. */

void ios_set_mat_row(IOS *ios, int i, int len, int ind[], double val[])
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_mat_row: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_set_mat_row: i = %d; row number out of range", i);
      if (!(0 <= len && len <= ios_get_num_cols(ios)))
         fault("ios_set_mat_row: i = %d; len = %d; invalid row length",
            i, len);
      iet_set_mat_row(ios->iet, i, len, ind, val);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_mat_col - replace column of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_mat_col(IOS *ios, int j, int len, int ind[],
--    double val[]);
--
-- *Description*
--
-- The routine ios_set_mat_col replaces the contents of j-th column of
-- the constraint matrix of the current subproblem.
--
-- Row indices and numeric values of new column elements must be placed
-- in locations ind[1], ..., ind[len] and val[1], ..., val[len], resp.,
-- where 0 <= len <= m is new length of j-th column, and m is the number
-- of rows in the current subproblem. Note that zero elements as well as
-- elements with identical row indices are not allowed.
--
-- If the parameter len is zero, the both parameters ind and val can be
-- specified as NULL. */

void ios_set_mat_col(IOS *ios, int j, int len, int ind[], double val[])
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_mat_col: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_mat_col: j = %d; column number out of range",
            j);
      if (!(0 <= len && len <= ios_get_num_rows(ios)))
         fault("ios_set_mat_col: j = %d; len = %d; invalid column lengt"
            "h", j, len);
      iet_set_mat_col(ios->iet, j, len, ind, val);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_row_stat - set row status.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_row_stat(IOS *ios, int i, int stat);
--
-- *Description*
--
-- The routine ios_set_row_stat sets (changes) the status of i-th row
-- of the current subproblem as specified by the parameter stat:
--
-- IOS_BS   - make the row basic (make the constraint inactive);
-- IOS_NL   - make the row non-basic (make the constraint active);
-- IOS_NU   - make the row non-basic and set it to the upper bound; if
--            the row is not double-bounded, this status is equivalent
--            to IOS_NL (only for this routine);
-- IOS_NF   - the same as IOS_NL (only for this routine);
-- IOS_NS   - the same as IOS_NL (only for this routine). */

void ios_set_row_stat(IOS *ios, int i, int stat)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_row_stat: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_set_row_stat: i = %d; row number out of range", i);
      if (!(stat == IOS_BS || stat == IOS_NL || stat == IOS_NU ||
            stat == IOS_NF || stat == IOS_NS))
         fault("ios_set_row_stat: i = %d; stat = %d; invalid row status"
            , i, stat);
      iet_set_row_stat(ios->iet, i, stat);
      return;
}

/*----------------------------------------------------------------------
-- ios_set_col_stat - set column status.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_set_col_stat(IOS *ios, int j, int stat);
--
-- *Description*
--
-- The routine ios_set_col_stat sets (changes) the status of j-th column
-- of the current subproblem as specified by the parameter stat:
--
-- IOS_BS   - make the column basic;
-- IOS_NL   - make the column non-basic;
-- IOS_NU   - make the column non-basic and set it to the upper bound;
--            if the column is not of double-bounded type, this status
--            is the same as IOS_NL (only for this routine);
-- IOS_NF   - the same as IOS_NL (only for this routine);
-- IOS_NS   - the same as IOS_NL (only for this routine). */

void ios_set_col_stat(IOS *ios, int j, int stat)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_set_col_stat: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_set_col_stat: j = %d; column number out of range",
            j);
      if (!(stat == IOS_BS || stat == IOS_NL || stat == IOS_NU ||
            stat == IOS_NF || stat == IOS_NS))
         fault("ios_set_col_stat: j = %d; stat = %d; invalid column sta"
            "tus", j, stat);
      iet_set_col_stat(ios->iet, j, stat);
      return;
}

/*----------------------------------------------------------------------
-- ios_del_rows - delete specified rows from current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_del_rows(IOS *ios, int nrs, int num[]);
--
-- *Description*
--
-- The routine ios_del_rows deletes specified rows from the current
-- subproblem. Ordinal numbers of rows to be deleted should be placed
-- in locations num[1], num[2], ..., num[nrs], where nrs > 0.
--
-- Note that deleting rows involves changing ordinal numbers of other
-- rows remaining in the current subproblem. New ordinal numbers of the
-- remaining rows can be determined with the assumption that the order
-- of rows is not changed. */

void ios_del_rows(IOS *ios, int nrs, int num[])
{     int m, i, t;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_del_rows: current subproblem does not exist");
      if (nrs < 1)
         fault("ios_del_rows: nrs = %d; invalid number of rows", nrs);
      /* detach extensions from local descriptors of specified rows to
         be deleted */
      m = ios_get_num_rows(ios);
      for (t = 1; t <= nrs; t++)
      {  i = num[t];
         if (!(1 <= i && i <= m))
            fault("ios_del_rows: num[%d] = %d; row number out of range",
               t, i);
         if (ios_get_row_ptr(ios, i) == NULL)
            fault("ios_del_rows: num[%d] = %d; duplicate row numbers ar"
               "e not allowed", t, i);
         ios_detach_row(ios, i);
      }
      /* delete specified rows from the current subproblem */
      iet_del_rows(ios->iet, nrs, num);
      return;
}

/*----------------------------------------------------------------------
-- ios_del_cols - delete specified columns from current subproblem.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void ios_del_cols(IOS *ios, int ncs, int num[]);
--
-- *Description*
--
-- The routine ios_del_cols deletes specified columns from the current
-- subproblem. Ordinal numbers of columns to be deleted should be placed
-- in locations num[1], num[2], ..., num[ncs], where ncs > 0.
--
-- Note that deleting columns involves changing ordinal numbers of other
-- columns remaining in the current subproblem. New ordinal numbers of
-- the remaining columns can be determined with the assumption that the
-- order of columns is not changed. */

void ios_del_cols(IOS *ios, int ncs, int num[])
{     int n, j, t;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_del_cols: current subproblem does not exist");
      if (ncs < 1)
         fault("ios_del_cols: ncs = %d; invalid number of columns",
            ncs);
      /* detach extensions from local descriptors of specified columns
         to be deleted */
      n = ios_get_num_cols(ios);
      for (t = 1; t <= ncs; t++)
      {  j = num[t];
         if (!(1 <= j && j <= n))
            fault("ios_del_cols: num[%d] = %d; column number out of ran"
               "ge", t, j);
         if (ios_get_col_ptr(ios, j) == NULL)
            fault("ios_del_cols: num[%d] = %d; duplicate column indices"
               " are not allowed", t, j);
         ios_detach_col(ios, j);
      }
      /* delete specified columns from the current subproblem */
      iet_del_cols(ios->iet, ncs, num);
      return;
}

/**********************************************************************/
/* * *                SUBPROBLEM QUERYING ROUTINES                * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_get_obj_dir - determine optimization direction flag.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_obj_dir(IOS *ios);
--
-- *Returns*
--
-- The routine ios_get_obj_dir returns the optimization direction flag
-- as follows:
--
-- IOS_MIN means minimization;
-- IOS_MAX means maximization. */

int ios_get_obj_dir(IOS *ios)
{     int dir = ios->dir;
      return dir;
}

/*----------------------------------------------------------------------
-- ios_get_num_rows - determine number of rows.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_num_rows(IOS *ios);
--
-- *Returns*
--
-- The routine ios_get_num_rows returns the number of rows in the
-- current subproblem. */

int ios_get_num_rows(IOS *ios)
{     int m;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_num_rows: current subproblem does not exist");
      m = iet_get_num_rows(ios->iet);
      return m;
}

/*----------------------------------------------------------------------
-- ios_get_num_cols - determine number of columns.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_num_cols(IOS *ios);
--
-- *Returns*
--
-- The routine ios_get_num_cols returns the number of columns in the
-- current subproblem. */

int ios_get_num_cols(IOS *ios)
{     int n;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_num_cols: current subproblem does not exist");
      n = iet_get_num_cols(ios->iet);
      return n;
}

/*----------------------------------------------------------------------
-- ios_get_num_nz - determine number of constraint coefficients.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_num_nz(IOS *ios);
--
-- *Returns*
--
-- The routine ios_get_num_nz returns the number of (non-zero) elements
-- in the constraint matrix of the current subproblem. */

int ios_get_num_nz(IOS *ios)
{     int nz;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_num_nz: current subproblem does not exist");
      nz = iet_get_num_nz(ios->iet);
      return nz;
}

/*----------------------------------------------------------------------
-- ios_get_row_name - obtain row name.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- char *ios_get_row_name(IOS *ios, int i);
--
-- *Returns*
--
-- The routine ios_get_row_name returns a pointer to the symbolic name
-- assigned to i-th row of the current subproblem. However, if the row
-- has no symbolic name assigned, the routine returns NULL.
--
-- If this routine is called on processing the event IOS_V_DELROW, the
-- parameter i must be zero, in which case a pointer to the name of the
-- row being deleted is returned. */

char *ios_get_row_name(IOS *ios, int i)
{     char *name;
      if (i == 0 && ios->event == IOS_V_DELROW)
      {  name = ios->hook_name;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_name: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_name: i = %d; row number out of range", i);
      name = iet_get_row_name(ios->iet, i);
done: return name;
}

/*----------------------------------------------------------------------
-- ios_get_col_name - obtain column name.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- char *ios_get_col_name(IOS *ios, int j);
--
-- *Returns*
--
-- The routine ios_get_col_name returns a pointer to the symbolic name
-- assigned to j-th column of the current subproblem. However, if the
-- column has no symbolic name assigned, the routine returns NULL.
--
-- If this routine is called on processing the event IOS_V_DELCOL, the
-- parameter j must be zero, in which case a pointer to the name of the
-- column being deleted is returned. */

char *ios_get_col_name(IOS *ios, int j)
{     char *name;
      if (j == 0 && ios->event == IOS_V_DELCOL)
      {  name = ios->hook_name;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_name: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_name: j = %d; column number out of range",
            j);
      name = iet_get_col_name(ios->iet, j);
done: return name;
}

/*----------------------------------------------------------------------
-- ios_get_row_mark - obtain row mark.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_row_mark(IOS *ios, int i);
--
-- *Returns*
--
-- The routine ios_get_row_mark returns the mark (application specific
-- attribute) assigned to i-th row of the current subproblem.
--
-- If this routine is called on processing the event IOS_V_DELROW, the
-- parameter i must be zero, in which case the mark attribute assigned
-- to the row being deleted is returned. */

int ios_get_row_mark(IOS *ios, int i)
{     int mark;
      if (i == 0 && ios->event == IOS_V_DELROW)
      {  insist(ios->hook_link.rgd != NULL);
         mark = ios->hook_link.rgd->mark;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_mark: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_mark: i = %d; row number out of range", i);
      mark = ios_get_rgd_ptr(ios, i)->mark;
done: return mark;
}

/*----------------------------------------------------------------------
-- ios_get_row_link - obtain row link.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void *ios_get_row_link(IOS *ios, int i);
--
-- *Returns*
--
-- The routine ios_get_row_link returns the link (application specific
-- attribute) assigned to i-th row of the current subproblem.
--
-- If this routine is called on processing the event IOS_V_DELROW, the
-- parameter i must be zero, in which case the mark attribute assigned
-- to the row being deleted is returned. */

void *ios_get_row_link(IOS *ios, int i)
{     void *link;
      if (i == 0 && ios->event == IOS_V_DELROW)
      {  insist(ios->hook_link.rgd != NULL);
         link = ios->hook_link.rgd->link;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_link: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_link: i = %d; row number out of range", i);
      link = ios_get_rgd_ptr(ios, i)->link;
done: return link;
}

/*----------------------------------------------------------------------
-- ios_get_col_mark - obtain column mark.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_col_mark(IOS *ios, int j);
--
-- *Returns*
--
-- The routine ios_get_col_mark returns the mark (application specific
-- attribute) assigned to j-th column of the current subproblem.
--
-- If this routine is called on processing the event IOS_V_DELCOL, the
-- parameter j must be zero, in which case the mark attribute assigned
-- to the column being deleted is returned. */

int ios_get_col_mark(IOS *ios, int j)
{     int mark;
      if (j == 0 && ios->event == IOS_V_DELCOL)
      {  insist(ios->hook_link.cgd != NULL);
         mark = ios->hook_link.cgd->mark;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_mark: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_mark: j = %d; column number out of range",
            j);
      mark = ios_get_cgd_ptr(ios, j)->mark;
done: return mark;
}

/*----------------------------------------------------------------------
-- ios_get_col_link - obtain column link.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- void *ios_get_col_link(IOS *ios, int j);
--
-- *Returns*
--
-- The routine ios_get_col_link returns the link (application specific
-- attribute) assigned to j-th column of the current subproblem.
--
-- If this routine is called on processing the event IOS_V_DELCOL, the
-- parameter j must be zero, in which case the link attribute assigned
-- to the column being deleted is returned. */

void *ios_get_col_link(IOS *ios, int j)
{     void *link;
      if (j == 0 && ios->event == IOS_V_DELCOL)
      {  insist(ios->hook_link.cgd != NULL);
         link = ios->hook_link.cgd->link;
         goto done;
      }
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_link: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_link: j = %d; column number out of range",
            j);
      link = ios_get_cgd_ptr(ios, j)->link;
done: return link;
}

/*----------------------------------------------------------------------
-- ios_get_col_kind - determine column kind.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_col_kind(IOS *ios, int j);
--
-- *Returns*
--
-- The routine ios_get_col_kind returns the kind of structural variable
-- corresponding to j-th column of the current subproblem:
--
-- IOS_NUM - continuous variable;
-- IOS_INT - integer variable. */

int ios_get_col_kind(IOS *ios, int j)
{     int kind;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_kind: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_kind: j = %d; column number out of range",
            j);
      kind = ios_get_cgd_ptr(ios, j)->kind;
      return kind;
}

/*----------------------------------------------------------------------
-- ios_get_row_bnds - determine row type and bounds.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_row_bnds(IOS *ios, int i, double *lb, double *ub);
--
-- *Description*
--
-- The routine ios_get_row_bnds determines the type, lower bound, and
-- upper bound of i-th row of the current subproblem.
--
-- The lower and upper bounds are stored to locations specified by the
-- parameters lb and ub, respectively. If some of these parameters is
-- NULL, the corresponding value is not stored. The type is returned by
-- the routine on exit.
--
-- The type and bounds have the following meaning:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IOS_FR   -inf <  x <  +inf   free variable
--    IOS_LO     lb <= x <  +inf   lower bound
--    IOS_UP   -inf <  x <=  ub    upper bound
--    IOS_DB     lb <= x <=  ub    double bound
--    IOS_FX           x  =  lb    fixed variable
--
-- where x is the auxiliary variable associated with i-th row.
--
-- If the row has no lower bound, *lb is set to zero. If the row has no
-- upper bound, *ub is set to zero. If the row is of fixed type, *lb and
-- *ub are set to the same value.
--
-- *Returns*
--
-- The routine ios_get_row_bnds returns the type of the row. */

int ios_get_row_bnds(IOS *ios, int i, double *lb, double *ub)
{     int type;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_bnds: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_bnds: i = %d; row number out of range", i);
      type = iet_get_row_bnds(ios->iet, i, lb, ub);
      return type;
}

/*----------------------------------------------------------------------
-- ios_get_col_bnds - determine column type and bounds.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_col_bnds(IOS *ios, int j, double *lb, double *ub);
--
-- *Description*
--
-- The routine ios_get_col_bnds determines the type, lower bound, and
-- upper bound of j-th column of the current subproblem.
--
-- The lower and upper bounds are stored to locations specified by the
-- parameters lb and ub, respectively. If some of these parameters is
-- NULL, the corresponding value is not stored. The type is returned by
-- the routine on exit.
--
-- The type and bounds have the following meaning:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IOS_FR   -inf <  x <  +inf   free variable
--    IOS_LO     lb <= x <  +inf   lower bound
--    IOS_UP   -inf <  x <=  ub    upper bound
--    IOS_DB     lb <= x <=  ub    double bound
--    IOS_FX           x  =  lb    fixed variable
--
-- where x is the structural variable associated with j-th column.
--
-- If the column has no lower bound, *lb is set to zero. If the column
-- has no upper bound, *ub is set to zero. If the column is of fixed
-- type, *lb and *ub are set to the same value.
--
-- *Returns*
--
-- The routine ios_get_col_bnds returns the type of the column. */

int ios_get_col_bnds(IOS *ios, int j, double *lb, double *ub)
{     int type;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_bnds: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_bnds: j = %d; column number out of range",
            j);
      type = iet_get_col_bnds(ios->iet, j, lb, ub);
      return type;
}

/*----------------------------------------------------------------------
-- ios_get_obj_coef - determine objective coefficient.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- double ios_get_obj_coef(IOS *ios, int j);
--
-- *Returns*
--
-- The routine ios_get_obj_coef returns objective coefficient at j-th
-- column of the current subproblem. If j is 0, the routine returns
-- constant term of the objective function of the current subproblem. */

double ios_get_obj_coef(IOS *ios, int j)
{     double coef;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_obj_coef: current subproblem does not exist");
      if (!(0 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_obj_coef: j = %d; column number out of range",
            j);
      coef = iet_get_obj_coef(ios->iet, j);
      return coef;
}

/*----------------------------------------------------------------------
-- ios_get_mat_row - obtain row of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_mat_row(IOS *ios, int i, int ind[], double val[]);
--
-- *Description*
--
-- The routine ios_get_mat_row scans the contents of i-th row of the
-- constraint matrix of the current subproblem and stores column indices
-- and numeric values of corresponding (non-zero) elements to locations
-- ind[1], ..., ind[len] and val[1], ..., val[len], respectively, where
-- 0 <= len <= n is the number of non-zero elements in i-th row, n is
-- the number of columns in the current subproblem. If the parameter ind
-- or val is NULL, the corresponding information is not stored.
--
-- *Returns*
--
-- The routine ios_get_mat_row returns the number of non-zero elements
-- in i-th row (len). */

int ios_get_mat_row(IOS *ios, int i, int ind[], double val[])
{     int len;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_mat_row: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_mat_row: i = %d; row number out of range", i);
      len = iet_get_mat_row(ios->iet, i, ind, val);
      return len;
}

/*----------------------------------------------------------------------
-- ios_get_mat_col - obtain column of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_mat_col(IOS *ios, int j, int ind[], double val[]);
--
-- *Description*
--
-- The routine ios_get_mat_col scans the contents of j-th column of the
-- constraint matrix of the current subproblem and stores row indices
-- and numeric values of corresponding (non-zero) elements to locations
-- ind[1], ..., ind[len] and val[1], ..., val[len], respectively, where
-- 0 <= len <= m is the number of non-zero elements in j-th column, m is
-- the number of rows in the current subproblem. If the parameter ind or
-- val is NULL, the corresponding information is not stored.
--
-- The routine ios_get_mat_col returns the number of non-zero elements
-- in j-th column (len). */

int ios_get_mat_col(IOS *ios, int j, int ind[], double val[])
{     int len;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_mat_col: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_mat_col: j = %d; column number out of range",
            j);
      len = iet_get_mat_col(ios->iet, j, ind, val);
      return len;
}

/**********************************************************************/
/* * *              BASIC SOLUTION QUERYING ROUTINES              * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- ios_p_status - determine primal status of basic solution.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_p_status(IOS *ios);
--
-- *Returns*
--
-- The routine ios_p_status reports the primal status of basic solution
-- to LP relaxaton of the current subproblem:
--
-- IOS_UNDEF  - primal status is undefined;
-- IOS_FEAS   - solution is primal feasible;
-- IOS_INFEAS - solution is primal infeasible;
-- IOS_NOFEAS - LP relaxation has no primal feasible solution. */

int ios_p_status(IOS *ios)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_p_status: current subproblem does not exist");
      return ios->p_stat;
}

/*----------------------------------------------------------------------
-- ios_d_status - determine dual status of basic solution.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_d_status(IOS *ios);
--
-- *Returns*
--
-- The routine ios_d_status reports the dual status of basic solution
-- to LP relaxaton of the current subproblem:
--
-- IOS_UNDEF  - dual status is undefined;
-- IOS_FEAS   - solution is dual feasible;
-- IOS_INFEAS - solution is dual infeasible;
-- IOS_NOFEAS - LP relaxation has no dual feasible solution. */

int ios_d_status(IOS *ios)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_d_status: current subproblem does not exist");
      return ios->d_stat;
}

/*----------------------------------------------------------------------
-- ios_get_row_soln - obtain basic solution for given row.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_row_soln(IOS *ios, int i, double *prim, double *dual);
--
-- *Description*
--
-- The routine ios_get_row_soln obtains basic solution information for
-- i-th row of the current subproblem.
--
-- The primal and dual values of the auxiliary variable are stored in
-- locations specified by the parameters prim and dual, respectively.
-- If some of these parameters is NULL, the corresponding value is not
-- stored.
--
-- *Returns*
--
-- The routine ios_get_row_soln returns the status of i-th row:
--
-- IET_BS - basic row (inactive constraint);
-- IET_NL - non-basic row on its lower bound (active constraint);
-- IET_NU - non-basic row on its upper bound (active constraint);
-- IET_NF - non-basic free (unbounded) row;
-- IET_NS - non-basic fixed row (equality constraint). */

int ios_get_row_soln(IOS *ios, int i, double *prim, double *dual)
{     IOSROW *row;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_soln: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_soln: i = %d; row number out of range", i);
      row = ios_get_row_ptr(ios, i);
      if (prim != NULL) *prim = row->prim;
      if (dual != NULL) *dual = row->dual;
      return iet_get_row_stat(ios->iet, i);
}

/*----------------------------------------------------------------------
-- ios_get_col_soln - obtain basic solution for given column.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_get_col_soln(IOS *ios, int j, double *prim, double *dual);
--
-- *Description*
--
-- The routine ios_get_col_soln obtains basic solution information for
-- j-th column of the current subproblem.
--
-- The primal and dual values of the structural variable are stored in
-- locations specified by the parameters prim and dual, respectively.
-- If some of these parameters is NULL, the corresponding value is not
-- stored.
--
-- *Returns*
--
-- The routine ios_get_col_soln returns the status of j-th column:
--
-- IET_BS - basic column;
-- IET_NL - non-basic column on its lower bound;
-- IET_NU - non-basic column on its upper bound;
-- IET_NF - non-basic free (unbounded) column;
-- IET_NS - non-basic fixed column. */

int ios_get_col_soln(IOS *ios, int j, double *prim, double *dual)
{     IOSCOL *col;
      if (ios_get_curr_node(ios) == 0)
         fault("ios_get_col_soln: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_get_col_soln: j = %d; column number out of range",
            j);
      col = ios_get_col_ptr(ios, j);
      if (prim != NULL) *prim = col->prim;
      if (dual != NULL) *dual = col->dual;
      return iet_get_col_stat(ios->iet, j);
}

/*----------------------------------------------------------------------
-- ios_get_row_pi - determine Lagrange multiplier for given row.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- double ios_get_row_pi(IOS *ios, int i);
--
-- *Returns*
--
-- The routine ios_get_row_pi returns Lagrange multiplier for i-th row
-- (constraint) needed to generate columns.
--
-- *Comments*
--
-- Using Lagrange multipliers returned by the routine ios_get_row_pi
-- the application procedure can compute reduced costs of columns which
-- are missing in the current subproblem as follows:
--
--    d = (if p_status is IOS_FEAS then c else 0) +
--                                                                   (1)
--      + a[1]*pi[1] + a[2]*pi[2] + ... + a[m]*pi[m],
--
-- where d is the reduced cost to be computed, p_status is the primal
-- status of basic solution to LP relaxation as reported by the routine
-- ios_p_status, c is the objective coefficient, a[1], ..., a[m] are
-- constraint coeffcients for rows of LP relaxation, pi[1], ..., pi[m]
-- are Lagrange multipliers.
--
-- If the basic solution is optimal (i.e. primal and dual feasible),
-- Lagrange multipliers correspond to the original objective function.
-- If LP relaxation has no primal feasible solution, the multipliers
-- correspond to the sum of primal infeasibilities. In the latter case,
-- if the original objective function has to be maximized, the sum of
-- primal infeasibilities is taken with the minus sign in order to keep
-- the original objective sense. */

double ios_get_row_pi(IOS *ios, int i)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_get_row_pi: current subproblem does not exist");
      if (!(1 <= i && i <= ios_get_num_rows(ios)))
         fault("ios_get_row_pi: i = %d; row number out of range", i);
      return ios_get_row_ptr(ios, i)->pi;
}

/*----------------------------------------------------------------------
-- ios_is_col_frac - check if specified column has fractional value.
--
-- *Synopsis*
--
-- #include "glpios.h"
-- int ios_is_col_frac(IOS *ios, int j);
--
-- *Returns*
--
-- If j-th column of the current subproblem is of integer kind and has
-- a fractional primal value in basic solution of the LP relaxation, the
-- routine returns non-zero. Otherwise, if the column is continuous or
-- its value is integral within a given tolerance, zero is returned. */

int ios_is_col_frac(IOS *ios, int j)
{     if (ios_get_curr_node(ios) == 0)
         fault("ios_is_col_frac: current subproblem does not exist");
      if (!(1 <= j && j <= ios_get_num_cols(ios)))
         fault("ios_is_col_frac: j = %d; column number out of range",
            j);
      return ios_get_col_ptr(ios, j)->frac;
}

/* eof */

/* glpiet.c */

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

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "glpiet.h"
#include "glplib.h"

/**********************************************************************/
/* * *                  TREE MANAGEMENT ROUTINES                  * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- iet_create_tree - create implicit enumeration tree.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- IET *iet_create_tree(void);
--
-- The routine iet_create_tree creates the implicit enumeration tree.
-- Being created the tree consists of the only root subproblem which is
-- empty, i.e. has no rows and columns. The reference number assigned
-- to the root subproblem is always 1. Note that the root subproblem is
-- initially in the frozen state and therefore needs to be revived.
--
-- *Returns*
--
-- The routine returns a pointer to the tree created. This pointer must
-- be used in all subsequent operations. */

IET *iet_create_tree(void)
{     IET *iet;
      IETNPD *node;
      int p;
      /* create implicit enumeration tree */
      iet = umalloc(sizeof(IET));
      iet->npd_pool = dmp_create_pool(sizeof(IETNPD));
      iet->rgd_pool = dmp_create_pool(sizeof(IETRGD));
      iet->cgd_pool = dmp_create_pool(sizeof(IETCGD));
      iet->dqe_pool = dmp_create_pool(sizeof(IETDQE));
      iet->bqe_pool = dmp_create_pool(sizeof(IETBQE));
      iet->cqe_pool = dmp_create_pool(sizeof(IETCQE));
      iet->aqe_pool = dmp_create_pool(sizeof(IETAQE));
      iet->aij_pool = dmp_create_pool(sizeof(IETAIJ));
      iet->sqe_pool = dmp_create_pool(sizeof(IETSQE));
      iet->row_pool = dmp_create_pool(sizeof(IETROW));
      iet->col_pool = dmp_create_pool(sizeof(IETCOL));
      iet->str_pool = create_str_pool();
      iet->str_buf = ucalloc(255+1, sizeof(char));
      iet->nslots = 20;
      iet->avail = 0;
      iet->slot = ucalloc(1+iet->nslots, sizeof(IETNPS));
      iet->head = NULL;
      iet->tail = NULL;
      iet->a_cnt = 0;
      iet->n_cnt = 0;
      iet->t_cnt = 0;
      iet->hook = NULL;
      iet->info = NULL;
      iet->curr = NULL;
      iet->m_max = 50;
      iet->n_max = 100;
      iet->m = 0;
      iet->n = 0;
      iet->nz = 0;
      iet->c0 = 0.0;
      iet->old_c0 = 0.0;
      iet->row = ucalloc(1+iet->m_max, sizeof(IETROW *));
      iet->col = ucalloc(1+iet->n_max, sizeof(IETCOL *));
      /* initialize the stack of free slots */
      for (p = iet->nslots; p >= 1; p--)
      {  iet->slot[p].node = NULL;
         iet->slot[p].next = iet->avail;
         iet->avail = p;
      }
      /* pull a free slot for the root node */
      p = iet->avail;
      insist(p == 1);
      iet->avail = iet->slot[p].next;
      insist(iet->slot[p].node == NULL);
      iet->slot[p].next = 0;
      /* create the root subproblem */
      iet->slot[p].node = node = dmp_get_atom(iet->npd_pool);
      node->p = p;
      node->up = NULL;
      node->level = 0;
      node->count = 0;
      node->r_add = NULL;
      node->c_add = NULL;
      node->r_del = NULL;
      node->c_del = NULL;
      node->r_bnds = NULL;
      node->c_bnds = NULL;
      node->c_obj = NULL;
      node->r_mat = NULL;
      node->c_mat = NULL;
      node->r_stat = NULL;
      node->c_stat = NULL;
      node->link = NULL;
      node->temp = NULL;
      node->prev = NULL;
      node->next = NULL;
      /* add the root subproblem to the active list */
      iet->head = iet->tail = node;
      iet->a_cnt++;
      iet->n_cnt++;
      iet->t_cnt++;
      return iet;
}

/*----------------------------------------------------------------------
-- iet_install_hook - install higher-level hook routine.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_install_hook(IET *iet, void (*hook)(void *info, int what,
--    char *name, void *link), void *info);
--
-- *Description*
--
-- The routine iet_install_hook installs a higher-level hook routine.
--
-- The hook routine is called whenever a subproblem, row, or column is
-- being deleted from the tree. The purpose of the hook routine is to
-- delete additional higher-level information which might be associated
-- with corresponding object.
--
-- The parameter info is a transitional pointer which is passed to the
-- hook routine.
--
-- The parameter what specifies what object is being deleted:
--
-- IOS_ND - subproblem;
-- IOS_RD - row;
-- IOS_CD - column.
--
-- If a subproblem is being deleted, the parameter name is NULL. If a
-- row or column is being deleted, the parameter name is the symbolic
-- name assigned to corresponding object.
--
-- The parameter link is the link to a global higher-level extension of
-- corresponding object. */

void iet_install_hook(IET *iet, void (*hook)(void *info, int what,
      char *name, void *link), void *info)
{     iet->hook = hook;
      iet->info = info;
      return;
}

/*----------------------------------------------------------------------
-- iet_revive_node - revive specified subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_revive_node(IET *iet, int p);
--
-- *Description*
--
-- The routine iet_revive_node revives the specified subproblem, whose
-- reference number is p, and thereby makes it the current subproblem.
-- Note that the specified subproblem must be active. Besides, if the
-- current subproblem already exists, it must be frozen before reviving
-- another subproblem. */

void iet_revive_node(IET *iet, int p)
{     IETNPD *node, *root;
      IETRGD *rgd, *r_head, *r_tail;
      IETCGD *cgd, *c_head, *c_tail;
      IETDQE *dqe;
      IETBQE *bqe;
      IETCQE *cqe;
      IETAQE *aqe;
      IETAIJ *aij;
      IETSQE *sqe;
      IETROW *row;
      IETCOL *col;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_revive_node: p = %d; invalid subproblem reference n"
            "umber", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* the specified subproblem must be active */
      if (node->count != 0)
         fault("iet_revive_node: p = %d; reviving inactive subproblem n"
            "ot allowed", p);
      /* the current subproblem must not exist */
      if (iet->curr != NULL)
         fault("iet_revive_node: current subproblem already exists");
      /* the specified subproblem becomes current */
      iet->curr = node;
      insist(iet->m == 0);
      insist(iet->n == 0);
      insist(iet->nz == 0);
      iet->c0 = 0.0;
      iet->old_c0 = iet->c0;
      /* obtain pointer to the root subproblem */
      root = iet->slot[1].node;
      insist(root != NULL);
      /* build the path from the root to the given node */
      node->temp = NULL;
      for (node = node; node != NULL; node = node->up)
      {  if (node->up == NULL)
            insist(node == root);
         else
            node->up->temp = node;
      }
      /* walk from the root to the given node and build the ordered
         linked lists of rows and columns which have been introduced in
         each subproblem */
      r_head = r_tail = NULL; /* the list of rows */
      c_head = c_tail = NULL; /* the list of columns */
      for (node = root; node != NULL; node = node->temp)
      {  /* add own rows of this subproblem to the row list */
         for (rgd = node->r_add; rgd != NULL; rgd = rgd->next)
         {  /* check that the row is not yet in the row list */
            insist(rgd->i == 0);
            /* mark the row that now it is there */
            rgd->i = +1;
            /* and add it to the end of the row list */
            rgd->temp = NULL;
            if (r_head == NULL)
               r_head = rgd;
            else
               r_tail->temp = rgd;
            r_tail = rgd;
         }
         /* add own columns of this subproblem to the column list */
         for (cgd = node->c_add; cgd != NULL; cgd = cgd->next)
         {  /* check that the column is not yet in the column list */
            insist(cgd->j == 0);
            /* mark the column that now it is there */
            cgd->j = +1;
            /* and add it to the end of the column list */
            cgd->temp = NULL;
            if (c_head == NULL)
               c_head = cgd;
            else
               c_tail->temp = cgd;
            c_tail = cgd;
         }
         /* mark rows which have been deleted from this subproblem */
         for (dqe = node->r_del; dqe != NULL; dqe = dqe->next)
         {  /* obtain pointer to global descriptor of the row */
            rgd = dqe->u.row;
            /* this can be only non-own row, so it must be already in
               the row list */
            insist(rgd->i == +1);
            /* mark the row to be deleted */
            rgd->i = -1;
         }
         /* mark columns which have been deleted from this subproblem */
         for (dqe = node->c_del; dqe != NULL; dqe = dqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = dqe->u.col;
            /* this can be only non-own column, so it must be already
               in the column list */
            insist(cgd->j == +1);
            /* mark the column to be deleted */
            cgd->j = -1;
         }
      }
      /* build the set of rows of the subproblem to be revived */
      insist(iet->m == 0);
      for (rgd = r_head; rgd != NULL; rgd = rgd->temp)
      {  /* if the row is marked to be deleted, skip it */
         if (rgd->i == -1)
         {  rgd->i = 0;
            continue;
         }
         /* assign ordinal number to the row */
         rgd->i = ++(iet->m);
         insist(rgd->i <= iet->m_max);
         /* create local descriptor for the row */
         iet->row[rgd->i] = row = dmp_get_atom(iet->row_pool);
         row->glob = rgd;
         row->type = IET_FR;
         row->lb = 0.0;
         row->ub = 0.0;
         row->set_by = rgd->host;
         row->ptr = NULL;
         row->stat = IET_BS;
         row->old_type = row->type;
         row->old_lb = row->lb;
         row->old_ub = row->ub;
         row->old_stat = row->stat;
         row->link = NULL;
      }
      /* build the set of columns of the subproblem to be revived */
      insist(iet->n == 0);
      for (cgd = c_head; cgd != NULL; cgd = cgd->temp)
      {  /* if the column is marked to be deleted, skip it */
         if (cgd->j == -1)
         {  cgd->j = 0;
            continue;
         }
         /* assign ordinal number to the column */
         cgd->j = ++(iet->n);
         insist(cgd->j <= iet->n_max);
         /* create local descriptor for the column */
         iet->col[cgd->j] = col = dmp_get_atom(iet->col_pool);
         col->glob = cgd;
         col->type = IET_FX;
         col->lb = 0.0;
         col->ub = 0.0;
         col->coef = 0.0;
         col->set_by = cgd->host;
         col->ptr = NULL;
         col->stat = IET_NS;
         col->old_type = col->type;
         col->old_lb = col->lb;
         col->old_ub = col->ub;
         col->old_coef = col->coef;
         col->old_stat = col->stat;
         col->link = NULL;
      }
      /* walk from the root to the given node and update the attribute
         set_by for rows and columns of the subproblem to be revived */
      for (node = root; node != NULL; node = node->temp)
      {  /* go through the list of rows which have been replaced by
            iet_set_mat_row in this subproblem */
         for (aqe = node->r_mat; aqe != NULL; aqe = aqe->next)
         {  /* obtain pointer to global descriptor of the row */
            rgd = aqe->u.row;
            /* skip the row missing in the revived subproblem */
            if (rgd->i == 0) continue;
            /* change the attribute set_by, if necessary */
            row = iet->row[rgd->i];
            if (row->set_by->level < node->level) row->set_by = node;
         }
         /* go through the list of columns which have been replaced by
            iet_set_mat_col in this subproblem */
         for (aqe = node->c_mat; aqe != NULL; aqe = aqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = aqe->u.col;
            /* skip the column missing in the revived subproblem */
            if (cgd->j == 0) continue;
            /* change the attribute set_by, if necessary */
            col = iet->col[cgd->j];
            if (col->set_by->level < node->level) col->set_by = node;
         }
      }
      /* walk from the root to the revived node and restore all other
         attributes of rows and columns of the revived subproblem */
      for (node = root; node != NULL; node = node->temp)
      {  /* if the given node has been reached, save attributes of rows
            and columns which currently correspond to the parent of the
            revived subproblem */
         if (node->temp == NULL)
         {  int i, j;
            iet->old_c0 = iet->c0;
            for (i = 1; i <= iet->m; i++)
            {  row = iet->row[i];
               row->old_type = row->type;
               row->old_lb = row->lb;
               row->old_ub = row->ub;
               row->old_stat = row->stat;
            }
            for (j = 1; j <= iet->n; j++)
            {  col = iet->col[j];
               col->old_type = col->type;
               col->old_lb = col->lb;
               col->old_ub = col->ub;
               col->old_coef = col->coef;
               col->old_stat = col->stat;
            }
         }
         /* restore types and bounds of rows */
         for (bqe = node->r_bnds; bqe != NULL; bqe = bqe->next)
         {  /* obtain pointer to global descriptor of the row */
            rgd = bqe->u.row;
            /* skip the row missing in the revived subproblem */
            if (rgd->i == 0) continue;
            /* obtain pointer to local descriptor of the row */
            row = iet->row[rgd->i];
            /* set new type and bounds of the row */
            row->type = bqe->type;
            row->lb = bqe->lb;
            row->ub = bqe->ub;
         }
         /* restore type and bounds of columns */
         for (bqe = node->c_bnds; bqe != NULL; bqe = bqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = bqe->u.col;
            /* skip the column missing in the revived subproblem */
            if (cgd->j == 0) continue;
            /* obtain pointer to local descriptor of the column */
            col = iet->col[cgd->j];
            /* set new type and bounds of the column */
            col->type = bqe->type;
            col->lb = bqe->lb;
            col->ub = bqe->ub;
         }
         /* restore objective coefficients of columns */
         for (cqe = node->c_obj; cqe != NULL; cqe = cqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = cqe->col;
            /* null pointer means new constant term of the objective
               function */
            if (cgd == NULL)
            {  iet->c0 = cqe->coef;
               continue;
            }
            /* skip the column missing in the revived subproblem */
            if (cgd->j == 0) continue;
            /* obtain pointer to local descriptor of the column */
            col = iet->col[cgd->j];
            /* set new objective coefficient */
            col->coef = cqe->coef;
         }
         /* restore rows of the constraint matrix */
         for (aqe = node->r_mat; aqe != NULL; aqe = aqe->next)
         {  /* obtain pointer to global descriptor of the row */
            rgd = aqe->u.row;
            /* skip the row missing in the current subproblem */
            if (rgd->i == 0) continue;
            /* obtain pointer to local descriptor of the row */
            row = iet->row[rgd->i];
            /* skip the row if it has been replaced in a subsequent
               subproblem of higher level */
            if (row->set_by->level > node->level) continue;
            /* add constraint coefficients to the row */
            for (aij = aqe->ptr; aij != NULL; aij = aij->link)
            {  insist(aij->row == rgd);
               /* obtain pointer to global descriptor of corresponding
                  column */
               cgd = aij->col;
               /* skip the coefficient if it is placed in the column
                  which is missing in the revived subproblem */
               if (cgd->j == 0) continue;
               /* obtain pointer to local descriptor of the column */
               col = iet->col[cgd->j];
               /* skip the coefficient if it is placed in the column
                  which has been replaced in another subproblem */
               if (col->set_by->level > node->level) continue;
               /* add the coefficient to the constraint matrix */
               aij->r_prev = NULL;
               aij->r_next = row->ptr;
               aij->c_prev = NULL;
               aij->c_next = col->ptr;
               if (row->ptr != NULL) row->ptr->r_prev = aij;
               if (col->ptr != NULL) col->ptr->c_prev = aij;
               row->ptr = col->ptr = aij;
               /* increase the non-zero count */
               iet->nz++;
            }
         }
         /* restore columns of the constraint matrix */
         for (aqe = node->c_mat; aqe != NULL; aqe = aqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = aqe->u.col;
            /* skip the column missing in the revived subproblem */
            if (cgd->j == 0) continue;
            /* obtain pointer to local descriptor of the column */
            col = iet->col[cgd->j];
            /* skip the column if it has been replaced in a subsequent
               subproblem of higher level */
            if (col->set_by->level > node->level) continue;
            /* add constraint coefficients to the column */
            for (aij = aqe->ptr; aij != NULL; aij = aij->link)
            {  insist(aij->col == cgd);
               /* obtain pointer to global descriptor of corresponding
                  row */
               rgd = aij->row;
               /* skip the coefficient if it is placed in the row which
                  is missing in the revived subproblem */
               if (rgd->i == 0) continue;
               /* obtain pointer to local descriptor of the row */
               row = iet->row[rgd->i];
               /* skip the coefficient if it is placed in the row which
                  has been replaced in another subproblem */
               if (row->set_by->level > node->level) continue;
               /* add the coefficient to the constraint matrix */
               aij->r_prev = NULL;
               aij->r_next = row->ptr;
               aij->c_prev = NULL;
               aij->c_next = col->ptr;
               if (row->ptr != NULL) row->ptr->r_prev = aij;
               if (col->ptr != NULL) col->ptr->c_prev = aij;
               row->ptr = col->ptr = aij;
               /* increase the non-zero count */
               iet->nz++;
            }
         }
         /* restore statuses of rows */
         for (sqe = node->r_stat; sqe != NULL; sqe = sqe->next)
         {  /* obtain pointer to global descriptor of the row */
            rgd = sqe->u.row;
            /* skip the row missing in the revived subproblem */
            if (rgd->i == 0) continue;
            /* obtain pointer to local descriptor of the row */
            row = iet->row[rgd->i];
            /* set new status of the row */
            row->stat = sqe->stat;
         }
         /* restore statuses of columns */
         for (sqe = node->c_stat; sqe != NULL; sqe = sqe->next)
         {  /* obtain pointer to global descriptor of the column */
            cgd = sqe->u.col;
            /* skip the column missing in the revived subproblem */
            if (cgd->j == 0) continue;
            /* obtain pointer to local descriptor of the column */
            col = iet->col[cgd->j];
            /* set new status of the column */
            col->stat = sqe->stat;
         }
      }
      /* the specified subproblem has been revived; delete its change
         lists */
      node = iet->curr;
      /* delete the row types/bounds change list */
      while (node->r_bnds != NULL)
      {  bqe = node->r_bnds;
         node->r_bnds = bqe->next;
         dmp_free_atom(iet->bqe_pool, bqe);
      }
      /* delete the column types/bounds change list */
      while (node->c_bnds != NULL)
      {  bqe = node->c_bnds;
         node->c_bnds = bqe->next;
         dmp_free_atom(iet->bqe_pool, bqe);
      }
      /* delete the column objective coefficients change list */
      while (node->c_obj != NULL)
      {  cqe = node->c_obj;
         node->c_obj = cqe->next;
         dmp_free_atom(iet->cqe_pool, cqe);
      }
      /* delete the change list for rows of the constraint matrix */
      while (node->r_mat != NULL)
      {  aqe = node->r_mat;
         node->r_mat = aqe->next;
         dmp_free_atom(iet->aqe_pool, aqe);
      }
      /* delete the change list for columns of the constraint matrix */
      while (node->c_mat != NULL)
      {  aqe = node->c_mat;
         node->c_mat = aqe->next;
         dmp_free_atom(iet->aqe_pool, aqe);
      }
      /* delete the row statuses change list */
      while (node->r_stat != NULL)
      {  sqe = node->r_stat;
         node->r_stat = sqe->next;
         dmp_free_atom(iet->sqe_pool, sqe);
      }
      /* delete the column statuses change list */
      while (node->c_stat != NULL)
      {  sqe = node->c_stat;
         node->c_stat = sqe->next;
         dmp_free_atom(iet->sqe_pool, sqe);
      }
      return;
}

/*----------------------------------------------------------------------
-- iet_freeze_node - freeze current subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_freeze_node(IET *iet);
--
-- *Description*
--
-- The routine iet_freeze_node freezes the current subproblem. */

void iet_freeze_node(IET *iet)
{     IETNPD *node;
      IETBQE *bqe;
      IETCQE *cqe;
      IETAQE *aqe;
      IETAIJ *aij;
      IETSQE *sqe;
      IETROW *row;
      IETCOL *col;
      int i, j;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_freeze_node: current subproblem does not exist");
      /* build change lists for rows */
      insist(node->r_bnds == NULL);
      insist(node->r_mat == NULL);
      insist(node->r_stat == NULL);
      for (i = 1; i <= iet->m; i++)
      {  /* obtain pointer to local descriptor of i-th row */
         row = iet->row[i];
         /* save type and bounds of i-th row, if changed */
         if (!(row->type == row->old_type && row->lb == row->old_lb &&
               row->ub == row->old_ub))
         {  /* create type/bounds change entry */
            bqe = dmp_get_atom(iet->bqe_pool);
            bqe->u.row = row->glob;
            bqe->type = row->type;
            bqe->lb = row->lb;
            bqe->ub = row->ub;
            bqe->next = node->r_bnds;
            node->r_bnds = bqe;
         }
         /* save i-th row of the constraint matrix, if replaced */
         if (row->set_by == node)
         {  /* create constraint matrix change entry */
            aqe = dmp_get_atom(iet->aqe_pool);
            aqe->u.row = row->glob;
            aqe->ptr = NULL;
            aqe->next = node->r_mat;
            node->r_mat = aqe;
            /* attach all constraint coefficients in i-th row to the
               matrix change entry */
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
            {  aij->link = aqe->ptr;
               aqe->ptr = aij;
            }
         }
         /* save status of i-th row, if changed */
         if (row->stat != row->old_stat)
         {  /* create status change entry */
            sqe = dmp_get_atom(iet->sqe_pool);
            sqe->u.row = row->glob;
            sqe->stat = row->stat;
            sqe->next = node->r_stat;
            node->r_stat = sqe;
         }
      }
      /* build change lists for columns */
      insist(node->c_bnds == NULL);
      insist(node->c_obj == NULL);
      insist(node->c_mat == NULL);
      insist(node->c_stat == NULL);
      for (j = 1; j <= iet->n; j++)
      {  /* obtain pointer to local descriptor of j-th column */
         col = iet->col[j];
         /* save type and bounds of j-th column, if changed */
         if (!(col->type == col->old_type && col->lb == col->old_lb &&
               col->ub == col->old_ub))
         {  /* create type/bounds change entry */
            bqe = dmp_get_atom(iet->bqe_pool);
            bqe->u.col = col->glob;
            bqe->type = col->type;
            bqe->lb = col->lb;
            bqe->ub = col->ub;
            bqe->next = node->c_bnds;
            node->c_bnds = bqe;
         }
         /* save objective coefficient of j-th column, if changed */
         if (col->coef != col->old_coef)
         {  /* create objective coefficient change entry */
            cqe = dmp_get_atom(iet->cqe_pool);
            cqe->col = col->glob;
            cqe->coef = col->coef;
            cqe->next = node->c_obj;
            node->c_obj = cqe;
         }
         /* save j-th column of the constraint matrix, if replaced */
         if (col->set_by == node)
         {  /* create constraint matrix change entry */
            aqe = dmp_get_atom(iet->aqe_pool);
            aqe->u.col = col->glob;
            aqe->ptr = NULL;
            aqe->next = node->c_mat;
            node->c_mat = aqe;
            /* attach all constraint coefficients in j-th column to the
               matrix change entry except those coefficients which have
               been already attached to entries for rows replaced */
            for (aij = col->ptr; aij != NULL; aij = aij->c_next)
            {  if (iet->row[aij->row->i]->set_by != node)
               {  aij->link = aqe->ptr;
                  aqe->ptr = aij;
               }
            }
         }
         /* save status of j-th column, if changed */
         if (col->stat != col->old_stat)
         {  /* create status change entry */
            sqe = dmp_get_atom(iet->sqe_pool);
            sqe->u.col = col->glob;
            sqe->stat = col->stat;
            sqe->next = node->c_stat;
            node->c_stat = sqe;
         }
      }
      /* save constant term of the objective function, if changed */
      if (iet->c0 != iet->old_c0)
      {  /* create objective coefficient change entry */
         cqe = dmp_get_atom(iet->cqe_pool);
         cqe->col = NULL;
         cqe->coef = iet->c0;
         cqe->next = node->c_obj;
         node->c_obj = cqe;
      }
      /* delete row local descriptors */
      for (i = 1; i <= iet->m; i++)
      {  row = iet->row[i];
         insist(row->glob->i == i);
         row->glob->i = 0;
         dmp_free_atom(iet->row_pool, row);
      }
      /* delete column local descriptors */
      for (j = 1; j <= iet->n; j++)
      {  col = iet->col[j];
         insist(col->glob->j == j);
         col->glob->j = 0;
         dmp_free_atom(iet->col_pool, col);
      }
      /* the current subproblem has been frozen */
      iet->curr = NULL;
      iet->m = 0;
      iet->n = 0;
      iet->nz = 0;
      return;
}

/*----------------------------------------------------------------------
-- iet_clone_node - clone specified subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_clone_node(IET *iet, int p, int nnn);
--
-- *Description*
--
-- The routine iet_clone_node clones the specified subproblem, whose
-- reference number is p, creating its nnn exact copies. Note that the
-- specified subproblem must be active and must be in the frozen state
-- (i.e. it must not be the current subproblem).
--
-- Each clone, an exact copy of the specified subproblem, becomes a new
-- active subproblem added to the end of the active list. After cloning
-- the specified subproblem becomes inactive. */

void iet_clone_node(IET *iet, int p, int nnn)
{     IETNPD *node, *orig;
      /* obtain pointer to the subproblem to be cloned */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_clone_node: p = %d; invalid subproblem reference nu"
            "mber", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* the specified subproblem must be active */
      if (node->count != 0)
         fault("iet_clone_node: p = %d; cloning inactive subproblem not"
            " allowed", p);
      /* and must be in the frozen state */
      if (iet->curr == node)
         fault("iet_clone_node: p = %d; cloning current subproblem not "
            "allowed", p);
      /* remove the specified subproblem from the active list, because
         it becomes inactive */
      if (node->prev == NULL)
         iet->head = node->next;
      else
         node->prev->next = node->next;
      if (node->next == NULL)
         iet->tail = node->prev;
      else
         node->next->prev = node->prev;
      node->prev = node->next = NULL;
      iet->a_cnt--;
      /* set the child count of the specified subproblem, which is the
         number of clones to be created */
      if (nnn < 1)
         fault("iet_clone_node: nnn = %d; invalid number of clone subpr"
            "oblems", nnn);
      node->count = nnn;
      /* save pointer to the specified subproblem */
      orig = node;
      /* create clone subproblems */
      while (nnn-- > 0)
      {  /* if no free slots are available, increase the room */
         if (iet->avail == 0)
         {  int nslots = iet->nslots;
            IETNPS *save = iet->slot;
            iet->nslots = nslots + nslots;
            insist(iet->nslots > nslots);
            iet->slot = ucalloc(1+iet->nslots, sizeof(IETNPS));
            memcpy(&iet->slot[1], &save[1], nslots * sizeof(IETNPS));
            /* push more free slots into the stack */
            for (p = iet->nslots; p > nslots; p--)
            {  iet->slot[p].node = NULL;
               iet->slot[p].next = iet->avail;
               iet->avail = p;
            }
            ufree(save);
         }
         /* pull a free slot from the stack */
         p = iet->avail;
         iet->avail = iet->slot[p].next;
         insist(iet->slot[p].node == NULL);
         iet->slot[p].next = 0;
         /* create descriptor for new subproblem */
         iet->slot[p].node = node = dmp_get_atom(iet->npd_pool);
         node->p = p;
         node->up = orig;
         node->level = orig->level + 1;
         node->count = 0;
         node->r_add = NULL;
         node->c_add = NULL;
         node->r_del = NULL;
         node->c_del = NULL;
         node->r_bnds = NULL;
         node->c_bnds = NULL;
         node->c_obj = NULL;
         node->r_mat = NULL;
         node->c_mat = NULL;
         node->r_stat = NULL;
         node->c_stat = NULL;
         node->link = NULL;
         node->temp = NULL;
         node->prev = iet->tail;
         node->next = NULL;
         /* add the new subproblem to the end of the active list */
         if (iet->head == NULL)
            iet->head = node;
         else
            iet->tail->next = node;
         iet->tail = node;
         iet->a_cnt++;
         iet->n_cnt++;
         iet->t_cnt++;
      }
      return;
}

/*----------------------------------------------------------------------
-- iet_set_node_link - set link to subproblem extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_node_link(IET *iet, int p, void *link);
--
-- *Description*
--
-- The routine iet_set_node_link sets the link to a global extension of
-- the subproblem whose reference number is p. */

void iet_set_node_link(IET *iet, int p, void *link)
{     IETNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_set_node_link: p = %d; invalid subproblem reference"
            " number", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* set the link to subproblem extension */
      node->link = link;
      return;
}

/*----------------------------------------------------------------------
-- iet_delete_node - delete specified subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_delete_node(IET *iet, int p);
--
-- *Description*
--
-- The routine iet_delete_node deletes the specified subproblem, whose
-- reference number is p. The subproblem must be active and must be in
-- the frozen state (i.e. it must not be the current subproblem).
--
-- Note that deletion is performed recursively, i.e. if a subproblem to
-- be deleted is the only child of its parent, the parent subproblem is
-- also deleted, etc. */

void iet_delete_node(IET *iet, int p)
{     IETNPD *node, *temp;
      IETRGD *rgd;
      IETCGD *cgd;
      IETDQE *dqe;
      IETBQE *bqe;
      IETCQE *cqe;
      IETAQE *aqe;
      IETAIJ *aij;
      IETSQE *sqe;
      /* obtain pointer to the subproblem to be deleted */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_delete_node: p = %d; invalid subproblem reference n"
            "umber", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* the specified subproblem must be active */
      if (node->count != 0)
         fault("iet_delete_node: p = %d; deleting inactive subproblem n"
            "ot allowed", p);
      /* and must be in the frozen state */
      if (iet->curr == node)
         fault("iet_delete_node: p = %d; deleting current subproblem no"
            "t allowed", p);
      /* remove the specified subproblem from the active list, because
         it is gone from the tree */
      if (node->prev == NULL)
         iet->head = node->next;
      else
         node->prev->next = node->next;
      if (node->next == NULL)
         iet->tail = node->prev;
      else
         node->next->prev = node->prev;
      node->prev = node->next = NULL;
      iet->a_cnt--;
loop: /* recursive deletion starts here */
      /* delete global descriptors of own rows */
      while (node->r_add != NULL)
      {  rgd = node->r_add;
         node->r_add = rgd->next;
         if (iet->hook != NULL)
            iet->hook(iet->info, IET_RD, rgd->name == NULL ? NULL :
               get_str(iet->str_buf, rgd->name), rgd->link);
         if (rgd->name != NULL) delete_str(rgd->name);
         dmp_free_atom(iet->rgd_pool, rgd);
      }
      /* delete global descriptors of own columns */
      while (node->c_add != NULL)
      {  cgd = node->c_add;
         node->c_add = cgd->next;
         if (iet->hook != NULL)
            iet->hook(iet->info, IET_CD, cgd->name == NULL ? NULL :
               get_str(iet->str_buf, cgd->name), cgd->link);
         if (cgd->name != NULL) delete_str(cgd->name);
         dmp_free_atom(iet->cgd_pool, cgd);
      }
      /* delete the rows deletion list */
      while (node->r_del != NULL)
      {  dqe = node->r_del;
         node->r_del = dqe->next;
         dmp_free_atom(iet->dqe_pool, dqe);
      }
      /* delete the columns deletion list */
      while (node->c_del != NULL)
      {  dqe = node->c_del;
         node->c_del = dqe->next;
         dmp_free_atom(iet->dqe_pool, dqe);
      }
      /* delete the row types/bounds change list */
      while (node->r_bnds != NULL)
      {  bqe = node->r_bnds;
         node->r_bnds = bqe->next;
         dmp_free_atom(iet->bqe_pool, bqe);
      }
      /* delete the column types/bounds change list */
      while (node->c_bnds != NULL)
      {  bqe = node->c_bnds;
         node->c_bnds = bqe->next;
         dmp_free_atom(iet->bqe_pool, bqe);
      }
      /* delete the column objective coefficients change list */
      while (node->c_obj != NULL)
      {  cqe = node->c_obj;
         node->c_obj = cqe->next;
         dmp_free_atom(iet->cqe_pool, cqe);
      }
      /* delete the change list for rows of the constraint matrix */
      while (node->r_mat != NULL)
      {  aqe = node->r_mat;
         node->r_mat = aqe->next;
         while (aqe->ptr != NULL)
         {  aij = aqe->ptr;
            aqe->ptr = aij->link;
            dmp_free_atom(iet->aij_pool, aij);
         }
         dmp_free_atom(iet->aqe_pool, aqe);
      }
      /* delete the change list for columns of the constraint matrix */
      while (node->c_mat != NULL)
      {  aqe = node->c_mat;
         node->c_mat = aqe->next;
         while (aqe->ptr != NULL)
         {  aij = aqe->ptr;
            aqe->ptr = aij->link;
            dmp_free_atom(iet->aij_pool, aij);
         }
         dmp_free_atom(iet->aqe_pool, aqe);
      }
      /* delete the row statuses change list */
      while (node->r_stat != NULL)
      {  sqe = node->r_stat;
         node->r_stat = sqe->next;
         dmp_free_atom(iet->sqe_pool, sqe);
      }
      /* delete the column statuses change list */
      while (node->c_stat != NULL)
      {  sqe = node->c_stat;
         node->c_stat = sqe->next;
         dmp_free_atom(iet->sqe_pool, sqe);
      }
      /* free the corresponding node slot */
      p = node->p;
      insist(iet->slot[p].node == node);
      iet->slot[p].node = NULL;
      iet->slot[p].next = iet->avail;
      iet->avail = p;
      /* save pointer to the parent subproblem */
      temp = node->up;
      /* delete the subproblem descriptor */
      if (iet->hook != NULL)
         iet->hook(iet->info, IET_ND, NULL, node->link);
      dmp_free_atom(iet->npd_pool, node);
      iet->n_cnt--;
      /* take pointer to the parent subproblem */
      node = temp;
      if (node != NULL)
      {  /* the parent subproblem exists; decrease the number of its
            child subproblems */
         insist(node->count > 0);
         node->count--;
         /* if now the parent subproblem has no childs, it also must be
            deleted */
         if (node->count == 0) goto loop;
      }
      return;
}

/*----------------------------------------------------------------------
-- iet_delete_tree - delete implicit enumeration tree.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_delete_tree(IET *iet);
--
-- *Description*
--
-- The routine iet_delete_tree deletes the implicit enumeration tree,
-- which the parameter iet points to, freeing all the memory allocated
-- to this object. Note that if the current subproblem exists, it must
-- be frozen before deleting the tree. */

void iet_delete_tree(IET *iet)
{     /* the current subproblem must not exist */
      if (iet->curr != NULL)
         fault("iet_delete_tree: current subproblem still exists");
      /* delete all subproblems which are still active */
      while (iet->tail != NULL) iet_delete_node(iet, iet->tail->p);
      /* now the tree must be empty */
      insist(iet->a_cnt == 0);
      insist(iet->n_cnt == 0);
      /* and all atoms must return to their memory pools */
      insist(iet->npd_pool->count == 0);
      insist(iet->rgd_pool->count == 0);
      insist(iet->cgd_pool->count == 0);
      insist(iet->dqe_pool->count == 0);
      insist(iet->bqe_pool->count == 0);
      insist(iet->cqe_pool->count == 0);
      insist(iet->aqe_pool->count == 0);
      insist(iet->aij_pool->count == 0);
      insist(iet->sqe_pool->count == 0);
      insist(iet->row_pool->count == 0);
      insist(iet->col_pool->count == 0);
      insist(iet->str_pool->count == 0);
      /* free all the memory */
      dmp_delete_pool(iet->npd_pool);
      dmp_delete_pool(iet->rgd_pool);
      dmp_delete_pool(iet->cgd_pool);
      dmp_delete_pool(iet->dqe_pool);
      dmp_delete_pool(iet->bqe_pool);
      dmp_delete_pool(iet->cqe_pool);
      dmp_delete_pool(iet->aqe_pool);
      dmp_delete_pool(iet->aij_pool);
      dmp_delete_pool(iet->sqe_pool);
      dmp_delete_pool(iet->row_pool);
      dmp_delete_pool(iet->col_pool);
      dmp_delete_pool(iet->str_pool);
      ufree(iet->str_buf);
      ufree(iet->slot);
      ufree(iet->row);
      ufree(iet->col);
      ufree(iet);
      return;
}

/**********************************************************************/
/* * *                  TREE EXPLORING ROUTINES                   * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- iet_get_tree_size - determine current size of the tree.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_get_tree_size(IET *iet, int *a_cnt, int *n_cnt, int *t_cnt);
--
-- *Description*
--
-- The routine iet_get_tree_size stores the following three counts which
-- characterize the current size of the tree:
--
-- a_cnt is the current number of active nodes, i.e. the current size of
--       the active list;
--
-- n_cnt is the current number of all (active and inactive) nodes;
--
-- t_cnt is the total number of nodes including those which have been
--       already removed from the tree. This count is increased whenever
--       a new node appears in the tree and never decreased.
--
-- If some of the parameters a_cnt, n_cnt, t_cnt is a null pointer, the
-- corresponding count is not stored. */

void iet_get_tree_size(IET *iet, int *a_cnt, int *n_cnt, int *t_cnt)
{     if (a_cnt != NULL) *a_cnt = iet->a_cnt;
      if (n_cnt != NULL) *n_cnt = iet->n_cnt;
      if (t_cnt != NULL) *t_cnt = iet->t_cnt;
      return;
}

/*----------------------------------------------------------------------
-- iet_get_curr_node - determine current active subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_curr_node(IET *iet);
--
-- *Returns*
--
-- The routine iet_get_curr_node returns the reference number of the
-- current active subproblem. However, if the current subproblem does
-- not exist, the routine returns zero. */

int iet_get_curr_node(IET *iet)
{     IETNPD *node;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      /* return its reference number */
      return node == NULL ? 0 : node->p;
}

/*----------------------------------------------------------------------
-- iet_get_next_node - determine next active subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_next_node(IET *iet, int p);
--
-- *Returns*
--
-- If the parameter p is zero, the routine iet_get_next_node returns
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

int iet_get_next_node(IET *iet, int p)
{     IETNPD *node;
      if (p == 0)
      {  /* obtain pointer to the first active subproblem */
         node = iet->head;
      }
      else
      {  /* obtain pointer to the specified subproblem */
         if (!(1 <= p && p <= iet->nslots))
err:        fault("iet_get_next_node: p = %d; invalid subproblem refere"
               "nce number", p);
         node = iet->slot[p].node;
         if (node == NULL) goto err;
         /* the specified subproblem must be active */
         if (node->count != 0)
            fault("iet_get_next_node: p = %d; subproblem not in the act"
               "ive list", p);
         /* obtain pointer to the next active subproblem */
         node = node->next;
      }
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/*----------------------------------------------------------------------
-- iet_get_prev_node - determine previous active subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_prev_node(IET *iet, int p);
--
-- *Returns*
--
-- If the parameter p is zero, the routine iet_get_prev_node returns
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

int iet_get_prev_node(IET *iet, int p)
{     IETNPD *node;
      if (p == 0)
      {  /* obtain pointer to the last active subproblem */
         node = iet->tail;
      }
      else
      {  /* obtain pointer to the specified subproblem */
         if (!(1 <= p && p <= iet->nslots))
err:        fault("iet_get_prev_node: p = %d; invalid subproblem refere"
               "nce number", p);
         node = iet->slot[p].node;
         if (node == NULL) goto err;
         /* the specified subproblem must be active */
         if (node->count != 0)
            fault("iet_get_prev_node: p = %d; subproblem not in the act"
               "ive list", p);
         /* obtain pointer to the previous active subproblem */
         node = node->prev;
      }
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/*----------------------------------------------------------------------
-- iet_get_up_node - determine parent subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_up_node(IET *iet, int p);
--
-- *Returns*
--
-- The parameter p must specify the reference number of some (active or
-- inactive) subproblem, in which case the routine iet_get_up_node
-- returns the reference number of its parent subproblem. However, if
-- the specified subproblem is the root of the tree and therefore has no
-- parent, the routine returns zero. */

int iet_get_up_node(IET *iet, int p)
{     IETNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_get_up_node: p = %d; invalid subproblem reference n"
            "umber", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* obtain pointer to the parent subproblem */
      node = node->up;
      /* return the reference number */
      return node == NULL ? 0 : node->p;
}

/*----------------------------------------------------------------------
-- iet_get_node_lev - determine subproblem level.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_node_lev(IET *iet, int p);
--
-- *Returns*
--
-- The routine iet_get_node_lev returns the level which the subproblem
-- whose reference number is p has in the tree. (The root subproblem
-- has the level 0, and the level of any other subproblem is the level
-- of its parent plus one.) However, if the parameter p is not a valid
-- subproblem reference number, the routine returns negative value. */

int iet_get_node_lev(IET *iet, int p)
{     IETNPD *node;
      /* obtain pointer to the specified subproblem */
      node = (1 <= p && p <= iet->nslots ? iet->slot[p].node : NULL);
      /* return the subproblem level */
      return node == NULL ? -1 : node->level;
}

/*----------------------------------------------------------------------
-- iet_get_node_cnt - determine number of child subproblems.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_node_cnt(IET *iet, int p);
--
-- *Returns*
--
-- The routine iet_get_node_cnt returns the number of child subproblems
-- for the subproblem whose reference number is p.
--
-- Zero means the subproblem p is active and therefore has no childs.
--
-- Positive value means the subproblem p is inactive where the value is
-- the number of its childs.
--
-- Negative value means the reference number p is invalid. */

int iet_get_node_cnt(IET *iet, int p)
{     IETNPD *node;
      /* obtain pointer to the specified subproblem */
      node = (1 <= p && p <= iet->nslots ? iet->slot[p].node : NULL);
      /* return the number of child subproblems */
      return node == NULL ? -1 : node->count;
}

/*----------------------------------------------------------------------
-- iet_get_node_link - obtain link to subproblem extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void *iet_get_node_link(IET *iet, int p);
--
-- *Returns*
--
-- The routine iet_get_node_link returns the link to a global extension
-- of the subproblem whose reference number is p. */

void *iet_get_node_link(IET *iet, int p)
{     IETNPD *node;
      /* obtain pointer to the specified subproblem */
      if (!(1 <= p && p <= iet->nslots))
err:     fault("iet_get_node_link: p = %d; invalid subproblem reference"
            " number", p);
      node = iet->slot[p].node;
      if (node == NULL) goto err;
      /* return the link to subproblem extension */
      return node->link;
}

/*----------------------------------------------------------------------
-- iet_pseudo_root - find pseudo-root of the tree.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_pseudo_root(IET *iet);
--
-- *Description*
--
-- The routine iet_pseudo_root finds so-called pseudo-root of the tree,
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
-- *Returns*
--
-- The routine iet_pseudo_root returns the subproblem reference number
-- which corresponds to the pseudo-root. However, if the tree is empty,
-- zero is returned. */

int iet_pseudo_root(IET *iet)
{     IETNPD *root, *node;
      /* obtain pointer to the root node of the entire tree */
      root = iet->slot[1].node;
      /* if the tree is empty, there is no pseudo-root node */
      if (root == NULL) goto done;
      /* obtain pointer to any active node */
      node = iet->head;
      insist(node != NULL);
      /* build the path from the root to the active node */
      node->temp = NULL;
      for (node = node; node != NULL; node = node->up)
      {  if (node->up == NULL)
            insist(node == root);
         else
            node->up->temp = node;
      }
      /* walk from the root to the active node and find the pseudo-root
         of the tree */
      for (root = root; root != NULL; root = root->temp)
         if (root->count != 1) break;
      insist(root != NULL);
done: /* return the reference number of the pseudo-root found */
      return root == NULL ? 0 : root->p;
}

/**********************************************************************/
/* * *               SUBPROBLEM MODIFYING ROUTINES                * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- iet_add_rows - add new rows to current subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_add_rows(IET *iet, int nrs);
--
-- *Description*
--
-- The routine iet_add_rows adds nrs rows (constraints) to the current
-- subproblem. New rows are always added to the end of the row list, so
-- the ordinal numbers assigned to existing rows are not changed.
--
-- Each new row is initially free (unbounded) and has empty list of the
-- constraint coefficients. */

void iet_add_rows(IET *iet, int nrs)
{     IETNPD *node;
      IETRGD *rgd;
      IETROW *row;
      int m_new, i;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_add_rows: current subproblem does not exist");
      /* determine new number of rows */
      if (nrs < 1)
         fault("iet_add_rows: nrs = %d; invalid parameter", nrs);
      m_new = iet->m + nrs;
      insist(m_new > 0);
      /* increase the room, if necessary */
      if (iet->m_max < m_new)
      {  IETROW **save = iet->row;
         while (iet->m_max < m_new)
         {  iet->m_max += iet->m_max;
            insist(iet->m_max > 0);
         }
         iet->row = ucalloc(1+iet->m_max, sizeof(IETROW *));
         memcpy(&iet->row[1], &save[1], iet->m * sizeof(IETROW *));
         ufree(save);
      }
      /* add new rows to the end of the row list */
      for (i = iet->m+1; i <= m_new; i++)
      {  /* create global descriptor */
         rgd = dmp_get_atom(iet->rgd_pool);
         rgd->host = node;
         rgd->name = NULL;
         rgd->i = i;
         rgd->link = NULL;
         rgd->temp = NULL;
         rgd->next = NULL;
         /* and add it to the end of the linked list of own rows of the
            current subproblem */
         if (node->r_add == NULL)
         {  /* it is the first own row */
            node->r_add = rgd;
         }
         else
         {  /* it is not the first own row */
            insist(i > 1);
            row = iet->row[i-1];
            insist(row->glob->host == node);
            insist(row->glob->next == NULL);
            row->glob->next = rgd;
         }
         /* create local descriptor */
         iet->row[i] = row = dmp_get_atom(iet->row_pool);
         row->glob = rgd;
         row->type = IET_FR;
         row->lb = 0.0;
         row->ub = 0.0;
         row->set_by = node;
         row->ptr = NULL;
         row->stat = IET_BS;
         row->old_type = row->type;
         row->old_lb = row->lb;
         row->old_ub = row->ub;
         row->old_stat = row->stat;
         row->link = NULL;
      }
      /* set new number of rows */
      iet->m = m_new;
      return;
}

/*----------------------------------------------------------------------
-- iet_add_cols - add new columns to current subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_add_cols(IET *iet, int ncs);
--
-- *Description*
--
-- The routine iet_add_cols adds ncs columns (structural variables) to
-- the current subproblem. New columns are always added to the end of
-- the column list, so the ordinal numbers assigned to existng columns
-- are not changed.
--
-- Each new column is initially fixed at zero and has empty list of the
-- constraint coefficients. */

void iet_add_cols(IET *iet, int ncs)
{     IETNPD *node;
      IETCGD *cgd;
      IETCOL *col;
      int n_new, j;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_add_cols: current subproblem does not exist");
      /* determine new number of columns */
      if (ncs < 1)
         fault("iet_add_cols: ncs = %d; invalid parameter", ncs);
      n_new = iet->n + ncs;
      insist(n_new > 0);
      /* increase the room, if necessary */
      if (iet->n_max < n_new)
      {  IETCOL **save = iet->col;
         while (iet->n_max < n_new)
         {  iet->n_max += iet->n_max;
            insist(iet->n_max > 0);
         }
         iet->col = ucalloc(1+iet->n_max, sizeof(IETCOL *));
         memcpy(&iet->col[1], &save[1], iet->n * sizeof(IETCOL *));
         ufree(save);
      }
      /* add new columns to the end of the column list */
      for (j = iet->n+1; j <= n_new; j++)
      {  /* create global descriptor */
         cgd = dmp_get_atom(iet->cgd_pool);
         cgd->host = node;
         cgd->name = NULL;
         cgd->j = j;
         cgd->link = NULL;
         cgd->temp = NULL;
         cgd->next = NULL;
         /* and add it to the end of the linked list of own columns of
            the current subproblem */
         if (node->c_add == NULL)
         {  /* it is the first own column */
            node->c_add = cgd;
         }
         else
         {  /* it is not the first own column */
            insist(j > 1);
            col = iet->col[j-1];
            insist(col->glob->host == node);
            insist(col->glob->next == NULL);
            col->glob->next = cgd;
         }
         /* create local descriptor */
         iet->col[j] = col = dmp_get_atom(iet->col_pool);
         col->glob = cgd;
         col->type = IET_FX;
         col->lb = 0.0;
         col->ub = 0.0;
         col->coef = 0.0;
         col->set_by = node;
         col->ptr = NULL;
         col->stat = IET_NS;
         col->old_type = col->type;
         col->old_lb = col->lb;
         col->old_ub = col->ub;
         col->old_coef = col->coef;
         col->old_stat = col->stat;
         col->link = NULL;
      }
      /* set new number of columns */
      iet->n = n_new;
      return;
}

/*----------------------------------------------------------------------
-- iet_check_name - check correctness of symbolic name.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_check_name(IET *iet, char *name);
--
-- *Description*
--
-- The routine iet_check_name checks if given symbolic name is correct
-- (a name is correct if it contains 1 up to 255 graphic characters).
--
-- *Returns*
--
-- If the symbolic name is correct, the routine returns zero. Otherwise
-- non-zero is returned. */

int iet_check_name(IET *iet, char *name)
{     int t;
      insist(iet == iet);
      if (name[0] == '\0') return 1; /* empty */
      for (t = 0; name[t] != '\0'; t++)
      {  if (t == 255) return 1; /* too long */
         if (!isgraph((unsigned char)name[t])) return 1; /* ugly */
      }
      return 0; /* nice */
}

/*----------------------------------------------------------------------
-- iet_set_row_name - assign symbolic name to row.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_row_name(IET *iet, int i, char *name);
--
-- *Description*
--
-- The routine iet_set_row_name assigns a given symbolic name to i-th
-- row of the current subproblem.
--
-- If the parameter name is NULL, the routine just erases existing name
-- of the i-th row.
--
-- NOTE: Changing the row name has the global effect. */

void iet_set_row_name(IET *iet, int i, char *name)
{     IETRGD *rgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_row_name: current subproblem does not exist");
      /* obtain pointer to global descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_row_name: i = %d; row number out of range", i);
      rgd = iet->row[i]->glob;
      /* assign given symbolic name to i-th row */
      if (name == NULL)
      {  if (rgd->name != NULL)
         {  delete_str(rgd->name);
            rgd->name = NULL;
         }
      }
      else
      {  if (iet_check_name(iet, name))
            fault("iet_set_row_name: i = %d; invalid name", i);
         if (rgd->name == NULL)
            rgd->name = create_str(iet->str_pool);
         set_str(rgd->name, name);
      }
      return;
}

/*----------------------------------------------------------------------
-- iet_set_col_name - assign symbolic name to column.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_col_name(IET *iet, int j, char *name);
--
-- *Description*
--
-- The routine iet_set_col_name assigns a given symbolic name to j-th
-- column of the current subproblem.
--
-- If the parameter name is NULL, the routine just erases existing name
-- of the j-th column.
--
-- NOTE: Changing the column name has the global effect. */

void iet_set_col_name(IET *iet, int j, char *name)
{     IETCGD *cgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_col_name: current subproblem does not exist");
      /* obtain pointer to global descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_col_name: j = %d; column number out of range",
            j);
      cgd = iet->col[j]->glob;
      /* assign given symbolic name to j-th column */
      if (name == NULL)
      {  if (cgd->name != NULL)
         {  delete_str(cgd->name);
            cgd->name = NULL;
         }
      }
      else
      {  if (iet_check_name(iet, name))
            fault("iet_set_col_name: j = %d; invalid name", j);
         if (cgd->name == NULL)
            cgd->name = create_str(iet->str_pool);
         set_str(cgd->name, name);
      }
      return;
}

/*----------------------------------------------------------------------
-- iet_set_row_link - set link to row global extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_row_link(IET *iet, int i, void *link);
--
-- *Description*
--
-- The routine iet_set_row_link sets the link to a global extension of
-- i-th row of the current subproblem. */

void iet_set_row_link(IET *iet, int i, void *link)
{     IETRGD *rgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_row_link: current subproblem does not exist");
      /* obtain pointer to global descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_row_link: i = %d; row number out of range", i);
      rgd = iet->row[i]->glob;
      /* set the link to global extension */
      rgd->link = link;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_col_link - set link to column global extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_col_link(IET *iet, int j, void *link);
--
-- *Description*
--
-- The routine iet_set_col_link sets the link to a global extension of
-- j-th column of the current subproblem. */

void iet_set_col_link(IET *iet, int j, void *link)
{     IETCGD *cgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_col_link: current subproblem does not exist");
      /* obtain pointer to global descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_col_link: j = %d; column number out of range",
            j);
      cgd = iet->col[j]->glob;
      /* set the link to global extension */
      cgd->link = link;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_row_bnds - set row type and bounds.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_row_bnds(IET *iet, int i, int type, double lb,
--    double ub);
--
-- *Description*
--
-- The routine iet_set_row_bnds sets the type and bounds of i-th row of
-- the current subproblem.
--
-- The parameters type, lb, and ub specify the type, lower bound, and
-- upper bound, respectively, as shown below:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IET_FR   -inf <  x <  +inf   free variable
--    IET_LO     lb <= x <  +inf   lower bound
--    IET_UP   -inf <  x <=  ub    upper bound
--    IET_DB     lb <= x <=  ub    double bound
--    IET_FX           x  =  lb    fixed variable
--
-- where x is auxiliary variable associated with i-th row.
--
-- If the row has no lower bound, the parameter lb is ignored. If the
-- row has no upper bound, the parameter ub is ignored. If the row is
-- of fixed type, the parameter lb is used, and the parameter ub is
-- ignored. */

void iet_set_row_bnds(IET *iet, int i, int type, double lb, double ub)
{     IETROW *row;
      int stat;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_row_bnds: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_row_bnds: i = %d; row number out of range", i);
      row = iet->row[i];
      /* adjust bounds and determine non-basic status */
      switch (type)
      {  case IET_FR:
            lb = ub = 0.0; stat = IET_NF; break;
         case IET_LO:
            ub = 0.0; stat = IET_NL; break;
         case IET_UP:
            lb = 0.0; stat = IET_NU; break;
         case IET_DB:
            if (lb >= ub)
               fault("iet_set_row_bnds: i = %d; lb = %.*g; ub = %.*g; i"
                  "nvalid row bounds", i, DBL_DIG, lb, DBL_DIG, ub);
            stat = row->stat;
            if (!(stat == IET_NL || stat == IET_NU))
               stat = (fabs(lb) <= fabs(ub) ? IET_NL : IET_NU);
            break;
         case IET_FX:
            ub = lb; stat = IET_NS; break;
         default:
            fault("iet_set_row_bnds: i = %d; type = %d; invalid row typ"
               "e", i, type);
      }
      /* set new type and bounds */
      row->type = type;
      row->lb = lb;
      row->ub = ub;
      /* if the row is non-basic, adjust its status */
      if (row->stat != IET_BS) row->stat = stat;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_col_bnds - set column type and bounds.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_col_bnds(IET *iet, int j, int type, double lb,
--    double ub);
--
-- *Description*
--
-- The routine iet_set_col_bnds sets the type and bounds of j-th column
-- of the current subproblem.
--
-- The parameters type, lb, and ub specify the type, lower bound, and
-- upper bound, respectively, as shown below:
--
--     Type          Bounds            Note
--    -------------------------------------------
--    IET_FR   -inf <  x <  +inf   free variable
--    IET_LO     lb <= x <  +inf   lower bound
--    IET_UP   -inf <  x <=  ub    upper bound
--    IET_DB     lb <= x <=  ub    double bound
--    IET_FX           x  =  lb    fixed variable
--
-- where x is structural variable associated with j-th column.
--
-- If the column has no lower bound, the parameter lb is ignored. If
-- the column has no upper bound, the parameter ub is ignored. If the
-- column is of fixed type, the parameter lb is used, and the parameter
-- ub is ignored. */

void iet_set_col_bnds(IET *iet, int j, int type, double lb, double ub)
{     IETCOL *col;
      int stat;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_col_bnds: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_col_bnds: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* adjust bounds and determine non-basic status */
      switch (type)
      {  case IET_FR:
            lb = ub = 0.0; stat = IET_NF; break;
         case IET_LO:
            ub = 0.0; stat = IET_NL; break;
         case IET_UP:
            lb = 0.0; stat = IET_NU; break;
         case IET_DB:
            if (lb >= ub)
               fault("iet_set_col_bnds: j = %d; lb = %.*g; ub = %.*g; i"
                  "nvalid column bounds", j, DBL_DIG, lb, DBL_DIG, ub);
            stat = col->stat;
            if (!(stat == IET_NL || stat == IET_NU))
               stat = (fabs(lb) <= fabs(ub) ? IET_NL : IET_NU);
            break;
         case IET_FX:
            ub = lb; stat = IET_NS; break;
         default:
            fault("iet_set_col_bnds: j = %d; type = %d; invalid column "
               "type", j, type);
      }
      /* set new type and bounds */
      col->type = type;
      col->lb = lb;
      col->ub = ub;
      /* if the column is non-basic, adjust its status */
      if (col->stat != IET_BS) col->stat = stat;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_obj_coef - set objective coefficient or constant term.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_obj_coef(IET *iet, int j, double coef);
--
-- *Description*
--
-- The routine iet_set_obj_coef sets the objective coefficient at j-th
-- column of the current subproblem.
--
-- If the parameter j is 0, the routine sets the constant term (shift)
-- of the objective function of the current subproblem. */

void iet_set_obj_coef(IET *iet, int j, double coef)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_obj_coef: current subproblem does not exist");
      if (j == 0)
      {  /* set the constant term */
         iet->c0 = coef;
         goto done;
      }
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_obj_coef: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* set the objective coefficient */
      col->coef = coef;
done: return;
}

/*----------------------------------------------------------------------
-- iet_set_mat_row - replace row of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_mat_row(IET *iet, int i, int len, int ind[],
--    double val[]);
--
-- *Description*
--
-- The routine iet_set_mat_row replaces the contents of i-th row of the
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

void iet_set_mat_row(IET *iet, int i, int len, int ind[], double val[])
{     IETNPD *node;
      IETROW *row;
      IETCOL *col;
      IETAIJ *aij;
      int j, t;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_set_mat_row: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_mat_row: i = %d; row number out of range", i);
      row = iet->row[i];
      /* remove all existing elements from i-th row */
      while (row->ptr != NULL)
      {  /* take next element in i-th row */
         aij = row->ptr;
         /* remove the element from the row list */
         row->ptr = aij->r_next;
         /* determine number j of corresponding column */
         j = aij->col->j;
         /* obtain pointer to local descriptor of j-th column */
         insist(1 <= j && j <= iet->n);
         col = iet->col[j];
         /* remove the element from the column list */
         if (aij->c_prev == NULL)
            col->ptr = aij->c_next;
         else
            aij->c_prev->c_next = aij->c_next;
         if (aij->c_next == NULL)
            ;
         else
            aij->c_next->c_prev = aij->c_prev;
         /* check if the element is created by the current subproblem;
            if so, return it to the memory pool */
         if (row->set_by == node || col->set_by == node)
            dmp_free_atom(iet->aij_pool, aij);
         /* decrease the non-zero count */
         iet->nz--;
      }
      /* create new contents of i-th row */
      if (!(0 <= len && len <= iet->n))
         fault("iet_set_mat_row: i = %d; len = %d; invalid row length",
            i, len);
      for (t = 1; t <= len; t++)
      {  /* take number j of corresponding column */
         j = ind[t];
         /* obtain pointer to local descriptor of j-th column */
         if (!(1 <= j && j <= iet->n))
            fault("iet_set_mat_row: i = %d; ind[%d] = %d; column index "
               "out of range", i, t, j);
         col = iet->col[j];
         /* if there is element with the same column index, it can only
            be found in the beginning of the list of j-th column */
         if (col->ptr != NULL && col->ptr->row->i == i)
            fault("iet_set_mat_row: i = %d; ind[%d] = %d; duplicate col"
               "umn indices not allowed", i, t, j);
         /* create new element */
         aij = dmp_get_atom(iet->aij_pool);
         aij->row = row->glob;
         aij->col = col->glob;
         if (val[t] == 0.0)
            fault("iet_set_mat_row: i = %d; ind[%d] = %d; zero element "
               "not allowed", i, t, j);
         aij->val = val[t];
         aij->link = NULL;
         aij->r_prev = NULL;
         aij->r_next = row->ptr;
         aij->c_prev = NULL;
         aij->c_next = col->ptr;
         /* add the new element to the beginning of the lists of i-th
            row and j-th column */
         if (row->ptr != NULL) row->ptr->r_prev = aij;
         if (col->ptr != NULL) col->ptr->c_prev = aij;
         row->ptr = col->ptr = aij;
         /* increase the non-zero count */
         iet->nz++;
      }
      /* i-th row has been replaced by the current subproblem */
      row->set_by = node;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_mat_col - replace column of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_mat_col(IET *iet, int j, int len, int ind[],
--    double val[]);
--
-- *Description*
--
-- The routine iet_set_mat_col replaces the contents of j-th column of
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

void iet_set_mat_col(IET *iet, int j, int len, int ind[], double val[])
{     IETNPD *node;
      IETROW *row;
      IETCOL *col;
      IETAIJ *aij;
      int i, t;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_set_mat_col: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_mat_col: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* remove all existing elements from j-th column */
      while (col->ptr != NULL)
      {  /* take next element in j-th column */
         aij = col->ptr;
         /* remove the element from the column list */
         col->ptr = aij->c_next;
         /* determine number i of corresponding row */
         i = aij->row->i;
         /* obtain pointer to local descriptor of i-th row */
         insist(1 <= i && i <= iet->m);
         row = iet->row[i];
         /* remove the element from the row list */
         if (aij->r_prev == NULL)
            row->ptr = aij->r_next;
         else
            aij->r_prev->r_next = aij->r_next;
         if (aij->r_next == NULL)
            ;
         else
            aij->r_next->r_prev = aij->r_prev;
         /* check if the element is created by the current subproblem;
            if so, return it to the memory pool */
         if (row->set_by == node || col->set_by == node)
            dmp_free_atom(iet->aij_pool, aij);
         /* decrease the non-zero count */
         iet->nz--;
      }
      /* create new contents of j-th column */
      if (!(0 <= len && len <= iet->m))
         fault("iet_set_mat_col: j = %d; len = %d; invalid column lengt"
            "h", j, len);
      for (t = 1; t <= len; t++)
      {  /* take number i of corresponding row */
         i = ind[t];
         /* obtain pointer to local descriptor of i-th row */
         if (!(1 <= i && i <= iet->m))
            fault("iet_set_mat_col: j = %d; ind[%d] = %d; row index out"
               " of range", j, t, i);
         row = iet->row[i];
         /* if there is element with the same row index, it can only be
            found in the beginning of the list of i-th row */
         if (row->ptr != NULL && row->ptr->col->j == j)
            fault("iet_set_mat_col: j = %d; ind[%d] = %d; duplicate row"
               " indices now allowed", j, t, i);
         /* create new element */
         aij = dmp_get_atom(iet->aij_pool);
         aij->row = row->glob;
         aij->col = col->glob;
         if (val[t] == 0.0)
            fault("iet_set_mat_col: j = %d; ind[%d] = %d; zero element "
               "not allowed", j, t, i);
         aij->val = val[t];
         aij->link = NULL;
         aij->r_prev = NULL;
         aij->r_next = row->ptr;
         aij->c_prev = NULL;
         aij->c_next = col->ptr;
         /* add the new element to the beginning of the lists of i-th
            row and j-th column */
         if (row->ptr != NULL) row->ptr->r_prev = aij;
         if (col->ptr != NULL) col->ptr->c_prev = aij;
         row->ptr = col->ptr = aij;
         /* increase the non-zero count */
         iet->nz++;
      }
      /* j-th column has been replaced by the current subproblem */
      col->set_by = node;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_row_stat - set row status.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_row_stat(IET *iet, int i, int stat);
--
-- *Description*
--
-- The routine iet_set_row_stat sets (changes) the status of i-th row
-- of the current subproblem as specified by the parameter stat:
--
-- IET_BS   - make the row basic (make the constraint inactive);
-- IET_NL   - make the row non-basic (make the constraint active);
-- IET_NU   - make the row non-basic and set it to the upper bound; if
--            the row is not double-bounded, this status is equivalent
--            to IET_NL (only for this routine);
-- IET_NF   - the same as IET_NL (only for this routine);
-- IET_NS   - the same as IET_NL (only for this routine). */

void iet_set_row_stat(IET *iet, int i, int stat)
{     IETROW *row;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_row_stat: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_row_stat: i = %d; row number out of range", i);
      row = iet->row[i];
      /* adjust non-basic status depending on the row type */
      if (!(stat == IET_BS || stat == IET_NL || stat == IET_NU ||
            stat == IET_NF || stat == IET_NS))
         fault("iet_set_row_stat: i = %d; stat = %d; invalid row status"
            , i, stat);
      if (stat != IET_BS)
      {  switch (row->type)
         {  case IET_FR:
               stat = IET_NF; break;
            case IET_LO:
               stat = IET_NL; break;
            case IET_UP:
               stat = IET_NU; break;
            case IET_DB:
               if (!(stat == IET_NL || stat == IET_NU)) stat = IET_NL;
               break;
            case IET_FX:
               stat = IET_NS; break;
            default:
               insist(row != row);
         }
      }
      /* set row status */
      row->stat = stat;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_col_stat - set column status.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_col_stat(IET *iet, int j, int stat);
--
-- *Description*
--
-- The routine iet_set_col_stat sets (changes) the status of j-th column
-- of the current subproblem as specified by the parameter stat:
--
-- IET_BS   - make the column basic;
-- IET_NL   - make the column non-basic;
-- IET_NU   - make the column non-basic and set it to the upper bound;
--            if the column is not of double-bounded type, this status
--            is the same as IET_NL (only for this routine);
-- IET_NF   - the same as IET_NL (only for this routine);
-- IET_NS   - the same as IET_NL (only for this routine). */

void iet_set_col_stat(IET *iet, int j, int stat)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_col_stat: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_col_stat: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* adjust non-basic status depending on the column type */
      if (!(stat == IET_BS || stat == IET_NL || stat == IET_NU ||
            stat == IET_NF || stat == IET_NS))
         fault("iet_set_col_stat: j = %d; stat = %d; invalid column sta"
            "tus", j, stat);
      if (stat != IET_BS)
      {  switch (col->type)
         {  case IET_FR:
               stat = IET_NF; break;
            case IET_LO:
               stat = IET_NL; break;
            case IET_UP:
               stat = IET_NU; break;
            case IET_DB:
               if (!(stat == IET_NL || stat == IET_NU)) stat = IET_NL;
               break;
            case IET_FX:
               stat = IET_NS; break;
            default:
               insist(col != col);
         }
      }
      /* set column status */
      col->stat = stat;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_row_locl - set link to row local extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_row_locl(IET *iet, int i, void *link);
--
-- *Description*
--
-- The routine iet_set_row_locl sets the link to an local extension of
-- i-th row of the current subproblem. */

void iet_set_row_locl(IET *iet, int i, void *link)
{     IETROW *row;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_row_locl: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_set_row_locl: i = %d; row number out of range", i);
      row = iet->row[i];
      /* set the link to local extension */
      row->link = link;
      return;
}

/*----------------------------------------------------------------------
-- iet_set_col_locl - set link to column local extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_set_col_locl(IET *iet, int j, void *link);
--
-- *Description*
--
-- The routine iet_set_col_locl sets the link to an local extension of
-- j-th column of the current subproblem. */

void iet_set_col_locl(IET *iet, int j, void *link)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_set_col_locl: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_set_col_locl: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* set the link to local extension */
      col->link = link;
      return;
}

/*----------------------------------------------------------------------
-- iet_del_rows - delete specified rows from current subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_del_rows(IET *iet, int nrs, int num[]);
--
-- *Description*
--
-- The routine iet_del_rows deletes specified rows from the current
-- subproblem. Ordinal numbers of rows to be deleted should be placed
-- in locations num[1], num[2], ..., num[nrs], where nrs > 0.
--
-- Note that deleting rows involves changing ordinal numbers of other
-- rows remaining in the current subproblem. New ordinal numbers of the
-- remaining rows can be determined with the assumption that the order
-- of rows is not changed. */

void iet_del_rows(IET *iet, int nrs, int num[])
{     IETNPD *node;
      IETRGD *rgd;
      int i, t, m_new;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_del_rows: current subproblem does not exist");
      /* mark rows to be deleted */
      if (nrs < 1)
         fault("iet_del_rows: nrs = %d; invalid parameter", nrs);
      for (t = 1; t <= nrs; t++)
      {  /* take number i of row to be deleted */
         i = num[t];
         /* obtain pointer to global descriptor of i-th row */
         if (!(1 <= i && i <= iet->m))
            fault("iet_del_rows: num[%d] = %d; row number out of range",
               t, i);
         rgd = iet->row[i]->glob;
         /* check that i-th row is not marked yet */
         if (rgd->i == 0)
            fault("iet_del_rows: num[%d] = %d; duplicate row numbers no"
               "t allowed", t, i);
         insist(rgd->i == i);
         /* clear i-th row of the constraint matrix */
         iet_set_mat_row(iet, i, 0, NULL, NULL);
         /* mark i-th row to be deleted */
         rgd->i = 0;
      }
      /* delete all marked rows and rebuild the linked list of own rows
         remaining in the current subproblem */
      m_new = 0;
      node->r_add = NULL;
      for (i = 1; i <= iet->m; i++)
      {  /* obtain pointer to global descriptor of i-th row */
         rgd = iet->row[i]->glob;
         /* check whether i-th row is marked or not */
         if (rgd->i == 0)
         {  /* i-th row is marked and should be deleted */
            if (rgd->host == node)
            {  /* it is own row of the current subproblem; delete its
                  global descriptor */
               if (iet->hook != NULL)
                  iet->hook(iet->info, IET_RD, rgd->name == NULL ? NULL
                     : get_str(iet->str_buf, rgd->name), rgd->link);
               if (rgd->name != NULL) delete_str(rgd->name);
               dmp_free_atom(iet->rgd_pool, rgd);
            }
            else
            {  /* it is non-own row inherited from some subproblem; add
                  deletion entry to the deletion list */
               IETDQE *dqe;
               dqe = dmp_get_atom(iet->dqe_pool);
               dqe->u.row = rgd;
               dqe->next = node->r_del;
               node->r_del = dqe;
            }
            /* delete local descriptor of the row */
            dmp_free_atom(iet->row_pool, iet->row[i]);
         }
         else
         {  /* i-th row is not marked and should be kept */
            /* assign new ordinal number to the row */
            rgd->i = ++m_new;
            iet->row[rgd->i] = iet->row[i];
            /* if it is own row, add it to the end of the linked list of
               own rows of the current subproblem */
            if (rgd->host == node)
            {  if (node->r_add == NULL)
               {  /* it is the first own row */
                  node->r_add = rgd;
               }
               else
               {  /* it is not the first own row */
                  IETROW *row;
                  insist(m_new > 1);
                  row = iet->row[m_new-1];
                  insist(row->glob->host == node);
                  insist(row->glob->next == NULL);
                  row->glob->next = rgd;
               }
               /* currently it is the last own row */
               rgd->next = NULL;
            }
         }
      }
      /* set new number of rows */
      iet->m = m_new;
      return;
}

/*----------------------------------------------------------------------
-- iet_del_cols - delete specified columns from current subproblem.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void iet_del_cols(IET *iet, int ncs, int num[]);
--
-- *Description*
--
-- The routine iet_del_cols deletes specified columns from the current
-- subproblem. Ordinal numbers of columns to be deleted should be placed
-- in locations num[1], num[2], ..., num[ncs], where ncs > 0.
--
-- Note that deleting columns involves changing ordinal numbers of other
-- columns remaining in the current subproblem. New ordinal numbers of
-- the remaining columns can be determined with the assumption that the
-- order of columns is not changed. */

void iet_del_cols(IET *iet, int ncs, int num[])
{     IETNPD *node;
      IETCGD *cgd;
      int j, t, n_new;
      /* obtain pointer to the current subproblem */
      node = iet->curr;
      if (node == NULL)
         fault("iet_del_cols: current subproblem does not exist");
      /* mark columns to be deleted */
      if (ncs < 1)
         fault("iet_del_cols: ncs = %d; invalid parameter", ncs);
      for (t = 1; t <= ncs; t++)
      {  /* take number j of column to be deleted */
         j = num[t];
         /* obtain pointer to global descriptor of j-th column */
         if (!(1 <= j && j <= iet->n))
            fault("iet_del_cols: num[%d] = %d; column number out of ran"
               "ge", t, j);
         cgd = iet->col[j]->glob;
         /* check that j-th column is not marked yet */
         if (cgd->j == 0)
            fault("iet_del_cols: num[%d] = %d; duplicate column numbers"
               " not allowed", t, j);
         insist(cgd->j == j);
         /* clear j-th column of the constraint matrix */
         iet_set_mat_col(iet, j, 0, NULL, NULL);
         /* mark j-th column to be deleted */
         cgd->j = 0;
      }
      /* delete all marked columns and rebuild the linked list of own
         columns remaining in the current subproblem */
      n_new = 0;
      node->c_add = NULL;
      for (j = 1; j <= iet->n; j++)
      {  /* obtain pointer to global descriptor of j-th column */
         cgd = iet->col[j]->glob;
         /* check whether j-th column is marked or not */
         if (cgd->j == 0)
         {  /* j-th column is marked and should be deleted */
            if (cgd->host == node)
            {  /* it is own column of the current subproblem; delete its
                  global descriptor */
               if (iet->hook != NULL)
                  iet->hook(iet->info, IET_CD, cgd->name == NULL ? NULL
                     : get_str(iet->str_buf, cgd->name), cgd->link);
               if (cgd->name != NULL) delete_str(cgd->name);
               dmp_free_atom(iet->cgd_pool, cgd);
            }
            else
            {  /* it is non-own column inherited from some subproblem;
                  add deletion entry to the deletion list */
               IETDQE *dqe;
               dqe = dmp_get_atom(iet->dqe_pool);
               dqe->u.col = cgd;
               dqe->next = node->c_del;
               node->c_del = dqe;
            }
            /* delete local descriptor of the column */
            dmp_free_atom(iet->col_pool, iet->col[j]);
         }
         else
         {  /* j-th column is marked and should be kept */
            /* assign new ordinal number to the column */
            cgd->j = ++n_new;
            iet->col[cgd->j] = iet->col[j];
            /* if it is own column, add it to the end of the linked list
               of own columns of the current subproblem */
            if (cgd->host == node)
            {  if (node->c_add == NULL)
               {  /* it is the first own column */
                  node->c_add = cgd;
               }
               else
               {  /* it is not the first own column */
                  IETCOL *col;
                  insist(n_new > 1);
                  col = iet->col[n_new-1];
                  insist(col->glob->host == node);
                  insist(col->glob->next == NULL);
                  col->glob->next = cgd;
               }
               /* currently it is the last own column */
               cgd->next = NULL;
            }
         }
      }
      /* set new number of columns */
      iet->n = n_new;
      return;
}

/**********************************************************************/
/* * *                SUBPROBLEM QUERYING ROUTINES                * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- iet_get_num_rows - determine number of rows.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_num_rows(IET *iet);
--
-- *Returns*
--
-- The routine iet_get_num_rows returns the number of rows in the
-- current subproblem. */

int iet_get_num_rows(IET *iet)
{     if (iet->curr == NULL)
         fault("iet_get_num_rows: current subproblem does not exist");
      return iet->m;
}

/*----------------------------------------------------------------------
-- iet_get_num_cols - determine number of columns.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_num_cols(IET *iet);
--
-- *Returns*
--
-- The routine iet_get_num_cols returns the number of columns in the
-- current subproblem. */

int iet_get_num_cols(IET *iet)
{     if (iet->curr == NULL)
         fault("iet_get_num_cols: current subproblem does not exist");
      return iet->n;
}

/*----------------------------------------------------------------------
-- iet_get_num_nz - determine number of constraint coefficients.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_num_nz(IET *iet);
--
-- *Returns*
--
-- The routine iet_get_num_nz returns the number of (non-zero) elements
-- in the constraint matrix of the current subproblem. */

int iet_get_num_nz(IET *iet)
{     if (iet->curr == NULL)
         fault("iet_get_num_nz: current subproblem does not exist");
      return iet->nz;
}

/*----------------------------------------------------------------------
-- iet_get_row_name - obtain row name.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- char *iet_get_row_name(IET *iet, int i);
--
-- *Returns*
--
-- The routine iet_get_row_name returns a pointer to the symbolic name
-- assigned to i-th row of the current subproblem. However, if the row
-- has no symbolic name assigned, the routine returns NULL. */

char *iet_get_row_name(IET *iet, int i)
{     IETRGD *rgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_row_name: current subproblem does not exist");
      /* obtain pointer to global descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_row_name: i = %d; row number out of range", i);
      rgd = iet->row[i]->glob;
      /* return pointer to the row name */
      return
         rgd->name == NULL ? NULL : get_str(iet->str_buf, rgd->name);
}

/*----------------------------------------------------------------------
-- iet_get_col_name - obtain column name.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- char *iet_get_col_name(IET *iet, int j);
--
-- *Returns*
--
-- The routine iet_get_col_name returns a pointer to the symbolic name
-- assigned to j-th column of the current subproblem. However, if the
-- column has no symbolic name assigned, the routine returns NULL. */

char *iet_get_col_name(IET *iet, int j)
{     IETCGD *cgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_col_name: current subproblem does not exist");
      /* obtain pointer to global descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_col_name: j = %d; column number out of range",
            j);
      cgd = iet->col[j]->glob;
      /* return pointer to the column name */
      return
         cgd->name == NULL ? NULL : get_str(iet->str_buf, cgd->name);
}

/*----------------------------------------------------------------------
-- iet_get_row_link - obtain link to row global extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void *iet_get_row_link(IET *iet, int i);
--
-- *Returns*
--
-- The routine iet_get_row_link returns the link to a global extension
-- of i-th row of the current subproblem. */

void *iet_get_row_link(IET *iet, int i)
{     IETRGD *rgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_row_link: current subproblem does not exist");
      /* obtain pointer to global descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_row_link: i = %d; row number out of range", i);
      rgd = iet->row[i]->glob;
      /* return the link to global extension */
      return rgd->link;
}

/*----------------------------------------------------------------------
-- iet_get_col_link - obtain link to column global extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void *iet_get_col_link(IET *iet, int j);
--
-- *Returns*
--
-- The routine iet_get_col_link returns the link to a global extension
-- of j-th column of the current subproblem. */

void *iet_get_col_link(IET *iet, int j)
{     IETCGD *cgd;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_col_link: current subproblem does not exist");
      /* obtain pointer to global descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_col_link: j = %d; column number out of range",
            j);
      cgd = iet->col[j]->glob;
      /* return the link to global extension */
      return cgd->link;
}

/*----------------------------------------------------------------------
-- iet_get_row_bnds - determine row type and bounds.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_row_bnds(IET *iet, int i, double *lb, double *ub);
--
-- *Description*
--
-- The routine iet_get_row_bnds determines the type, lower bound, and
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
--    IET_FR   -inf <  x <  +inf   free variable
--    IET_LO     lb <= x <  +inf   lower bound
--    IET_UP   -inf <  x <=  ub    upper bound
--    IET_DB     lb <= x <=  ub    double bound
--    IET_FX           x  =  lb    fixed variable
--
-- where x is the auxiliary variable associated with i-th row.
--
-- If the row has no lower bound, *lb is set to zero. If the row has no
-- upper bound, *ub is set to zero. If the row is of fixed type, *lb and
-- *ub are set to the same value.
--
-- *Returns*
--
-- The routine iet_get_row_bnds returns the type of the row. */

int iet_get_row_bnds(IET *iet, int i, double *lb, double *ub)
{     IETROW *row;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_row_bnds: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_row_bnds: i = %d; row number out of range", i);
      row = iet->row[i];
      /* store bounds of the row */
      if (lb != NULL) *lb = row->lb;
      if (ub != NULL) *ub = row->ub;
      /* and return the row type */
      return row->type;
}

/*----------------------------------------------------------------------
-- iet_get_col_bnds - determine column type and bounds.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_col_bnds(IET *iet, int j, double *lb, double *ub);
--
-- *Description*
--
-- The routine iet_get_col_bnds determines the type, lower bound, and
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
--    IET_FR   -inf <  x <  +inf   free variable
--    IET_LO     lb <= x <  +inf   lower bound
--    IET_UP   -inf <  x <=  ub    upper bound
--    IET_DB     lb <= x <=  ub    double bound
--    IET_FX           x  =  lb    fixed variable
--
-- where x is the structural variable associated with j-th column.
--
-- If the column has no lower bound, *lb is set to zero. If the column
-- has no upper bound, *ub is set to zero. If the column is of fixed
-- type, *lb and *ub are set to the same value.
--
-- *Returns*
--
-- The routine iet_get_col_bnds returns the type of the column. */

int iet_get_col_bnds(IET *iet, int j, double *lb, double *ub)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_col_bnds: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_col_bnds: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* store bounds of the column */
      if (lb != NULL) *lb = col->lb;
      if (ub != NULL) *ub = col->ub;
      /* and return the column type */
      return col->type;
}

/*----------------------------------------------------------------------
-- iet_get_obj_coef - determine objective coefficient.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- double iet_get_obj_coef(IET *iet, int j);
--
-- *Returns*
--
-- The routine iet_get_obj_coef returns objective coefficient at j-th
-- column of the current subproblem. If j is 0, the routine returns
-- constant term of the objective function of the current subproblem. */

double iet_get_obj_coef(IET *iet, int j)
{     IETCOL *col;
      double coef;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_obj_coef: current subproblem does not exist");
      if (j == 0)
      {  /* get the constant term */
         coef = iet->c0;
         goto done;
      }
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_obj_coef: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* get the objective coefficient */
      coef = col->coef;
done: return coef;
}

/*----------------------------------------------------------------------
-- iet_get_mat_row - obtain row of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_mat_row(IET *iet, int i, int ind[], double val[]);
--
-- *Description*
--
-- The routine iet_get_mat_row scans the contents of i-th row of the
-- constraint matrix of the current subproblem and stores column indices
-- and numeric values of corresponding (non-zero) elements to locations
-- ind[1], ..., ind[len] and val[1], ..., val[len], respectively, where
-- 0 <= len <= n is the number of non-zero elements in i-th row, n is
-- the number of columns in the current subproblem. If the parameter ind
-- or val is NULL, the corresponding information is not stored.
--
-- *Returns*
--
-- The routine iet_get_mat_row returns the number of non-zero elements
-- in i-th row (len). */

int iet_get_mat_row(IET *iet, int i, int ind[], double val[])
{     IETROW *row;
      IETAIJ *aij;
      int len;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_mat_row: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_mat_row: i = %d; row number out of range", i);
      row = iet->row[i];
      /* obtain i-th row of the constraint matrix */
      len = 0;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  len++;
         if (ind != NULL) ind[len] = aij->col->j;
         if (val != NULL) val[len] = aij->val;
      }
      insist(len <= iet->n);
      return len;
}

/*----------------------------------------------------------------------
-- iet_get_mat_col - obtain column of constraint matrix.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_mat_col(IET *iet, int j, int ind[], double val[]);
--
-- *Description*
--
-- The routine iet_get_mat_col scans the contents of j-th column of the
-- constraint matrix of the current subproblem and stores row indices
-- and numeric values of corresponding (non-zero) elements to locations
-- ind[1], ..., ind[len] and val[1], ..., val[len], respectively, where
-- 0 <= len <= m is the number of non-zero elements in j-th column, m is
-- the number of rows in the current subproblem. If the parameter ind or
-- val is NULL, the corresponding information is not stored.
--
-- The routine iet_get_mat_col returns the number of non-zero elements
-- in j-th column (len). */

int iet_get_mat_col(IET *iet, int j, int ind[], double val[])
{     IETCOL *col;
      IETAIJ *aij;
      int len;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_mat_col: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_mat_col: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* obtain j-th column of the constraint matrix */
      len = 0;
      for (aij = col->ptr; aij != NULL; aij = aij->c_next)
      {  len++;
         if (ind != NULL) ind[len] = aij->row->i;
         if (val != NULL) val[len] = aij->val;
      }
      insist(len <= iet->m);
      return len;
}

/*----------------------------------------------------------------------
-- iet_get_row_stat - obtain row status.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_row_stat(IET *iet, int i);
--
-- *Returns*
--
-- The routine iet_get_row_stat returns the status of i-th row in the
-- current subproblem:
--
-- IET_BS - basic row (inactive constraint);
-- IET_NL - non-basic row on its lower bound (active constraint);
-- IET_NU - non-basic row on its upper bound (active constraint);
-- IET_NF - non-basic free (unbounded) row;
-- IET_NS - non-basic fixed row (equality constraint). */

int iet_get_row_stat(IET *iet, int i)
{     IETROW *row;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_row_stat: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_row_stat: i = %d; row number out of range", i);
      row = iet->row[i];
      /* return the row status */
      return row->stat;
}

/*----------------------------------------------------------------------
-- iet_get_col_stat - obtain column status.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- int iet_get_col_stat(IET *iet, int j);
--
-- *Returns*
--
-- The routine iet_get_col_stat returns the status of j-th column in the
-- current subproblem:
--
-- IET_BS - basic column;
-- IET_NL - non-basic column on its lower bound;
-- IET_NU - non-basic column on its upper bound;
-- IET_NF - non-basic free (unbounded) column;
-- IET_NS - non-basic fixed column. */

int iet_get_col_stat(IET *iet, int j)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_col_stat: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_col_stat: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* return the column status */
      return col->stat;
}

/*----------------------------------------------------------------------
-- iet_get_row_locl - obtain link to row local extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void *iet_get_row_locl(IET *iet, int i);
--
-- *Returns*
--
-- The routine iet_get_row_locl returns the link to an local extension
-- of i-th row of the current subproblem. */

void *iet_get_row_locl(IET *iet, int i)
{     IETROW *row;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_row_locl: current subproblem does not exist");
      /* obtain pointer to local descriptor of i-th row */
      if (!(1 <= i && i <= iet->m))
         fault("iet_get_row_locl: i = %d; row number out of range", i);
      row = iet->row[i];
      /* return the link to local extension */
      return row->link;
}

/*----------------------------------------------------------------------
-- iet_get_col_locl - obtain link to column local extension.
--
-- *Synopsis*
--
-- #include "glpiet.h"
-- void *iet_get_col_locl(IET *iet, int j);
--
-- *Returns*
--
-- The routine iet_get_col_locl returns the link to an local extension
-- of j-th column of the current subproblem. */

void *iet_get_col_locl(IET *iet, int j)
{     IETCOL *col;
      /* the current subproblem must exist */
      if (iet->curr == NULL)
         fault("iet_get_col_locl: current subproblem does not exist");
      /* obtain pointer to local descriptor of j-th column */
      if (!(1 <= j && j <= iet->n))
         fault("iet_get_col_locl: j = %d; column number out of range",
            j);
      col = iet->col[j];
      /* return the link to local extension */
      return col->link;
}

/* eof */

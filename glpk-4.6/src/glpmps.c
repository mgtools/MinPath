/* glpmps.c */

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
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "glpavl.h"
#include "glplib.h"
#include "glpmps.h"

/*----------------------------------------------------------------------
-- load_mps - load linear programming model in MPS format.
--
-- *Synopsis*
--
-- #include "glpmps.h"
-- MPS *load_mps(char *fname);
--
-- *Description*
--
-- The load_mps routine loads a linear programming model in the MPS
-- format from the text file whose name is the character string fname.
--
-- Detailed description of the MPS format can be found, for example,
-- in the following book:
--
-- B.A.Murtagh. Advanced Linear Programming: Computation and Practice.
-- McGraw-Hill, 1981.
--
-- *Returns*
--
-- The load_mps routine returns a pointer to the object of MPS type
-- that represents the loaded data in the MPS format. In case of error
-- the routine prints an appropriate error message and returns NULL. */

struct dsa
{     /* dynamic storage area */
      MPS *mps;
      /* pointer to MPS data block */
      AVLTREE *t_row, *t_col, *t_rhs, *t_rng, *t_bnd;
      /* symbol tables for rows, columns, right-hand sides, ranges, and
         bounds, respectively */
      char *fname;
      /* name of the input text file */
      FILE *fp;
      /* stream assigned to the input text file */
      int seqn;
      /* card sequential number */
      char card[80+1];
      /* card image buffer */
      char f1[2+1], f2[31+1], f3[31+1], f4[12+1], f5[31+1], f6[12+1];
      /* splitted card fields */
};

static char *unknown = "UNKNOWN";
/* default name */

static int read_card(struct dsa *dsa);
/* read the next card */

static int split_card(struct dsa *dsa);
/* split data card to separate fileds */

static int load_rows(struct dsa *dsa);
/* load ROWS section */

static int load_columns(struct dsa *dsa, AVLTREE *t_xxx);
/* load COLUMNS, RHS, or RANGES section */

static int load_bounds(struct dsa *dsa);
/* load BOUNDS section */

static int load_quadobj(struct dsa *dsa);
/* load QUADOBJ section */

MPS *load_mps(char *_fname)
{     AVLNODE *node;
      MPSCOL *col; MPSCQE *cptr, *cqe;
      MPSBND *bnd; MPSBQE *bptr, *bqe;
      MPSQFE *qptr, *qfe;
      struct dsa _dsa, *dsa = &_dsa;
      /* initialization */
      dsa->mps = NULL;
      dsa->t_row = dsa->t_col = dsa->t_rhs = dsa->t_rng = dsa->t_bnd =
         NULL;
      dsa->fname = _fname;
      print("load_mps: reading LP data from `%s'...", dsa->fname);
      dsa->fp = ufopen(dsa->fname, "r");
      if (dsa->fp == NULL)
      {  print("load_mps: unable to open `%s' - %s", dsa->fname,
            strerror(errno));
         goto fail;
      }
      dsa->seqn = 0;
      dsa->mps = umalloc(sizeof(MPS));
      memset(dsa->mps, 0, sizeof(MPS));
      dsa->mps->pool = dmp_create_pool(0);
      dsa->t_row = avl_create_tree(NULL, avl_strcmp);
      dsa->t_col = avl_create_tree(NULL, avl_strcmp);
      dsa->t_rhs = avl_create_tree(NULL, avl_strcmp);
      dsa->t_rng = avl_create_tree(NULL, avl_strcmp);
      dsa->t_bnd = avl_create_tree(NULL, avl_strcmp);
      /* process NAME indicator card */
      if (read_card(dsa)) goto fail;
      if (memcmp(dsa->card, "NAME ", 5) != 0)
      {  print("%s:%d: NAME indicator card missing", dsa->fname,
            dsa->seqn);
         goto fail;
      }
      memcpy(dsa->f3, dsa->card+14, 8);
      dsa->f3[8] = '\0'; strspx(dsa->f3);
      if (dsa->f3[0] == '\0') strcpy(dsa->f3, unknown);
      dsa->mps->name = dmp_get_atomv(dsa->mps->pool, strlen(dsa->f3)+1);
      strcpy(dsa->mps->name, dsa->f3);
      print("load_mps: name `%s'", dsa->mps->name);
      /* process ROWS section */
      if (read_card(dsa)) goto fail;
      if (memcmp(dsa->card, "ROWS ", 5) != 0)
      {  print("%s:%d: ROWS indicator card missing", dsa->fname,
            dsa->seqn);
         goto fail;
      }
      if (load_rows(dsa)) goto fail;
      dsa->mps->n_row = dsa->t_row->size;
      print("load_mps: %d rows", dsa->mps->n_row);
      /* process COLUMNS section */
      if (memcmp(dsa->card, "COLUMNS ", 8) != 0)
      {  print("%s:%d: COLUMNS indicator card missing", dsa->fname,
            dsa->seqn);
         goto fail;
      }
      if (load_columns(dsa, dsa->t_col)) goto fail;
      dsa->mps->n_col = dsa->t_col->size;
      print("load_mps: %d columns", dsa->mps->n_col);
      /* count non-zeros */
      {  int nz = 0;
         for (node = avl_find_next_node(dsa->t_col, NULL); node != NULL;
            node = avl_find_next_node(dsa->t_col, node))
         {  col = node->link;
            for (cqe = col->ptr; cqe != NULL; cqe = cqe->next) nz++;
         }
         print("load_mps: %d non-zeros", nz);
      }
      /* process RHS section */
      if (memcmp(dsa->card, "RHS ", 4) == 0)
         if (load_columns(dsa, dsa->t_rhs)) goto fail;
      dsa->mps->n_rhs = dsa->t_rhs->size;
      print("load_mps: %d right-hand side vector(s)", dsa->mps->n_rhs);
      /* process RANGES section */
      if (memcmp(dsa->card, "RANGES ", 7) == 0)
         if (load_columns(dsa, dsa->t_rng)) goto fail;
      dsa->mps->n_rng = dsa->t_rng->size;
      print("load_mps: %d range vector(s)", dsa->mps->n_rng);
      /* process BOUNDS section */
      if (memcmp(dsa->card, "BOUNDS ", 7) == 0)
         if (load_bounds(dsa)) goto fail;
      dsa->mps->n_bnd = dsa->t_bnd->size;
      print("load_mps: %d bound vector(s)", dsa->mps->n_bnd);
      /* process QUADOBJ section */
      dsa->mps->quad = NULL;
      if (memcmp(dsa->card, "QUADOBJ ", 8) == 0)
      {  int count = 0;
         if (load_quadobj(dsa)) goto fail;
         for (qfe = dsa->mps->quad; qfe != NULL; qfe = qfe->next)
            count++;
         print("load_mps: %d quadratic form elements", count);
      }
      /* process ENDATA indicator card */
      if (memcmp(dsa->card, "ENDATA ", 7) != 0)
      {  print("%s:%d: invalid indicator card", dsa->fname, dsa->seqn);
         goto fail;
      }
      print("load_mps: %d cards were read", dsa->seqn);
      ufclose(dsa->fp);
      /* build row list */
      dsa->mps->row = ucalloc(1+dsa->mps->n_row, sizeof(MPSROW *));
      for (node = avl_find_next_node(dsa->t_row, NULL); node != NULL;
         node = avl_find_next_node(dsa->t_row, node))
         dsa->mps->row[node->type] = node->link;
      avl_delete_tree(dsa->t_row);
      /* build column list and restore original order of elements */
      dsa->mps->col = ucalloc(1+dsa->mps->n_col, sizeof(MPSCOL *));
      for (node = avl_find_next_node(dsa->t_col, NULL); node != NULL;
         node = avl_find_next_node(dsa->t_col, node))
      {  col = node->link; cptr = NULL;
         while (col->ptr != NULL)
         {  cqe = col->ptr;
            col->ptr = cqe->next;
            cqe->next = cptr;
            cptr = cqe;
         }
         col->ptr = cptr;
         dsa->mps->col[node->type] = col;
      }
      avl_delete_tree(dsa->t_col);
      /* build rhs list and restore original order of elements */
      dsa->mps->rhs = ucalloc(1+dsa->mps->n_rhs, sizeof(MPSCOL *));
      for (node = avl_find_next_node(dsa->t_rhs, NULL); node != NULL;
         node = avl_find_next_node(dsa->t_rhs, node))
      {  col = node->link; cptr = NULL;
         while (col->ptr != NULL)
         {  cqe = col->ptr;
            col->ptr = cqe->next;
            cqe->next = cptr;
            cptr = cqe;
         }
         col->ptr = cptr;
         dsa->mps->rhs[node->type] = col;
      }
      avl_delete_tree(dsa->t_rhs);
      /* build ranges list and restore original order of elements */
      dsa->mps->rng = ucalloc(1+dsa->mps->n_rng, sizeof(MPSCOL *));
      for (node = avl_find_next_node(dsa->t_rng, NULL); node != NULL;
         node = avl_find_next_node(dsa->t_rng, node))
      {  col = node->link; cptr = NULL;
         while (col->ptr != NULL)
         {  cqe = col->ptr;
            col->ptr = cqe->next;
            cqe->next = cptr;
            cptr = cqe;
         }
         col->ptr = cptr;
         dsa->mps->rng[node->type] = col;
      }
      avl_delete_tree(dsa->t_rng);
      /* build bounds list and restore original order of elements */
      dsa->mps->bnd = ucalloc(1+dsa->mps->n_bnd, sizeof(MPSBND *));
      for (node = avl_find_next_node(dsa->t_bnd, NULL); node != NULL;
         node = avl_find_next_node(dsa->t_bnd, node))
      {  bnd = node->link; bptr = NULL;
         while (bnd->ptr != NULL)
         {  bqe = bnd->ptr;
            bnd->ptr = bqe->next;
            bqe->next = bptr;
            bptr = bqe;
         }
         bnd->ptr = bptr;
         dsa->mps->bnd[node->type] = bnd;
      }
      avl_delete_tree(dsa->t_bnd);
      /* restore original order of quadratic form elements */
      qptr = NULL;
      while (dsa->mps->quad != NULL)
      {  qfe = dsa->mps->quad;
         dsa->mps->quad = qfe->next;
         qfe->next = qptr;
         qptr = qfe;
      }
      dsa->mps->quad = qptr; 
      /* loading has been completed */
      return dsa->mps;
fail: /* something wrong in Danish kingdom */
      if (dsa->mps != NULL)
      {  if (dsa->mps->pool != NULL) dmp_delete_pool(dsa->mps->pool);
         if (dsa->mps->row != NULL) ufree(dsa->mps->row);
         if (dsa->mps->col != NULL) ufree(dsa->mps->col);
         if (dsa->mps->rhs != NULL) ufree(dsa->mps->rhs);
         if (dsa->mps->rng != NULL) ufree(dsa->mps->rng);
         if (dsa->mps->bnd != NULL) ufree(dsa->mps->bnd);
         ufree(dsa->mps);
      }
      if (dsa->t_row != NULL) avl_delete_tree(dsa->t_row);
      if (dsa->t_col != NULL) avl_delete_tree(dsa->t_col);
      if (dsa->t_rhs != NULL) avl_delete_tree(dsa->t_rhs);
      if (dsa->t_rng != NULL) avl_delete_tree(dsa->t_rng);
      if (dsa->t_bnd != NULL) avl_delete_tree(dsa->t_bnd);
      if (dsa->fp != NULL) ufclose(dsa->fp);
      return NULL;
}

/*----------------------------------------------------------------------
-- read_card - read the next card.
--
-- This routine reads the next 80-column card from the input text file
-- and places its image into the character string card. If the card was
-- read successfully, the routine returns zero, otherwise non-zero. */

static int read_card(struct dsa *dsa)
{     int k, c;
loop: dsa->seqn++;
      memset(dsa->card, ' ', 80), dsa->card[80] = '\0';
      k = 0;
      for (;;)
      {  c = fgetc(dsa->fp);
         if (ferror(dsa->fp))
         {  print("%s:%d: read error - %s", dsa->fname, dsa->seqn,
               strerror(errno));
            return 1;
         }
         if (feof(dsa->fp))
         {  if (k == 0)
               print("%s:%d: unexpected eof", dsa->fname, dsa->seqn);
            else
               print("%s:%d: missing final LF", dsa->fname, dsa->seqn);
            return 1;
         }
         if (c == '\r') continue;
         if (c == '\n') break;
         if (iscntrl(c))
         {  print("%s:%d: invalid control character 0x%02X", dsa->fname,
               dsa->seqn, c);
            return 1;
         }
         if (k == 80)
         {  print("%s:%d: card image too long", dsa->fname, dsa->seqn);
            return 1;
         }
         dsa->card[k++] = (char)c;
      }
      /* asterisk in the leftmost column means comment */
      if (dsa->card[0] == '*') goto loop;
      return 0;
}

/*----------------------------------------------------------------------
-- split_card - split data card to separate fields.
--
-- This routine splits the current data card to separate fileds f1, f2,
-- f3, f4, f5, and f6. If the data card has correct format, the routine
-- returns zero, otherwise non-zero. */

static int split_card(struct dsa *dsa)
{     /* col. 1: blank */
      if (memcmp(dsa->card+0, " ", 1))
fail: {  print("%s:%d: invalid data card", dsa->fname, dsa->seqn);
         return 1;
      }
      /* col. 2-3: field 1 (code) */
      memcpy(dsa->f1, dsa->card+1, 2);
      dsa->f1[2] = '\0'; strspx(dsa->f1);
      /* col. 4: blank */
      if (memcmp(dsa->card+3, " ", 1)) goto fail;
      /* col. 5-12: field 2 (name) */
      memcpy(dsa->f2, dsa->card+4, 8);
      dsa->f2[8] = '\0'; strspx(dsa->f2);
      /* col. 13-14: blanks */
      if (memcmp(dsa->card+12, "  ", 2)) goto fail;
      /* col. 15-22: field 3 (name) */
      memcpy(dsa->f3, dsa->card+14, 8);
      dsa->f3[8] = '\0'; strspx(dsa->f3);
      if (dsa->f3[0] == '$')
      {  /* from col. 15 to the end of the dsa->card is a comment */
         dsa->f3[0] = dsa->f4[0] = dsa->f5[0] = dsa->f6[0] = '\0';
         goto done;
      }
      /* col. 23-24: blanks */
      if (memcmp(dsa->card+22, "  ", 2)) goto fail;
      /* col. 25-36: field 4 (number) */
      memcpy(dsa->f4, dsa->card+24, 12);
      dsa->f4[12] = '\0'; strspx(dsa->f4);
      /* col. 37-39: blanks */
      if (memcmp(dsa->card+36, "   ", 3)) goto fail;
      /* col. 40-47: field 5 (name) */
      memcpy(dsa->f5, dsa->card+39,  8);
      dsa->f5[8]  = '\0'; strspx(dsa->f5);
      if (dsa->f5[0] == '$')
      {  /* from col. 40 to the end of the card is a comment */
         dsa->f5[0] = dsa->f6[0] = '\0';
         goto done;
      }
      /* col. 48-49: blanks */
      if (memcmp(dsa->card+47, "  ", 2)) goto fail;
      /* col. 50-61: field 6 (number) */
      memcpy(dsa->f6, dsa->card+49, 12);
      dsa->f6[12] = '\0'; strspx(dsa->f6);
      /* col. 62-71: blanks */
      if (memcmp(dsa->card+61, "          ", 10)) goto fail;
done: return 0;
}

/*----------------------------------------------------------------------
-- load_rows - load ROWS section.
--
-- The load_rows routine loads ROWS section reading data cards that
-- are placed after ROWS indicator card. If loading is ok, the routine
-- returns zero, otherwise non-zero. */

static int load_rows(struct dsa *dsa)
{     MPSROW *row;
      AVLNODE *node;
loop: /* process the next data card */
      if (read_card(dsa)) return 1;
      if (dsa->card[0] != ' ') goto done;
      if (split_card(dsa)) return 1;
      if (strcmp(dsa->f3, "'MARKER'") == 0)
      {  print("%s:%d: invalid use of marker in ROWS section",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f1[0] == '\0')
      {  print("%s:%d: missing row type in field 1", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (!(strchr("NGLE", dsa->f1[0]) != NULL && dsa->f1[1] == '\0'))
      {  print("%s:%d: unknown row type `%s' in field 1", dsa->fname,
            dsa->seqn, dsa->f1);
         return 1;
      }
      if (dsa->f2[0] == '\0')
      {  print("%s:%d: missing row name in field 2", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (dsa->f3[0] != '\0' || dsa->f4[0] != '\0' || dsa->f5[0] != '\0'
         || dsa->f6[0] != '\0')
      {  print("%s:%d: invalid data in fields 3-6", dsa->fname,
            dsa->seqn);
         return 1;
      }
      /* create new row */
      if (avl_find_by_key(dsa->t_row, dsa->f2) != NULL)
      {  print("%s:%d: row `%s' multiply specified", dsa->fname,
            dsa->seqn, dsa->f2);
         return 1;
      }
      row = dmp_get_atomv(dsa->mps->pool, sizeof(MPSROW));
      row->name = dmp_get_atomv(dsa->mps->pool, strlen(dsa->f2)+1);
      strcpy(row->name, dsa->f2);
      strcpy(row->type, dsa->f1);
      /* add row name to the symbol table */
      node = avl_insert_by_key(dsa->t_row, row->name);
      node->type = dsa->t_row->size;
      node->link = row;
      goto loop;
done: return 0;
}

/*----------------------------------------------------------------------
-- load_columns - load COLUMNS, RHS, or RANGES section.
--
-- The load_columns routine loads COLUMNS, RHS, or RANGES section (that
-- depends upon the parameter t_xxx) reading data cards that are placed
-- after the corresponding indicator card. If loading is ok, the routine
-- return zero, otherwise non-zero. */

static int load_columns(struct dsa *dsa, AVLTREE *t_xxx)
{     MPSCOL *col;
      MPSCQE *cqe;
      AVLNODE *node, *ref;
      char name[31+1]; double val;
      int flag = 0;
      strcpy(name, (t_xxx == dsa->t_col) ? "" : unknown);
loop: /* process the next data card */
      if (read_card(dsa)) return 1;
      if (dsa->card[0] != ' ') goto done;
      if (split_card(dsa)) return 1;
      /* process optional INTORG/INTEND markers */
      if (strcmp(dsa->f3, "'MARKER'") == 0)
      {  if (t_xxx != dsa->t_col)
         {  print("%s:%d): invalid use of marker in RHS or RANGES secti"
               "on", dsa->fname, dsa->seqn);
            return 1;
         }
         if (!(dsa->f1[0] == '\0' && dsa->f4[0] == '\0' && dsa->f6[0]
            == '\0'))
         {  print("%s:%d: invalid data in fields 1, 4, or 6",
               dsa->fname, dsa->seqn);
            return 1;
         }
         if (dsa->f2[0] == '\0')
         {  print("%s:%d: missing marker name in field 2", dsa->fname,
               dsa->seqn);
            return 1;
         }
         if (strcmp(dsa->f5, "'INTORG'") == 0)
            flag = 1;
         else if (strcmp(dsa->f5, "'INTEND'") == 0)
            flag = 0;
         else
         {  print("%s:%d: unknown marker in field 5", dsa->fname,
               dsa->seqn);
            return 1;
         }
         goto skip;
      }
      /* process the data card */
      if (dsa->f1[0] != '\0')
      {  print("%s:%d: invalid data in field 1", dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f2[0] == '\0') strcpy(dsa->f2, name);
      if (dsa->f2[0] == '\0')
      {  print("%s:%d: missing column name in field 2", dsa->fname,
            dsa->seqn);
         return 1;
      }
      strcpy(name, dsa->f2);
      /* search for column or vector specified in field 2 */
      node = avl_find_by_key(t_xxx, dsa->f2);
      if (node == NULL)
      {  /* not found; create new column or vector */
         col = dmp_get_atomv(dsa->mps->pool, sizeof(MPSCOL));
         col->name = dmp_get_atomv(dsa->mps->pool, strlen(dsa->f2)+1);
         strcpy(col->name, dsa->f2);
         col->flag = 0;
         col->ptr = NULL;
         /* add column or vector name to the symbol table */
         node = avl_insert_by_key(t_xxx, col->name);
         node->type = t_xxx->size;
         node->link = col;
      }
      col->flag = flag;
#if 1
      /* all elements of the same column or vector should be placed
         together as specified by MPS format, although such restriction
         is not essential for this routine */
      if (node->type < t_xxx->size)
      {  print("%s:%d: %s `%s' multiply specified",
            dsa->fname, dsa->seqn,
            t_xxx == dsa->t_col ? "column" :
            t_xxx == dsa->t_rhs ? "right-hand side vector" :
            t_xxx == dsa->t_rng ? "range vector" : "???", dsa->f2);
         return 1;
      }
#endif
      /* process the first row-element pair (fields 3 and 4) */
      if (dsa->f3[0] == '\0')
      {  print("%s:%d: missing row name in field 3", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (dsa->f4[0] == '\0')
      {  print("%s:%d: missing value in field 4", dsa->fname,
            dsa->seqn);
         return 1;
      }
      /* create new column or vector element */
      ref = avl_find_by_key(dsa->t_row, dsa->f3);
      if (ref == NULL)
      {  print("%s:%d: row `%s' not found", dsa->fname, dsa->seqn,
            dsa->f3);
         return 1;
      }
      if (str2dbl(dsa->f4, &val))
      {  print("%s:%d: invalid value `%s'", dsa->fname, dsa->seqn,
            dsa->f4);
         return 1;
      }
      cqe = dmp_get_atomv(dsa->mps->pool, sizeof(MPSCQE));
      cqe->ind = ref->type;
      cqe->val = val;
      cqe->next = col->ptr;
      col->ptr = cqe;
      /* process the second row-element pair (fields 5 and 6) */
      if (dsa->f5[0] == '\0' && dsa->f6[0] == '\0') goto skip;
      if (dsa->f5[0] == '\0')
      {  print("%s:%d: missing row name in field 5", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (dsa->f6[0] == '\0')
      {  print("%s:%d: missing value in filed 6", dsa->fname,
            dsa->seqn);
         return 1;
      }
      /* create new column or vector element */
      ref = avl_find_by_key(dsa->t_row, dsa->f5);
      if (ref == NULL)
      {  print("%s:%d: row `%s' not found", dsa->fname, dsa->seqn,
            dsa->f5);
         return 1;
      }
      if (str2dbl(dsa->f6, &val))
      {  print("%s:%d: invalid value `%s'", dsa->fname, dsa->seqn,
            dsa->f6);
         return 1;
      }
      cqe = dmp_get_atomv(dsa->mps->pool, sizeof(MPSCQE));
      cqe->ind = ref->type;
      cqe->val = val;
      cqe->next = col->ptr;
      col->ptr = cqe;
skip: goto loop;
done: return 0;
}

/*----------------------------------------------------------------------
-- load_bounds - load BOUNDS section.
--
-- The load_bounds routine loads BOUNDS section reading data cards that
-- are placed after BOUNDS indicator card. If loading is ok, the routine
-- returns zero, otherwise non-zero. */

static int load_bounds(struct dsa *dsa)
{     MPSBND *bnd;
      MPSBQE *bqe;
      AVLNODE *node, *ref;
      char name[31+1]; double val;
      strcpy(name, unknown);
loop: /* process the next data card */
      if (read_card(dsa)) return 1;
      if (dsa->card[0] != ' ') goto done;
      if (split_card(dsa)) return 1;
      if (strcmp(dsa->f3, "'MARKER'") == 0)
      {  print("%s:%d: invalid use of marker in BOUNDS section",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f1[0] == '\0')
      {  print("%s:%d: missing bound type in field 1", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (strcmp(dsa->f1, "LO") && strcmp(dsa->f1, "UP") &&
          strcmp(dsa->f1, "FX") && strcmp(dsa->f1, "FR") &&
          strcmp(dsa->f1, "MI") && strcmp(dsa->f1, "PL") &&
          strcmp(dsa->f1, "UI") && strcmp(dsa->f1, "BV"))
      {  print("%s:%d: unknown bound type `%s' in field 1", dsa->fname,
            dsa->seqn, dsa->f1);
         return 1;
      }
      if (dsa->f2[0] == '\0') strcpy(dsa->f2, name);
      strcpy(name, dsa->f2);
      if (dsa->f3[0] == '\0')
      {  print("%s:%d: missing column name in field 3", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (dsa->f4[0] == '\0')
      if (!strcmp(dsa->f1, "LO") || !strcmp(dsa->f1, "UP") ||
          !strcmp(dsa->f1, "FX") || !strcmp(dsa->f1, "UI"))
      {  print("%s:%d: missing value in field 4", dsa->fname,
            dsa->seqn);
         return 1;
      }
      if (dsa->f5[0] != '\0' && dsa->f6[0] != '\0')
      {  print("%s:%d: invalid data in field 5-6", dsa->fname,
            dsa->seqn);
         return 1;
      }
      /* search for bound vector specified in field 2 */
      node = avl_find_by_key(dsa->t_bnd, dsa->f2);
      if (node == NULL)
      {  /* not found; create new bound vector */
         bnd = dmp_get_atomv(dsa->mps->pool, sizeof(MPSBND));
         bnd->name = dmp_get_atomv(dsa->mps->pool, strlen(dsa->f2)+1);
         strcpy(bnd->name, dsa->f2);
         bnd->ptr = NULL;
         /* add vector name to the symbol table */
         node = avl_insert_by_key(dsa->t_bnd, bnd->name);
         node->type = dsa->t_bnd->size;
         node->link = bnd;
      }
#if 1
      /* all elements of the same bound vector should be placed
         together as specified by MPS format, although such restriction
         is not essential for this routine */
      if (node->type < dsa->t_bnd->size)
      {  print("%s:%d: bound vector `%s' multiply specified",
            dsa->fname, dsa->seqn, dsa->f2);
         return 1;
      }
#endif
      /* process column-element pair */
      ref = avl_find_by_key(dsa->t_col, dsa->f3);
      if (ref == NULL)
      {  print("%s:%d: column `%s' not found", dsa->fname, dsa->seqn,
            dsa->f3);
         return 1;
      }
      val = 0.0;
      if (dsa->f4[0] != '\0' && str2dbl(dsa->f4, &val))
      {  print("%s:%d: invalid value `%s'", dsa->fname, dsa->seqn,
            dsa->f4);
         return 1;
      }
      /* create new bound vector element */
      bqe = dmp_get_atomv(dsa->mps->pool, sizeof(MPSBQE));
      strcpy(bqe->type, dsa->f1);
      bqe->ind = ref->type;
      bqe->val = val;
      bqe->next = bnd->ptr;
      bnd->ptr = bqe;
      goto loop;
done: return 0;
}

/*----------------------------------------------------------------------
-- load_quadobj - load QUADOBJ section.
--
-- The load_quadobj routine loads QUADOBJ section reading data cards
-- that are placed after QUADOBJ indicator card. If loading is ok, the
-- routine returns zero, otherwise non-zero.
--
-- The QUADOBJ section specifies quadratic part x'*Q*x of the objective
-- function. Should note that this feature is non-standard extension of
-- MPS format. For detailed format description see:
--
-- I.Maros, C.Meszaros. A Repository of Convex Quadratic Programming
-- Problems. */

static int load_quadobj(struct dsa *dsa)
{     MPSQFE *qfe;
      AVLNODE *ref1, *ref2;
      double val;
loop: /* process the next data card */
      if (read_card(dsa)) return 1;
      if (dsa->card[0] != ' ') goto done;
      if (split_card(dsa)) return 1;
      if (strcmp(dsa->f3, "'MARKER'") == 0)
      {  print("%s:%d: invalid use of marker on QUADOBJ section",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (!(dsa->f1[0] == '\0' && dsa->f5[0] == '\0' &&
            dsa->f6[0] == '\0'))
      {  print("%s:%d: invalid data in fields 1, 5, or 6",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f2[0] == '\0')
      {  print("%s:%d: missing first column name in field 2",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f3[0] == '\0')
      {  print("%s:%d: missing second column name in field 3",
            dsa->fname, dsa->seqn);
         return 1;
      }
      if (dsa->f4[0] == '\0')
      {  print("%s:%d: missing value in field 4", dsa->fname,
            dsa->seqn);
         return 1;
      }
      ref1 = avl_find_by_key(dsa->t_col, dsa->f2);
      if (ref1 == NULL)
      {  print("%s:%d: column `%s' not found", dsa->fname, dsa->seqn,
            dsa->f2);
         return 1;
      }
      ref2 = avl_find_by_key(dsa->t_col, dsa->f3);
      if (ref2 == NULL)
      {  print("%s:%d: column `%s' not found", dsa->fname, dsa->seqn,
            dsa->f3);
         return 1;
      }
      val = 0.0;
      if (dsa->f4[0] != '\0' && str2dbl(dsa->f4, &val))
      {  print("%s:%d: invalid value `%s'", dsa->fname, dsa->seqn,
            dsa->f4);
         return 1;
      }
      /* create new quadratic form element */
      qfe = dmp_get_atomv(dsa->mps->pool, sizeof(MPSQFE));
      qfe->ind1 = ref1->type;
      qfe->ind2 = ref2->type;
      qfe->val = val;
      qfe->next = dsa->mps->quad;
      dsa->mps->quad = qfe;
      goto loop;
done: return 0;
}

/*----------------------------------------------------------------------
-- free_mps - free linear programming model in MPS format.
--
-- *Synopsis*
--
-- #include "glpmps.h"
-- void free_mps(MPS *mps);
--
-- *Description*
--
-- The free_mps routine frees all memory allocated to the object of MPS
-- type that represents linear programming model in the MPS format. */

void free_mps(MPS *mps)
{       dmp_delete_pool(mps->pool);
        ufree(mps->row);
        ufree(mps->col);
        ufree(mps->rhs);
        ufree(mps->rng);
        ufree(mps->bnd);
        ufree(mps);
        return;
}

/* eof */

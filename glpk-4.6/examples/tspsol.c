/* tspsol.c */

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
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "glpk.h"
#include "glptsp.h"

/*----------------------------------------------------------------------
-- This program is a stand-alone solver intended for solving Symmetric
-- Traveling Salesman Problem (TSP) using the branch-and-bound method.
--
-- Note that this program is only an illustrative example. It is *not*
-- a state-of-the-art code, therefore only TSP instances of small size
-- (perhaps not more than 100 cities) can be solved using this code.
--
-- To run this program use the following command:
--
--    tspsol tsp-file
--
-- where tsp-file specifies an input text file containing TSP data.
--
-- Detailed description of the input format recognized by this program
-- is given in the report: Gerhard Reinelt, "TSPLIB 95". This report as
-- well as TSPLIB, a library of sample TSP instances (and other related
-- problems), are freely available for research purposes at the webpage
-- <http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/>.
-- See also two examples sample.tsp and dantzig.tsp which are included
-- in the GLPK distribution (the latter is from TSPLIB).
--
-- Symmetric Traveling Salesman Problem
-- ------------------------------------
--
-- Let the complete undirected graph be given:
--
--    K = (V, E),                                                    (1)
--
-- where V = {1, ..., n} is the set of nodes, E = V cross V is the set
-- of edges. Let also each edge e = (i,j) be assigned a positive number
-- c[i,j], which is the length of e. The Symmetric Traveling Salesman
-- Problem (TSP) is to find a tour in K of minimal length.
--
-- Integer programming description of TSP
-- --------------------------------------
--
-- For a set of nodes W within V introduce the following notation:
--
--    d(W) = {(i,j):i in W and j not in W or i not in W and j in W}, (2)
--
-- i.e. d(W) is the set of edges which have exactly one endnode in W.
-- If W = {v}, i.e. W consists of the only node, we write simply d(v).
--
-- The integer programming description of TSP is the following:
--
--    minimize        sum c[i,j] * x[i,j]                            (3)
--                    i,j
--
--    subject to      sum      x[i,j]  = 2      for all v in V       (4)
--               (i,j) in d(v)
--
--                    sum      x[i,j] >= 2      for all W within V,  (5)
--               (i,j) in d(W)                  W != empty, W != V
--
--               x[i,j] in {0, 1}               for all i, j         (6)
--
-- The binary variables x[i,j] have conventional meaning: if x[i,j] = 1,
-- the edge (i,j) is included in the tour, otherwise, if x[i,j] = 0, the
-- edge is not included in the tour.
--
-- The constraints (4) are called degree constraints. They require that
-- for each node v in V there must be exactly two edges included in the
-- tour which are incident to v.
--
-- The constraints (5) are called subtour elimination constraints. They
-- are intended to forbid subtours. */

static char *in_file = NULL;
/* name of input text file in TSPLIB format */

static char *out_sol = NULL;
/* name of output text file to write basic solution of LP relaxation
   in plain text format whenever a better integer feasible solution has
   been found by the solver */

static int trace = 0;
/* if this flag is set, detailed output is produced */

static TSP *tsp;
/* TSP instance to be solved */

static int nodes;
/* total number of nodes in the problem */

/* Application extensions of rows and columns are the following:

   Mark           Link                 Origin
   ------------   ------------------   -------------------------------
   DEGREE_MARK    NULL                 Each row that corresponds to a
                                       degree constraint

   SUBTOUR_MARK   Pointer to SUBTOUR   Each row that corresponds to a
                                       subtour elimination constraint

   EDGE_MARK      Pointer to EDGE      Each column that corresponds to
                                       an edge of the graph */

#define DEGREE_MARK  'D'
#define SUBTOUR_MARK 'S'
#define EDGE_MARK    'E'

typedef struct SUBTOUR SUBTOUR;
typedef struct EDGE EDGE;

struct SUBTOUR
{     /* each subtour elimination constraint corresponds to some cut
         (W, V \ W), where W != empty, W != V; the complete list of
         nodes in W is stored as the linked list of SUBTOUR entries,
         and a pointer to the first entry is stored in the row link;
         this list is needed to determine constraint coefficients for
         columns having been generated */
      int i;
      /* number of node which belongs to set W */
      SUBTOUR *next;
      /* pointer to the next SUBTOUR entry for the same cut */
};

struct EDGE
{     /* each edge is identified by its two end-nodes i and j, where
         1 <= i < j <= nodes; pointer to the EDGE entry is stored in
         the column link */
      int i;
      /* number of one node incident to corresponding edge */
      int j;
      /* number of another node incident to corresponding edge */
};

static DMP *subtour_pool;
/* memory pool for SUBTOUR entries */

static DMP *edge_pool;
/* memory pool for EDGE entries */

/*----------------------------------------------------------------------
-- greedy_tour - find salesman's tour using greedy heuristic.
--
-- This routine finds a salesman's tour using greedy heuristic. Node
-- numbers are stored in locations tour[1], ..., tour[nn] in the order,
-- in which they have to be visited, i.e. the tour built by the routine
-- is: tour[1] -> tour[2] -> ... -> tour[nodes] -> tour[1], where nodes
-- is the number of nodes in the graph. The tour length is returned by
-- the routine on exit. */

static double greedy_tour(TSP *tsp, int tour[])
{     int i, j, k, d, dmin, *flag;
      double sum;
      /* flag[i] means that i-th node has been visited */
      flag = ucalloc(1+nodes, sizeof(int));
      for (i = 1; i <= nodes; i++) flag[i] = 0;
      tour[1] = 1, flag[1] = 1;
      for (i = 2; i <= nodes; i++)
      {  k = 0, dmin = INT_MAX;
         for (j = 2; j <= nodes; j++)
         {  if (!flag[j])
            {  d = tsp_distance(tsp, tour[i-1], j);
               if (dmin > d) k = j, dmin = d;
            }
         }
         insist(k != 0);
         tour[i] = k, flag[k] = 1;
      }
      ufree(flag);
      /* compute length of the tour found */
      sum = 0.0;
      for (i = 1; i <= nodes; i++)
      {  d = tsp_distance(tsp, tour[i], tour[i < nodes ? i+1 : 1]);
         sum += (double)d;
      }
      return sum;
}

/*----------------------------------------------------------------------
-- initialize - build root subproblem.
--
-- This routine builds the root subproblem, which initially includes
-- all degree constraints and edges (columns) for an initial tour found
-- by the greedy heuristic. */

static void initialize(IOS *ios)
{     EDGE *edge;
      int i, j, k, temp, *tour, ind[1+2];
      double sum, val[1+2];
      char name[255+1];
      /* add degree constraints */
      for (i = 1; i <= nodes; i++)
      {  ios_add_rows(ios, 1);
         sprintf(name, "(%d)", i);
         ios_set_row_name(ios, i, name);
         ios_set_row_attr(ios, i, DEGREE_MARK, NULL);
         ios_set_row_bnds(ios, i, IOS_FX, 2.0, 2.0);
      }
      /* find an initial tour using the greedy heuristic */
      tour = ucalloc(1+nodes, sizeof(int));
      sum = greedy_tour(tsp, tour);
      print("Initial tour length: %g", sum);
      /* set global upper bound found by the primal heuristic */
      ios->found = 1;
      ios->best = sum;
      /* add columns which correspond to edges in the initial tour */
      for (k = 1; k <= nodes; k++)
      {  ios_add_cols(ios, 1);
         i = tour[k];
         j = tour[k < nodes ? k+1 : 1];
         if (i > j) temp = i, i = j, j = temp;
         sprintf(name, "%d-%d*", i, j);
         ios_set_col_name(ios, k, name);
         edge = dmp_get_atom(edge_pool);
         edge->i = i, edge->j = j;
         ios_set_col_attr(ios, k, EDGE_MARK, edge);
         ios_set_col_kind(ios, k, IOS_INT);
         ios_set_col_bnds(ios, k, IOS_DB, 0.0, 1.0);
         ios_set_obj_coef(ios, k, (double)tsp_distance(tsp, i, j));
         /* initialy there are no subtour elimination constraints, so
            each column has exactly two constraint coefficients */
         ind[1] = i, val[1] = 1.0;
         ind[2] = j, val[2] = 1.0;
         ios_set_mat_col(ios, k, 2, ind, val);
      }
      ufree(tour);
      return;
}

/*----------------------------------------------------------------------
-- obtain_col - obtain column of the constraint matrix.
--
-- This routine builds the complete list of constraint coefficients for
-- a column which corresponds to given edge (i,j) and is missing in the
-- current subproblem.
--
-- Row indices and numeric values of constraint coefficients are stored
-- in locations ind[1], ..., ind[len] and val[1], ..., val[len], resp.,
-- where len is the total number of non-zero constraint coefficients in
-- specified column. */

static int obtain_col(IOS *ios, int i, int j, int ind[], double val[])
{     SUBTOUR *subt;
      int k, len, nrows, *cut;
      insist(1 <= i && i < j && j <= nodes);
      /* allocate and clear working array */
      cut = ucalloc(1+nodes, sizeof(int));
      for (k = 1; k <= nodes; k++) cut[k] = 0;
      /* first nodes rows always correspond to degree constraints, in
         which any column has exactly two constraint coefficients */
      ind[1] = i, val[1] = 1.0;
      ind[2] = j, val[2] = 1.0;
      len = 2;
      /* walk through other rows that correspond to subtour elimination
         constraints */
      nrows = ios_get_num_rows(ios);
      for (k = nodes+1; k <= nrows; k++)
      {  /* restore corresponding cut (W, V \ W) */
         insist(ios_get_row_mark(ios, k) == SUBTOUR_MARK);
         for (subt = (SUBTOUR *)ios_get_row_link(ios, k); subt != NULL;
            subt = subt->next) cut[subt->i] = 1;
         /* if edge (i,j) gets under the cut, the corresponding column
            has non-zero constrint coefficients in this row */
         if (cut[i] && !cut[j] || !cut[i] && cut[j])
         {  len++;
            ind[len] = k;
            val[len] = 1.0;
         }
         /* clear the array cut */
         for (subt = (SUBTOUR *)ios_get_row_link(ios, k); subt != NULL;
            subt = subt->next) cut[subt->i] = 0;
      }
      /* free working array */
      ufree(cut);
      return len;
}

/*----------------------------------------------------------------------
-- gen_edge_col - generate column to include it in current subproblem.
--
-- This routine walks through all edges for which corresponding columns
-- are missing in the current subproblem and checks their reduced costs
-- to choose a column which is able to improve the objective (if infeas
-- is 0) or the sum of primal infeasibilities (if infeas is 1). If such
-- column exists, the routine includes it in the current subproblem. */

static void gen_edge_col(IOS *ios, int infeas)
{     EDGE *edge;
      int i, j, k, best_i, best_j, nrows, ncols, len, *ind;
      double rc, best_rc, *val;
      char name[255+1];
      /* determine the number of rows and columns which are presented
         in the current subroblem */
      nrows = ios_get_num_rows(ios);
      ncols = ios_get_num_cols(ios);
      /* allocate working arrays */
      ind = ucalloc(1+nrows, sizeof(int));
      val = ucalloc(1+nrows, sizeof(double));
      /* nothing is chosen so far */
      best_rc = 0.0, best_i = 0, best_j = 0;
      /* walk through all edges of the graph */
      for (i = 1; i <= nodes; i++)
      {  for (j = i+1; j <= nodes; j++)
         {  /* if column that corresponds to edge (i,j) is already in
               the current subproblem, skip it */
            for (k = 1; k <= ncols; k++)
            {  insist(ios_get_col_mark(ios, k) == EDGE_MARK);
               edge = ios_get_col_link(ios, k);
               if (edge->i == i && edge->j == j) goto skip;
            }
            /* obtain constraint coefficients for this column */
            len = obtain_col(ios, i, j, ind, val);
            /* compute the reduced cost of this column */
            rc = (infeas ? 0.0 : (double)tsp_distance(tsp, i, j));
            for (k = 1; k <= len; k++)
               rc += ios_get_row_pi(ios, ind[k]) * val[k];
            /* if the reduced cost indicates no improvement, skip the
               column */
            if (rc > -1e-5) goto skip;
            /* choose the only column using Dantzig's pricing */
            if (rc < best_rc) best_rc = rc, best_i = i, best_j = j;
skip:       ;
         }
      }
      /* if a column has been chosen, generate it, i.e. include it in
         the current subproblem */
      if (best_rc != 0.0)
      {  i = best_i, j = best_j;
         ios_add_cols(ios, 1);
         k = ios_get_num_cols(ios);
         sprintf(name, "%d-%d", i, j);
         ios_set_col_name(ios, k, name);
         edge = dmp_get_atom(edge_pool);
         edge->i = i, edge->j = j;
         ios_set_col_attr(ios, k, EDGE_MARK, edge);
         ios_set_col_kind(ios, k, IOS_INT);
         ios_set_col_bnds(ios, k, IOS_DB, 0.0, 1.0);
         ios_set_obj_coef(ios, k, (double)tsp_distance(tsp, i, j));
         len = obtain_col(ios, i, j, ind, val);
         ios_set_mat_col(ios, k, len, ind, val);
      }
      /* free working arrays */
      ufree(ind);
      ufree(val);
      return;
}

/*----------------------------------------------------------------------
-- max_flow - find max flow with the simplex method.
--
-- This routine finds max flow in a given undirected network by means
-- of the simplex method.
--
-- The undirected capacitated network is specified by the parameters
-- nn, ne, beg, end, and cap. The parameter nn is number of vertices
-- (nodes), nn > 0, and the parameter ne is number of edges, ne >= 0.
-- k-th edge is specified by the triple (beg[k], end[k], cap[k]) for
-- k = 1, ..., ne, where beg[k] and end[k] are numbers of the first and
-- the second nodes of the k-th edge (it should be beg[k] < end[k]),
-- cap[k] > 0 is a capacity of the k-th edge. Loops and multiple edges
-- are not allowed.
--
-- The parameter s is the number of a source node, and the parameter t
-- is the number of a sink node.
--
-- On exit the routine computes elementary flows along edges and stores
-- their values to locations x[1], ..., x[ne]. Positive value of x[k]
-- means that the elementary flow goes from the node beg[k] to the node
-- end[k], and negative value means that the flow goes in the opposite
-- direction. A value returned by the routine is the total maximum flow
-- carried through the network. */

static double max_flow(int nn, int ne, int beg[], int end[],
      double cap[], int s, int t, double x[])
{     LPX *lp;
      int i, k, nz, *rn, *cn;
      double flow, *aa;
      /* some sanity checks */
      insist(nn > 0);
      insist(ne >= 0);
      insist(1 <= s && s <= nn);
      insist(1 <= t && t <= nn);
      insist(s != t);
      for (k = 1; k <= ne; k++)
      {  insist(1 <= beg[k] && beg[k] < end[k] && end[k] <= nn);
         insist(cap[k] > 0.0);
      }
      /* create LP problem instance */
      lp = lpx_create_prob();
      /* create LP rows; i-th row is the conservation condition of the
         flow in the i-th node, i = 1, ..., nn */
      lpx_add_rows(lp, nn);
      for (i = 1; i <= nn; i++)
         lpx_set_row_bnds(lp, i, LPX_FX, 0.0, 0.0);
      /* create LP columns; k-th column is the elementary flow carried
         along the k-th edge, k = 1, ..., ne; the last column with the
         number ne+1 is the total flow through the network, which goes
         along a dummy edge from the sink to the source */
      lpx_add_cols(lp, ne + 1);
      for (k = 1; k <= ne; k++)
      {  if (cap[k] <= 1e-8)
            lpx_set_col_bnds(lp, k, LPX_FX, 0.0, 0.0);
         else
            lpx_set_col_bnds(lp, k, LPX_DB, - cap[k], + cap[k]);
      }
      lpx_set_col_bnds(lp, ne + 1, LPX_FR, 0.0, 0.0);
      /* build the constraint matrix; structurally this matrix is an
         incidence matrix of the network, so each its column (including
         the last column for the dummy edge) has exactly two non-zero
         entries */
      rn = ucalloc(1 + 2 * (ne + 1), sizeof(int));
      cn = ucalloc(1 + 2 * (ne + 1), sizeof(int));
      aa = ucalloc(1 + 2 * (ne + 1), sizeof(double));
      nz = 0;
      for (k = 1; k <= ne; k++)
      {  /* x[k] > 0 means the elementary flow through the k-th edge
            goes from the node beg[k] into the node end[k] */
         nz++, rn[nz] = beg[k], cn[nz] = k, aa[nz] = -1.0;
         nz++, rn[nz] = end[k], cn[nz] = k, aa[nz] = +1.0;
      }
      /* the total flow through the network goes from the sink to the
         source along the dummy edge */
      nz++, rn[nz] = t, cn[nz] = ne + 1, aa[nz] = -1.0;
      nz++, rn[nz] = s, cn[nz] = ne + 1, aa[nz] = +1.0;
      /* check number of non-zero entries */
      insist(nz == 2 * (ne + 1));
      /* load the constraint matrix into the LP problem object */
      lpx_load_mat3(lp, nz, rn, cn, aa);
      ufree(rn);
      ufree(cn);
      ufree(aa);
      /* the objective function is the total flow through the network,
         which should be maximized */
      lpx_set_obj_dir(lp, LPX_MAX);
      lpx_set_col_coef(lp, ne + 1, 1.0);
      /* solve the LP problem */
      lpx_set_int_parm(lp, LPX_K_MSGLEV, 2);
      lpx_set_real_parm(lp, LPX_K_OUTDLY, 5.0);
      lpx_adv_basis(lp);
      insist(lpx_simplex(lp) == LPX_E_OK);
      insist(lpx_get_status(lp) == LPX_OPT);
      /* obtain the maximum flow through the network and the elementary
         flows through edges that correspond to the optimal solution */
      flow = lpx_get_obj_val(lp);
      for (k = 1; k <= ne; k++)
         lpx_get_col_info(lp, k, NULL, &x[k], NULL);
      /* delete LP problem instance */
      lpx_delete_prob(lp);
      /* return to the calling program */
      return flow;
}

/*----------------------------------------------------------------------
-- min_st_cut - find min (s,t)-cut for known max flow.
--
-- This routine finds min (s,t)-cut, which corresponds to a known max
-- flow in a given undirected network.
--
-- This routine should be called after the routine max_flow (see above)
-- with the same parameters, where the array x contains elementary flows
-- through edges of the network, which correspond to maximal flow found
-- by the routine max_flow.
--
-- The routine min_st_cut splits the set of nodes V of the network into
-- two non-empty subsets V(s) and V(t) = V \ V(s), where the source node
-- s belongs to V(s), the sink node t belongs to V(t), and all edges,
-- one node of which belongs to V(s) and other one belongs to V(t), are
-- saturated (i.e. x[k] = +cap[k] or x[k] = -cap[k]).
--
-- On exit the routine stores flags of the nodes v[i], i = 1, ..., nn,
-- to the locations cut[i], where cut[i] = 1 means v[i] belongs to V(s)
-- and cut[i] = 0 means v[i] belongs to V(t) = V \ V(s). The routine
-- also returns a value of the corresponding minimal (s,t)-cut, which is
-- the sum of capacities of all edges between V(s) and V(t). (Due to the
-- theorem of Ford and Fulkerson the value of the minimal cut should be
-- the same as the value of maximal flow.)
--
-- In order to determine the set V(s) the routine just finds all nodes,
-- which can be reached from the source node s via non-saturated edges.
-- The set V(t) is determined as the complement V \ V(s). */

static double min_st_cut(int nn, int ne, int beg[], int end[],
      double cap[], int s, int t, double x[], int cut[])
{     int i, j, k, p, q, *head1, *next1, *head2, *next2, *list;
      double temp;
      /* head1[i] points to the first edge with beg[k] = i
         next1[k] points to the next edge with the same beg[k]
         head2[i] points to the first edge with end[k] = i
         next2[k] points to the next edge with the same end[k] */
      head1 = ucalloc(1+nn, sizeof(int));
      head2 = ucalloc(1+nn, sizeof(int));
      next1 = ucalloc(1+ne, sizeof(int));
      next2 = ucalloc(1+ne, sizeof(int));
      for (i = 1; i <= nn; i++) head1[i] = head2[i] = 0;
      for (k = 1; k <= ne; k++)
      {  i = beg[k], next1[k] = head1[i], head1[i] = k;
         j = end[k], next2[k] = head2[j], head2[j] = k;
      }
      /* on constructing the set V(s) list[1], ..., list[p-1] contain
         nodes, which can be reached from the source node and have been
         visited, and list[p], ..., list[q] contain nodes, which can be
         reached from the source node but havn't been visited yet */
      list = ucalloc(1+nn, sizeof(int));
      for (i = 1; i <= nn; i++) cut[i] = 0;
      p = q = 1, list[1] = s, cut[s] = 1;
      while (p <= q)
      {  /* pick the next node, which is reachable from the source node
            and has not visited yet, and visit it */
         i = list[p++];
         /* walk through edges with beg[k] = i */
         for (k = head1[i]; k != 0; k = next1[k])
         {  j = end[k];
            insist(beg[k] == i);
            /* from v[i] we can reach v[j], if the elementary flow from
               v[i] to v[j] is non-saturated */
            if (cut[j] == 0 && x[k] <= + (cap[k] - 1e-7))
               list[++q] = j, cut[j] = 1;
         }
         /* walk through edges with end[k] = i */
         for (k = head2[i]; k != 0; k = next2[k])
         {  j = beg[k];
            insist(end[k] == i);
            /* from v[i] we can reach v[j], if the elementary flow from
               v[i] to v[j] is non-saturated */
            if (cut[j] == 0 && x[k] >= - (cap[k] - 1e-7))
               list[++q] = j, cut[j] = 1;
         }
      }
      /* the sink cannot belong to V(s) */
      insist(cut[t] == 0);
      /* free working arrays */
      ufree(head1);
      ufree(head2);
      ufree(next1);
      ufree(next2);
      ufree(list);
      /* compute value of the corresponding minimal (s,t)-cut */
      temp = 0.0;
      for (k = 1; k <= ne; k++)
      {  i = beg[k], j = end[k];
         if (cut[i] && !cut[j] || !cut[i] && cut[j]) temp += cap[k];
      }
      /* return to the calling program */
      return temp;
}

/*----------------------------------------------------------------------
-- stoer_wagner - find minimal cut with Stoer and Wagner algorithm.
--
-- This routine finds min cut in a given undirected network by means of
-- Stoer-Wagner algorithm.
--
-- The undirected capacitated network is specified by the parameters
-- nn, ne, beg, end, and cap. The parameter nn is number of vertices
-- (nodes), nn > 0, and the parameter ne is number of edges, ne >= 0.
-- k-th edge is specified by the triple (beg[k], end[k], cap[k]) for
-- k = 1, ..., ne, where beg[k] and end[k] are numbers of the first and
-- the second nodes of the k-th edge (it should be beg[k] < end[k]),
-- cap[k] > 0 is a capacity of the k-th edge. Loops and multiple edges
-- are not allowed.
--
-- Let V be the set of nodes of the network and let W be an arbitrary
-- non-empty proper subset of V. A cut associated with the subset W is
-- a subset of all the edges, one node of which belongs to W and other
-- node belongs to V \ W. The capacity of a cut (W, V \ W) is the sum
-- of the capacities of all edges, which belong to the cut. Minimal cut
-- is a cut, whose capacity is minimal.
--
-- On exit the routine stores flags of the nodes v[i], i = 1, ..., nn,
-- to the locations cut[i], where cut[i] = 1 means v[i] belongs to W
-- and cut[i] = 0 means v[i] belongs to V \ W, where W corresponds to
-- a minimal cut. The routine returns the capacity of the minimal cut.
--
-- The basic idea of Stoer-Wagner algorithm is the following. Let G be
-- a capacitated network and G(s,t) be a network, in which the nodes s
-- and t are merged into one new node, loops are deleted, but multuple
-- edges are retained. It is obvious that a minimum cut in G is the
-- minimum of two quantities: the minimum cut in G(s,t) and a minimum
-- cut that separates s and t. This allows to find a minimum cut in the
-- original network solving at most nn max flow problems.
--
-- M. Stoer, F. Wagner. A Simple Min Cut Algorithm. Algorithms, ESA'94
-- LNCS 855 (1994), pp. 141-47.
--
-- J. Cheriyan, R. Ravi. Approximation Algorithms for Network Problems.
-- Univ. of Waterloo (1998), p. 147. */

static double stoer_wagner(int nn, int ne, int beg[], int end[],
      double cap[], int cut[])
{     int i, j, k, *head1, *next1, *head2, *next2, I, J, K, S, T,
         DEG, NV, NE, *BEG, *END, *HEAD, *NEXT, *NUMB, *ADJ, *CUT;
      double min_cut, flow, temp, *X, *CAP, *SUM;
      /* some sanity checks */
      insist(nn > 0);
      insist(ne >= 0);
      for (k = 1; k <= ne; k++)
      {  insist(1 <= beg[k] && beg[k] < end[k] && end[k] <= nn);
         insist(cap[k] > 0.0);
      }
      /* head1[i] points to the first edge with beg[k] = i
         next1[k] points to the next edge with the same beg[k]
         head2[i] points to the first edge with end[k] = i
         next2[k] points to the next edge with the same end[k] */
      head1 = ucalloc(1+nn, sizeof(int));
      head2 = ucalloc(1+nn, sizeof(int));
      next1 = ucalloc(1+ne, sizeof(int));
      next2 = ucalloc(1+ne, sizeof(int));
      for (i = 1; i <= nn; i++) head1[i] = head2[i] = 0;
      for (k = 1; k <= ne; k++)
      {  i = beg[k], next1[k] = head1[i], head1[i] = k;
         j = end[k], next2[k] = head2[j], head2[j] = k;
      }
      /* an auxiliary network used in the algorithm is resulted from
         the original network by merging some nodes into one supernode;
         all variables and arrays related to this auxiliary network are
         denoted in caps */
      /* HEAD[I] points to the first node of the original network that
         belongs to the I-th supernode
         NEXT[i] points to the next node of the original network that
         belongs to the same supernode as the i-th node
         NUMB[i] is a supernode, which the i-th node belongs to */
      /* initially the auxiliary network is equivalent to the original
         network, i.e. each supernode consists of one node */
      NV = nn;
      HEAD = ucalloc(1+nn, sizeof(int));
      NEXT = ucalloc(1+nn, sizeof(int));
      NUMB = ucalloc(1+nn, sizeof(int));
      for (i = 1; i <= nn; i++) HEAD[i] = i, NEXT[i] = 0, NUMB[i] = i;
      /* number of edges in the auxiliary network is never greater than
         in the original one */
      BEG = ucalloc(1+ne, sizeof(int));
      END = ucalloc(1+ne, sizeof(int));
      CAP = ucalloc(1+ne, sizeof(double));
      X = ucalloc(1+ne, sizeof(double));
      /* allocate some auxiliary arrays */
      ADJ = ucalloc(1+nn, sizeof(int));
      SUM = ucalloc(1+nn, sizeof(double));
      CUT = ucalloc(1+nn, sizeof(int));
      /* currently no min cut is known so far */
      min_cut = DBL_MAX;
      /* main loop starts here */
      while (NV > 1)
      {  /* build the set of edges of the auxiliary network */
         NE = 0;
         /* multiple edges are not allowed in the max flow algorithm,
            so we can replace each multiple edge, which is the result
            of merging nodes into supernodes, by a single edge, whose
            capacity is the sum of capacities of particular edges;
            these summary capacities will be stored in the array SUM */
         for (I = 1; I <= NV; I++) SUM[I] = 0.0;
         for (I = 1; I <= NV; I++)
         {  /* DEG is number of single edges, which connects the I-th
               supernode and some J-th supernode, where I < J */
            DEG = 0;
            /* walk through nodes that belong to the I-th supernode */
            for (i = HEAD[I]; i != 0; i = NEXT[i])
            {  /* the i-th node belongs to the I-th supernode */
               /* walk through edges with beg[k] = i */
               for (k = head1[i]; k != 0; k = next1[k])
               {  j = end[k];
                  /* j-th node belongs to the J-th supernode */
                  J = NUMB[j];
                  /* ignore loops and edges with I > J */
                  if (I >= J) continue;
                  /* add an edge, which connects the I-th and the J-th
                     supernodes (if not added yet) */
                  if (SUM[J] == 0.0) ADJ[++DEG] = J;
                  /* sum up the capacity of the original edge */
                  insist(cap[k] > 0.0);
                  SUM[J] += cap[k];
               }
               /* walk through edges with end[k] = i */
               for (k = head2[i]; k != 0; k = next2[k])
               {  j = beg[k];
                  /* j-th node belongs to the J-th supernode */
                  J = NUMB[j];
                  /* ignore loops and edges with I > J */
                  if (I >= J) continue;
                  /* add an edge, which connects the I-th and the J-th
                     supernodes (if not added yet) */
                  if (SUM[J] == 0.0) ADJ[++DEG] = J;
                  /* sum up the capacity of the original edge */
                  insist(cap[k] > 0.0);
                  SUM[J] += cap[k];
               }
            }
            /* add single edges, which connect the I-th supernode to
               other supernodes, to the auxiliary network, and restore
               the array SUM for subsequent use */
            for (K = 1; K <= DEG; K++)
            {  NE++;
               insist(NE <= ne);
               J = ADJ[K];
               BEG[NE] = I, END[NE] = J, CAP[NE] = SUM[J];
               SUM[J] = 0.0;
            }
         }
         /* choose two arbitrary supernodes of the auxiliary network,
            one of which is the source and other is the sink */
         S = 1, T = NV;
         /* determine max flow from S to T */
         flow = max_flow(NV, NE, BEG, END, CAP, S, T, X);
         /* if the min cut, which separates the supernodes S and T, is
            less than the currently known, remember it */
         if (min_cut > flow)
         {  min_cut = flow;
            /* determine the min cut in the auxiliary network */
            temp = min_st_cut(NV, NE, BEG, END, CAP, S, T, X, CUT);
            /* check that Ford and Fulkerson are never mistaken :+) */
            insist(fabs(flow - temp) <= 1e-6 * (1.0 + fabs(flow)));
            /* determine the min cut in the original network */
            for (i = 1; i <= nn; i++) cut[i] = CUT[NUMB[i]];
            /* if the min cut is close to zero (that obviously means
               the network has unconnected components), the search can
               be prematurely terminated */
            if (min_cut <= 1e-8) break;
         }
         /* now merge all nodes of the original network, which belong
            to the supernodes S and T, into one new supernode; this is
            attained by carrying all nodes from T to S (for the sake of
            convenience T should be the last supernode) */
         insist(T == NV);
         /* assign new references to nodes from T */
         for (i = HEAD[T]; i != 0; i = NEXT[i]) NUMB[i] = S;
         /* find the last entry in the node list of S */
         i = HEAD[S];
         insist(i != 0);
         while (NEXT[i] != 0) i = NEXT[i];
         /* and attach to it the node list of T */
         NEXT[i] = HEAD[T];
         /* decrease number of nodes in the auxiliary network */
         NV--;
      }
      /* free working arrays */
      ufree(HEAD);
      ufree(NEXT);
      ufree(NUMB);
      ufree(BEG);
      ufree(END);
      ufree(CAP);
      ufree(X);
      ufree(ADJ);
      ufree(SUM);
      ufree(CUT);
      ufree(head1);
      ufree(head2);
      ufree(next1);
      ufree(next2);
      /* return to the calling program */
      return min_cut;
}

/*----------------------------------------------------------------------
-- gen_subt_row - generate subtour elimination constraint.
--
-- This routine is called from the application procedure in order to
-- generate a violated subtour elimination constraint.
--
-- Constraints of this class has the form:
--
--    sum x[i,j] >= 2, i in W, j in V \ W,
--
-- for all W, where W is a proper nonempty subset of V, V is the set of
-- nodes of the given graph.
--
-- In order to find a violated constraint of this class this routine
-- finds a min cut in a capacitated network, which has the same sets of
-- nodes and edges as the original graph, and where capacities of edges
-- are values of variables x[i,j] in a basic solution of the current
-- subproblem. */

static void gen_subt_row(IOS *ios)
{     SUBTOUR *head, *subt;
      EDGE *edge;
      int i, j, k, ncols, nn, ne, nz, len, *beg, *end, *cut, *ind;
      double x, min_cut, *cap, *val;
      /* the network has the same set of nodes as the original graph */
      nn = nodes;
      /* if some variable x[i,j] is zero in basic solution of the
         current subproblem, the capacity of corresponding edge in the
         network is zero, therefore such edge may not be included in
         the network; note that if a variable is not presented in the
         current subproblem, its value is zero by definition, so only
         variables presented in the subproblem should be considered */
      ncols = ios_get_num_cols(ios);
      /* count number of edges with non-zero capacity */
      ne = 0;
      for (k = 1; k <= ncols; k++)
      {  ios_get_col_soln(ios, k, &x, NULL);
         if (x >= 1e-5) ne++;
      }
      /* create the capacitated network */
      beg = ucalloc(1+ne, sizeof(int));
      end = ucalloc(1+ne, sizeof(int));
      cap = ucalloc(1+ne, sizeof(double));
      nz = 0;
      for (k = 1; k <= ncols; k++)
      {  ios_get_col_soln(ios, k, &x, NULL);
         if (x >= 1e-5)
         {  insist(ios_get_col_mark(ios, k) == EDGE_MARK);
            edge = ios_get_col_link(ios, k);
            i = edge->i, j = edge->j;
            insist(1 <= i && i < j && j <= nn);
            nz++;
            beg[nz] = i, end[nz] = j, cap[nz] = x;
         }
      }
      insist(nz == ne);
      /* find minimal cut in the capacitated network */
      cut = ucalloc(1+nn, sizeof(int));
      min_cut = stoer_wagner(nn, ne, beg, end, cap, cut);
      /* if the capacity of min cut is less than 2, the corresponding
         subtour elimination constraint is violated */
      if (min_cut <= 2.0 - 1e-5)
      {  /* build linked list of nodes which belongs to the set W to
            save it in the descriptor of new row */
         head = NULL;
         for (i = 1; i <= nn; i++)
         {  if (cut[i])
            {  subt = dmp_get_atom(subtour_pool);
               subt->i = i;
               subt->next = head;
               head = subt;
            }
         }
         /* determine constraint coefficients for the constraint to be
            generated; each edge that gets into the cut gives non-zero
            constraint coefficient at corresponding column */
         ind = ucalloc(1+ncols, sizeof(int));
         val = ucalloc(1+ncols, sizeof(double));
         len = 0;
         for (k = 1; k <= ncols; k++)
         {  insist(ios_get_col_mark(ios, k) == EDGE_MARK);
            edge = ios_get_col_link(ios, k);
            i = edge->i, j = edge->j;
            insist(1 <= i && i < j && j <= nn);
            if (cut[i] && !cut[j] || !cut[i] && cut[j])
            {  len++;
               ind[len] = k;
               val[len] = 1.0;
            }
         }
         /* include the new row in the current subproblem */
         ios_add_rows(ios, 1);
         k = ios_get_num_rows(ios);
         ios_set_row_attr(ios, k, SUBTOUR_MARK, head);
         ios_set_row_bnds(ios, k, IOS_LO, 2.0, 0.0);
         ios_set_mat_row(ios, k, len, ind, val);
         ufree(ind);
         ufree(val);
      }
      /* free working arrays */
      ufree(cut);
      ufree(beg);
      ufree(end);
      ufree(cap);
      return;
}

/*----------------------------------------------------------------------
-- branching - choose column to branch on.
--
-- This routine chooses a column (structural variable) to branch on.
--
-- That column is chosen whose primal value in basic solution of the
-- current LP relaxation is fractional and closest to 0.5. */

static void branching(IOS *ios)
{     int nv, j, j_min = 0;
      double x, d, d_min = DBL_MAX;
      /* walk through all columns */
      nv = ios_get_num_cols(ios);
      for (j = 1; j <= nv; j++)
      {  if (ios_is_col_frac(ios, j))
         {  /* primal value of j-th column is fractional */
            ios_get_col_soln(ios, j, &x, NULL);
            d = fabs(x - 0.5);
            if (d_min > d) j_min = j, d_min = d;
         }
      }
      /* branch on column j_min and solve the up branch next */
      ios_branch_on(ios, j_min, +1);
      return;
}

/*----------------------------------------------------------------------
-- select_node - select subproblem to continue the search.
--
-- This routine selects a subproblem from the active list to continue
-- the search.
--
-- That subproblem is selected which was added to the active list most
-- recently. This corresponds to the depth first search. */

static void select_node(IOS *ios)
{     int p;
      /* determine the reference number of the most recent subproblem;
         it is the last subproblem in the active list */
      p = ios_get_prev_node(ios, 0);
      /* select the subproblem p to continue the search */
      ios_select_node(ios, p);
      return;
}

/*----------------------------------------------------------------------
-- appl_proc - event-driven application procedure.
--
-- This routine is an event-driven application procedure called from
-- the IOS driver at certain points of the optimization process. */

static void appl_proc(IOS *ios, void *info)
{     insist(info == info);
      switch (ios->event)
      {  case IOS_V_INIT:
            /* create root subproblem */
            initialize(ios);
            /* enable column generation */
            ios->col_gen = 1;
            /* enable row generation */
            ios->row_gen = 1;
            /* the objective function is integral */
            ios->int_obj = 1;
            /* set message level */
            ios->msg_lev = (trace ? 3 : 2);
            break;
         case IOS_V_GENCOL:
            /* column generation is required */
            if (ios_p_status(ios) == IOS_FEAS)
            {  /* LP relaxation of the current subproblem has optimal
                  solution */
               gen_edge_col(ios, 0);
            }
            else
            {  /* LP relaxation of the current subproblem has no primal
                  feasible solution */
               gen_edge_col(ios, 1);
            }
            break;
         case IOS_V_GENROW:
            /* row generation required */
            gen_subt_row(ios);
            break;
         case IOS_V_BINGO:
            /* better integer feasible solution found */
            if (out_sol != NULL)
            {  /* extract LP relaxation of the current subproblem and
                  write its (integer) optimal solution to output file
                  in plain text format */
               LPX *lp = ios_extract_lp(ios);
               insist(lp != NULL);
               lpx_set_prob_name(lp, in_file);
               insist(lpx_warm_up(lp) == LPX_E_OK);
               insist(lpx_get_status(lp) == LPX_OPT);
               lpx_print_sol(lp, out_sol);
               lpx_delete_prob(lp);
            }
            break;
         case IOS_V_BRANCH:
            /* choose column to branch on */
            branching(ios);
            break;
         case IOS_V_SELECT:
            /* select active subproblem to continue the search */
            select_node(ios);
            break;
         case IOS_V_DELROW:
            /* some row is being deleted; free memory allocated to its
               application extension */
            switch (ios_get_row_mark(ios, 0))
            {  case DEGREE_MARK:
                  /* it is degree constraint */
                  break;
               case SUBTOUR_MARK:
                  /* it is subtour elimination constraint */
                  {  SUBTOUR *head, *subt;
                     head = ios_get_row_link(ios, 0);
                     while (head != NULL)
                     {  subt = head;
                        head = subt->next;
                        dmp_free_atom(subtour_pool, subt);
                     }
                  }
                  break;
               default:
                  insist(ios != ios);
            }
            break;
         case IOS_V_DELCOL:
            /* some column is being deleted; free memory allocated to
               its application extension */
            switch (ios_get_col_mark(ios, 0))
            {  case EDGE_MARK:
                  /* it is edge variable */
                  {  EDGE *edge;
                     edge = ios_get_col_link(ios, 0);
                     dmp_free_atom(edge_pool, edge);
                  }
                  break;
               default:
                  insist(ios != ios);
            }
            break;
         default:
            /* ignore other events */
            break;
      }
      /* return to the IOS driver */
      return;
}

/*----------------------------------------------------------------------
-- main - main program.
--
-- This main program parses command-line parameters and calls the IOS
-- driver to solve the Symmetric Traveling Salesman Problem (TSP). */

int main(int argc, char *argv[])
{     int k;
      /* parse command-line parameters */
#     define p(str) (strcmp(argv[k], str) == 0)
      for (k = 1; k < argc; k++)
      {  if (p("--help") || p("-h"))
         {  print("Usage: %s [options...] tsp-file", argv[0]);
            print("");
            print("Options:");
            print("   -o filename, --output filename");
            print("                     write solution to filename in p"
               "lain text format");
            print("   -t, --trace       produce detailed output");
            exit(EXIT_SUCCESS);
         }
         else if (p("--output") || p("-o"))
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
         else if (p("--trace") || p("-t"))
            trace = 1;
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
#     undef p
      /* remove output file specified in command-line */
      if (out_sol != NULL) remove(out_sol);
      /* read problem instance from input file */
      if (in_file == NULL)
      {  print("No input file specified; try %s --help", argv[0]);
         exit(EXIT_FAILURE);
      }
      tsp = tsp_read_data(in_file);
      if (tsp == NULL)
      {  print("TSP file processing error");
         exit(EXIT_FAILURE);
      }
      nodes = tsp->dimension;
      /* create dynamic memory pools */
      edge_pool = dmp_create_pool(sizeof(EDGE));
      subtour_pool = dmp_create_pool(sizeof(SUBTOUR));
      /* call the IOS driver */
      ios_driver(appl_proc, NULL);
      /* delete dynamic memory pools */
      insist(edge_pool->count == 0);
      dmp_delete_pool(edge_pool);
      insist(subtour_pool->count == 0);
      dmp_delete_pool(subtour_pool);
      /* delete problem instance */
      tsp_free_data(tsp);
      /* check that no memory blocks have been lost */
      insist(lib_env_ptr()->mem_total == 0);
      insist(lib_env_ptr()->mem_count == 0);
      /* return to the control program */
      return 0;
}

/* eof */

/*----------------------------------------------------------------------
 *
 * Copyright (C) 2001-2004, Nicolo' Giorgetti. All rights reserved. 
 * E-mail: <giorgetti@dii.unisi.it>.
 *
 * This file is part of GLPK (GNU Linear Programming Kit).
 *
 * GLPK is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * GLPK is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GLPK; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 *
 *-----------------------------------------------------------------------*/
 
static char *header= "glpkmex: MEX interface for the GLPK library\n";
static char *version="Version: 0.6.3\n";
static char *copyright="(C) 2001-2004, Nicolo' Giorgetti.\n";
static char *syntax="\nUsage: [xmin,fmin,status,extra]=glpkmex(sense,c,a,b,ctype,lb,ub,vartype,param,lpsolver,save)\n"; 
 
#define NIntP 17
#define NRealP 10


int lpxIntParam[NIntP]= {
   1,
   1,
   0,
   1,
   0,
   -1,
   0,
   200,
   1,
   2,
   0,
   1,
   0,
   0,
   2, 
   2,
   1
};

int IParam[NIntP]={
   LPX_K_MSGLEV,
   LPX_K_SCALE,
   LPX_K_DUAL,
   LPX_K_PRICE,
   LPX_K_ROUND,
   LPX_K_ITLIM,
   LPX_K_ITCNT,
   LPX_K_OUTFRQ,
   LPX_K_MPSINFO,
   LPX_K_MPSOBJ,
   LPX_K_MPSORIG,
   LPX_K_MPSWIDE,
   LPX_K_MPSFREE,
   LPX_K_MPSSKIP,
   LPX_K_BRANCH,
   LPX_K_BTRACK,
   LPX_K_PRESOL
};


double lpxRealParam[NRealP]={
   0.07,
   1e-7,
   1e-7,
   1e-9,
   -DBL_MAX,
   DBL_MAX,
   -1.0,
   0.0,
   1e-6,
   1e-7
};

int RParam[NRealP]={
   LPX_K_RELAX,
   LPX_K_TOLBND,
   LPX_K_TOLDJ,
   LPX_K_TOLPIV,
   LPX_K_OBJLL,
   LPX_K_OBJUL,
   LPX_K_TMLIM,
   LPX_K_OUTDLY,
   LPX_K_TOLINT,
   LPX_K_TOLOBJ
};

jmp_buf mark;              /* Address for long jump to jump to */
int     fperr;             /* Global error number */


// ---------------------------------------------------------------------
// Copyright (C) 2003 Andrew Makhorin <mao@mai2.rcnet.ru>, Department
// for Applied Informatics, Moscow Aviation Institute, Moscow, Russia.
// All rights reserved.
//
// Author: Ivo van Baren, Operations Research specialist,
//         i.van.baren@freeler.nl
//
// This file is a part of GLPK (GNU Linear Programming Kit).
//
// GLPK is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// GLPK is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GLPK; see the file COPYING. If not, write to the Free
// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
// --------------------------------------------------------------------


Unit Sample;

Interface

Uses
 Dialogs, SysUtils, Glpk4;

Procedure GLPKSolve;

Implementation

Procedure GLPKSolve;

  Var
   lp          : PLPX;
   rn,
   cn          : PIntArray;
   a           : PFloatArray;
   IntDummy    : Integer;
   Z,
   FloatDummy,
   x1, x2, x3  : Double;
   probname    : String;

 Begin
  GetMem(rn,10*SizeOf(Integer));
  GetMem(cn,10*SizeOf(Integer));
  GetMem(a,10*SizeOf(Double));

  lp := lpx_create_prob;

  lpx_set_prob_name(lp, PutZTString('Sample'));

  {* To demonstrate the use of strings;
     the := operator copies
     a PChar-string type to a String type  *}
  {*
  probname := lpx_get_prob_name(lp);
   ShowMessage(probname);
  *}

  lpx_add_rows(lp, 3);
  lpx_set_row_name(lp, 1, PutZTString('p'));
  lpx_set_row_bnds(lp, 1, LPX_UP, 0.0, 100.0);
  lpx_set_row_name(lp, 2, PutZTString('q'));
  lpx_set_row_bnds(lp, 2, LPX_UP, 0.0, 600.0);
  lpx_set_row_name(lp, 3, PutZTString('r'));
  lpx_set_row_bnds(lp, 3, LPX_UP, 0.0, 300.0);

  lpx_add_cols(lp, 3);
  lpx_set_col_name(lp, 1, PutZTString('x1'));
  lpx_set_col_bnds(lp, 1, LPX_LO, 0.0, 0.0);
  lpx_set_col_name(lp, 2, PutZTString('x2'));
  lpx_set_col_bnds(lp, 2, LPX_LO, 0.0, 0.0);
  lpx_set_col_name(lp, 3, PutZTString('x3'));
  lpx_set_col_bnds(lp, 3, LPX_LO, 0.0, 0.0);

  rn[1] := 1; cn[1] := 1; a[1] :=  1.0;
  rn[2] := 1; cn[2] := 2; a[2] :=  1.0;
  rn[3] := 1; cn[3] := 3; a[3] :=  1.0;
  rn[4] := 2; cn[4] := 1; a[4] := 10.0;
  rn[5] := 3; cn[5] := 1; a[5] :=  2.0;
  rn[6] := 2; cn[6] := 2; a[6] :=  4.0;
  rn[7] := 3; cn[7] := 2; a[7] :=  2.0;
  rn[8] := 2; cn[8] := 3; a[8] :=  5.0;
  rn[9] := 3; cn[9] := 3; a[9] :=  6.0;
  lpx_load_mat3(lp, 9, rn, cn, a);

  lpx_set_obj_dir(lp, LPX_MAX);

  lpx_set_col_coef(lp, 1, 10.0);
  lpx_set_col_coef(lp, 2, 6.0);
  lpx_set_col_coef(lp, 3, 4.0);

  lpx_simplex(lp);

  Z := lpx_get_obj_val(lp);

  lpx_get_col_info(lp, 1, IntDummy, x1, FloatDummy);
  lpx_get_col_info(lp, 2, IntDummy, x2, FloatDummy);
  lpx_get_col_info(lp, 3, IntDummy, x3, FloatDummy);

  ShowMessage(Format('Solution: Z = %g; x1 = %g; x2 = %g; x3 = %g', [Z, x1, x2, x3]));

  lpx_print_sol(lp,PutZTString('Sample_Solution.txt'));

  lpx_delete_prob(lp);

 End;




end.

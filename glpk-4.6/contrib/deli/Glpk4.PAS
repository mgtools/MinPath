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

Unit GLPK4;

Interface

Const
   LPX_LP          =  100;   {* linear programming (LP) *}
   LPX_MIP         =  101;   {* mixed integer programming (MIP) *}
   LPX_FR          =  110;   {* free variable:  -inf <  x[k] < +inf  *}
   LPX_LO          =  111;   {* lower bound:    l[k] <= x[k] < +inf  *}
   LPX_UP          =  112;   {* upper bound:    -inf <  x[k] <= u[k] *}
   LPX_DB          =  113;   {* double bound:   l[k] <= x[k] <= u[k] *}
   LPX_FX          =  114;   {* fixed variable: l[k]  = x[k]  = u[k] *}
   LPX_MIN         =  120;   {* minimization *}
   LPX_MAX         =  121;   {* maximization *}
   LPX_B_UNDEF     =  130;   {* current basis is undefined *}
   LPX_B_VALID     =  131;   {* current basis is valid *}
   LPX_P_UNDEF     =  132;   {* primal status is undefined *}
   LPX_P_FEAS      =  133;   {* solution is primal feasible *}
   LPX_P_INFEAS    =  134;   {* solution is primal infeasible *}
   LPX_P_NOFEAS    =  135;   {* no primal feasible solution exists *}
   LPX_D_UNDEF     =  136;   {* dual status is undefined *}
   LPX_D_FEAS      =  137;   {* solution is dual feasible *}
   LPX_D_INFEAS    =  138;   {* solution is dual infeasible *}
   LPX_D_NOFEAS    =  139;   {* no dual feasible solution exists *}
   LPX_BS          =  140;   {* basic variable *}
   LPX_NL          =  141;   {* non-basic variable on lower bound *}
   LPX_NU          =  142;   {* non-basic variable on upper bound *}
   LPX_NF          =  143;   {* non-basic free variable *}
   LPX_NS          =  144;   {* non-basic fixed variable *}
   LPX_T_UNDEF     =  150;   {* solution is undefined *}
   LPX_T_OPT       =  151;   {* solution is optimal *}
   LPX_CV          =  160;   {* continuous variable *}
   LPX_IV          =  161;   {* integer variable *}
   LPX_I_UNDEF     =  170;   {* integer status is undefined *}
   LPX_I_OPT       =  171;   {* solution is integer optimal *}
   LPX_I_FEAS      =  172;   {* solution is integer feasible *}
   LPX_I_NOFEAS    =  173;   {* no integer solution exists *}
{* status codes reported by the routine lpx_get_status: *}
   LPX_OPT         =  180;   {* optimal *}
   LPX_FEAS        =  181;   {* feasible *}
   LPX_INFEAS      =  182;   {* infeasible *}
   LPX_NOFEAS      =  183;   {* no feasible *}
   LPX_UNBND       =  184;   {* unbounded *}
   LPX_UNDEF       =  185;   {* undefined *}
{* exit codes returned by the simplex-based solver routines: *}
   LPX_E_OK        =  200;   {* success *}
   LPX_E_EMPTY     =  201;   {* empty problem *}
   LPX_E_BADB      =  202;   {* invalid initial basis *}
   LPX_E_INFEAS    =  203;   {* infeasible initial solution *}
   LPX_E_FAULT     =  204;   {* unable to start the search *}
   LPX_E_OBJLL     =  205;   {* objective lower limit reached *}
   LPX_E_OBJUL     =  206;   {* objective upper limit reached *}
   LPX_E_ITLIM     =  207;   {* iterations limit exhausted *}
   LPX_E_TMLIM     =  208;   {* time limit exhausted *}
   LPX_E_NOFEAS    =  209;   {* no feasible solution *}
   LPX_E_INSTAB    =  210;   {* numerical instability *}
   LPX_E_SING      =  211;   {* problems with basis matrix *}
   LPX_E_NOCONV    =  212;   {* no convergence (interior) *}
   LPX_E_NOPFS     =  213;   {* no primal feas. sol. (LP presolver) *}
   LPX_E_NODFS     =  214;   {* no dual feas. sol.   (LP presolver) *}
{* control parameter identifiers: *}
   LPX_K_MSGLEV    =  300;   {* lp->msg_lev *}
   LPX_K_SCALE     =  301;   {* lp->scale *}
   LPX_K_DUAL      =  302;   {* lp->dual *}
   LPX_K_PRICE     =  303;   {* lp->price *}
   LPX_K_RELAX     =  304;   {* lp->relax *}
   LPX_K_TOLBND    =  305;   {* lp->tol_bnd *}
   LPX_K_TOLDJ     =  306;   {* lp->tol_dj *}
   LPX_K_TOLPIV    =  307;   {* lp->tol_piv *}
   LPX_K_ROUND     =  308;   {* lp->round *}
   LPX_K_OBJLL     =  309;   {* lp->obj_ll *}
   LPX_K_OBJUL     =  310;   {* lp->obj_ul *}
   LPX_K_ITLIM     =  311;   {* lp->it_lim *}
   LPX_K_ITCNT     =  312;   {* lp->it_cnt *}
   LPX_K_TMLIM     =  313;   {* lp->tm_lim *}
   LPX_K_OUTFRQ    =  314;   {* lp->out_frq *}
   LPX_K_OUTDLY    =  315;   {* lp->out_dly *}
   LPX_K_BRANCH    =  316;   {* lp->branch *}
   LPX_K_BTRACK    =  317;   {* lp->btrack *}
   LPX_K_TOLINT    =  318;   {* lp->tol_int *}
   LPX_K_TOLOBJ    =  319;   {* lp->tol_obj *}
   LPX_K_MPSINFO   =  320;   {* lp->mps_info *}
   LPX_K_MPSOBJ    =  321;   {* lp->mps_obj *}
   LPX_K_MPSORIG   =  322;   {* lp->mps_orig *}
   LPX_K_MPSWIDE   =  323;   {* lp->mps_wide *}
   LPX_K_MPSFREE   =  324;   {* lp->mps_free *}
   LPX_K_MPSSKIP   =  325;   {* lp->mps_skip *}
   LPX_K_LPTORIG   =  326;   {* lp->lpt_orig *}
   LPX_K_PRESOL    =  327;   {* lp->presol *}

   MaxIntArray     = 1000000;
   MaxFloatArray   = 1000000;

Type
 PLPX     = Pointer;
 PLPXKKT  = Pointer;

 PIntArray      = ^IntArray;
 IntArray       = Array[0..MaxIntArray] Of Integer;

 PFloatArray    = ^FloatArray;
 FloatArray     = Array[0..MaxFloatArray] Of Double;


{* void lpx_add_cols(LPX *lp, int ncs); *}
{* add new columns to LP problem object *}
Procedure lpx_add_cols(lp : PLPX; ncs : Integer);
cdecl; external 'glpk-4.dll' index 177;

{* void lpx_add_rows(LPX *lp, int nrs); *}
{* add new rows to LP problem object *}
Procedure lpx_add_rows(lp : PLPX; nrs : Integer);
cdecl; external 'glpk-4.dll' index 178;

{* void lpx_adv_basis(LPX *lp); *}
{* build advanced initial basis *}
Procedure lpx_adv_basis(lp : PLPX);
cdecl; external 'glpk-4.dll' index 179;

{* void lpx_check_kkt(LPX *lp, int scaled, LPXKKT *kkt); *}
{* check Karush-Kuhn-Tucker conditions *}
Procedure lpx_check_kkt(lp : PLPX; scaled : Integer; kkt : PLPXKKT);
cdecl; external 'glpk-4.dll' index 180;

{* int lpx_check_name(char *name); *}
{* check correctness of symbolic name *}
Function lpx_check_name(name : String):Integer;
cdecl; external 'glpk-4.dll' index 181;

{* void lpx_clear_mat(LPX *lp); *}
{* clear rows and columns of the constraint matrix *}
Procedure lpx_clear_mat(lp : PLPX);
cdecl; external 'glpk-4.dll' index 182;

{* LPX *lpx_create_prob(void); *}
{* create LP problem object *}
Function lpx_create_prob:PLPX;
cdecl; external 'glpk-4.dll' index 183;

{* void lpx_del_items(LPX *lp); *}
{* remove rows and columns from LP problem object *}
Procedure lpx_del_items(lp : PLPX);
cdecl; external 'glpk-4.dll' index 184;

{* void lpx_delete_prob(LPX *lp); *}
{* delete LP problem object *}
Procedure lpx_delete_prob(lp : PLPX);
cdecl; external 'glpk-4.dll' index 185;

{* int lpx_dual_opt(LPX *lp); *}
{* find optimal solution (dual simplex) *}
Function lpx_dual_opt(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 186;

{* int lpx_dual_ratio_test(LPX *lp, int len, int ndx[], double val[], int how, double tol); *}
{* perform dual ratio test *}
Function lpx_dual_ratio_test(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray; how : Integer; tol : Double):Integer;
cdecl; external 'glpk-4.dll' index 187;

{* void lpx_estim_obj_chg(LPX *lp, int k, double dn_val, double up_val, double *dn_chg, double *up_chg, int ndx[], double val[]); *}
{* estimate obj. changes for down- and up-branch *}
Procedure lpx_estim_obj_chg(lp : PLPX; k : Integer; dn_val : Double; up_val : Double; Var dn_chg : Double; Var up_chg : Double; ndx : PIntArray; val : PFloatArray);
cdecl; external 'glpk-4.dll' index 188;

{* double lpx_eval_activity(LPX *lp, int len, int ndx[], double val[]); *}
{* compute activity of explicitly specified row *}
Function lpx_eval_activity(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray):Double;
cdecl; external 'glpk-4.dll' index 189;

{* double lpx_eval_red_cost(LPX *lp, int len, int ndx[], double val[]); *}
{* compute red. cost of explicitly specified column *}
Function lpx_eval_red_cost(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray):Double;
cdecl; external 'glpk-4.dll' index 190;

{* int lpx_eval_tab_col(LPX *lp, int k, int ndx[], double val[]); *}
{* compute column of the simplex table *}
Function lpx_eval_tab_col(lp : PLPX; k : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 191;

{* int lpx_eval_tab_row(LPX *lp, int k, int ndx[], double val[]); *}
{* compute row of the simplex table *}
Function lpx_eval_tab_row(lp : PLPX; k : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 192;

{* int lpx_get_class(LPX *lp); *}
{* query problem class *}
Function lpx_get_class(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 193;

{* void lpx_get_col_bnds(LPX *lp, int j, int *typx, double *lb, double *ub); *}
{* obtain column bounds *}
Procedure lpx_get_col_bnds(lp : PLPX; j : Integer; Var typx : Integer; Var lb : Double; Var ub : Double);
cdecl; external 'glpk-4.dll' index 194;

{* double lpx_get_col_coef(LPX *lp, int j); *}
{* obtain column objective coefficient *}
Function lpx_get_col_coef(lp : PLPX; j : Integer):Double;
cdecl; external 'glpk-4.dll' index 195;

{* void lpx_get_col_info(LPX *lp, int j, int *tagx, double *vx, double *dx); *}
{* obtain column solution information *}
Procedure lpx_get_col_info(lp : PLPX; j : Integer; Var tagx : Integer; Var vx : Double; Var dx : Double);
cdecl; external 'glpk-4.dll' index 196;

{* int lpx_get_col_kind(LPX *lp, int j); *}
{* query column kind *}
Function lpx_get_col_kind(lp : PLPX; j : Integer):Integer;
cdecl; external 'glpk-4.dll' index 197;

{* int lpx_get_col_mark(LPX *lp, int j); *}
{* determine column mark *}
Function lpx_get_col_mark(lp : PLPX; j : Integer):Integer;
cdecl; external 'glpk-4.dll' index 198;

{* char *lpx_get_col_name(LPX *lp, int j); *}
{* obtain column name *}
Function lpx_get_col_name(lp : PLPX; j : Integer):PChar;
cdecl; external 'glpk-4.dll' index 199;

{* int lpx_get_dual_stat(LPX *lp); *}
{* query dual status of basic solution *}
Function lpx_get_dual_stat(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 200;

{* int lpx_get_int_parm(LPX *lp, int parm); *}
{* query integer control parameter *}
Function lpx_get_int_parm(lp : PLPX; parm : Integer):Integer;
cdecl; external 'glpk-4.dll' index 201;

{* void lpx_get_ips_col(LPX *lp, int j, double *vx, double *dx); *}
{* obtain column interior point solution *}
Procedure lpx_get_ips_col(lp : PLPX; j : Integer; Var vx : Double; Var dx : Double);
cdecl; external 'glpk-4.dll' index 202;

{* double lpx_get_ips_obj(LPX *lp); *}
{* obtain interior point value of the obj. function *}
Function lpx_get_ips_obj(lp : PLPX):Double;
cdecl; external 'glpk-4.dll' index 203;

{* void lpx_get_ips_row(LPX *lp, int i, double *vx, double *dx); *}
{* obtain row interior point solution *}
Procedure lpx_get_ips_row(lp : PLPX; i : Integer; Var vx : Double; Var dx : Double);
cdecl; external 'glpk-4.dll' index 204;

{* int lpx_get_ips_stat(LPX *lp); *}
{* query status of interior point solution *}
Function lpx_get_ips_stat(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 205;

{* int lpx_get_mat_col(LPX *lp, int j, int ndx[], double val[]); *}
{* get column of the constraint matrix *}
Function lpx_get_mat_col(lp : PLPX; j : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 206;

{* int lpx_get_mat_row(LPX *lp, int i, int ndx[], double val[]); *}
{* get row of the constraint matrix *}
Function lpx_get_mat_row(lp : PLPX; i : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 207;

{* double lpx_get_mip_col(LPX *lp, int j); *}
{* obtain column activity for MIP solution *}
Function lpx_get_mip_col(lp : PLPX; j : Integer):Double;
cdecl; external 'glpk-4.dll' index 208;

{* double lpx_get_mip_obj(LPX *lp); *}
{* obtain value of the obj. func. for MIP solution *}
Function lpx_get_mip_obj(lp : PLPX):Double;
cdecl; external 'glpk-4.dll' index 209;

{* double lpx_get_mip_row(LPX *lp, int i); *}
{* obtain row activity for MIP solution *}
Function lpx_get_mip_row(lp : PLPX; i : Integer):Double;
cdecl; external 'glpk-4.dll' index 210;

{* int lpx_get_mip_stat(LPX *lp); *}
{* query status of MIP solution *}
Function lpx_get_mip_stat(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 211;

{* int lpx_get_num_bin(LPX *lp); *}
{* determine number of binary columns *}
Function lpx_get_num_bin(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 212;

{* int lpx_get_num_cols(LPX *lp); *}
{* determine number of columns *}
Function lpx_get_num_cols(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 213;

{* int lpx_get_num_int(LPX *lp); *}
{* determine number of integer columns *}
Function lpx_get_num_int(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 214;

{* int lpx_get_num_nz(LPX *lp); *}
{* determine number of constraint coefficients *}
Function lpx_get_num_nz(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 215;

{* int lpx_get_num_rows(LPX *lp); *}
{* determine number of rows *}
Function lpx_get_num_rows(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 216;

{* double lpx_get_obj_c0(LPX *lp); *}
{* obtain constant term of the obj. function *}
Function lpx_get_obj_c0(lp : PLPX):Double;
cdecl; external 'glpk-4.dll' index 217;

{* int lpx_get_obj_dir(LPX *lp); *}
{* determine optimization direction *}
Function lpx_get_obj_dir(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 218;

{* char *lpx_get_obj_name(LPX *lp); *}
{* obtain objective function name *}
Function lpx_get_obj_name(lp : PLPX):PChar;
cdecl; external 'glpk-4.dll' index 219;

{* double lpx_get_obj_val(LPX *lp); *}
{* obtain value of the objective function *}
Function lpx_get_obj_val(lp : PLPX):Double;
cdecl; external 'glpk-4.dll' index 220;

{* int lpx_get_prim_stat(LPX *lp); *}
{* query primal status of basic solution *}
Function lpx_get_prim_stat(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 221;

{* char *lpx_get_prob_name(LPX *lp); *}
{* obtain problem name *}
Function lpx_get_prob_name(lp : PLPX):PChar;
cdecl; external 'glpk-4.dll' index 222;

{* double lpx_get_real_parm(LPX *lp, int parm); *}
{* query real control parameter *}
Function lpx_get_real_parm(lp : PLPX; parm : Integer):Double;
cdecl; external 'glpk-4.dll' index 223;

{* void lpx_get_row_bnds(LPX *lp, int i, int *typx, double *lb, double *ub); *}
{* obtain row bounds *}
Procedure lpx_get_row_bnds(lp : PLPX; i : Integer; Var typx : Integer; Var lb : Double; Var ub : Double);
cdecl; external 'glpk-4.dll' index 224;

{* double lpx_get_row_coef(LPX *lp, int i); *}
{* obtain row objective coefficient *}
Function lpx_get_row_coef(lp : PLPX; i : Integer):Double;
cdecl; external 'glpk-4.dll' index 225;

{* void lpx_get_row_info(LPX *lp, int i, int *tagx, double *vx, double *dx); *}
{* obtain row solution information *}
Procedure lpx_get_row_info(lp : PLPX; i : Integer; Var tagx : Integer; Var vx : Double; Var dx : Double);
cdecl; external 'glpk-4.dll' index 226;

{* int lpx_get_row_mark(LPX *lp, int i); *}
{* determine row mark *}
Function lpx_get_row_mark(lp : PLPX; i : Integer):Integer;
cdecl; external 'glpk-4.dll' index 227;

{* char *lpx_get_row_name(LPX *lp, int i); *}
{* obtain row name *}
Function lpx_get_row_name(lp : PLPX; i : Integer):PChar;
cdecl; external 'glpk-4.dll' index 228;

{* int lpx_get_status(LPX *lp); *}
{* query basic solution status *}
Function lpx_get_status(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 229;

{* int lpx_integer(LPX *lp); *}
{* easy-to-use driver to the branch-and-bound method *}
Function lpx_integer(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 230;

{* int lpx_interior(LPX *lp); *}
{* easy-to-use driver to the interior point method *}
Function lpx_interior(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 231;

{* void lpx_load_mat3(LPX *lp, int nz, int rn[], int cn[], double a[]); *}
{* load constraint matrix *}
Procedure lpx_load_mat3(lp : PLPX; nz : Integer; rn : PIntArray; cn : PIntArray; a : PFloatArray);
cdecl; external 'glpk-4.dll' index 233;

{* void lpx_mark_col(LPX *lp, int j, int mark); *}
{* assign mark to column *}
Procedure lpx_mark_col(lp : PLPX; j : Integer; mark : Integer);
cdecl; external 'glpk-4.dll' index 234;

{* void lpx_mark_row(LPX *lp, int i, int mark); *}
{* assign mark to row *}
Procedure lpx_mark_row(lp : PLPX; i : Integer; mark : Integer);
cdecl; external 'glpk-4.dll' index 235;

{* int lpx_mixed_gomory(LPX *lp, int kind[], int len, int ndx[], double val[], double work[]); *}
{* generate Gomory's mixed integer cut *}
Function lpx_mixed_gomory(lp : PLPX; kind : PIntArray; len : Integer; ndx : PIntArray; val : PFloatArray; work : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 236;

{* int lpx_prim_art(LPX *lp); *}
{* find primal feasible solution (primal simplex) *}
Function lpx_prim_art(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 237;

{* int lpx_prim_opt(LPX *lp); *}
{* find optimal solution (primal simplex) *}
Function lpx_prim_opt(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 238;

{* int lpx_prim_ratio_test(LPX *lp, int len, int ndx[], double val[], int how, double tol); *}
{* perform primal ratio test *}
Function lpx_prim_ratio_test(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray; how : Integer; tol : Double):Integer;
cdecl; external 'glpk-4.dll' index 239;

{* int lpx_print_ips(LPX *lp, char *fname); *}
{* write interior point solution in printable format *}
Function lpx_print_ips(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 240;

{* int lpx_print_mip(LPX *lp, char *fname); *}
{* write MIP problem solution in printable format *}
Function lpx_print_mip(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 241;

{* int lpx_print_prob(LPX *lp, char *fname); *}
{* write problem data in plain text format *}
Function lpx_print_prob(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 242;

{* int lpx_print_sol(LPX *lp, char *fname); *}
{* write LP problem solution in printable format *}
Function lpx_print_sol(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 243;

{* int lpx_read_bas(LPX *lp, char *fname); *}
{* read predefined basis using MPS format *}
Function lpx_read_bas(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 244;

{* LPX *lpx_read_lpt(char *fname); *}
{* read LP/MIP problem using CPLEX LP format *}
Function lpx_read_lpt(fname : String):PLPX;
cdecl; external 'glpk-4.dll' index 245;

{* LPX *lpx_read_model(char *model, char *data, char *output); *}
{* read LP/MIP model written in GNU MathProg language *}
Function lpx_read_model(model : String; data : String; output : String):PLPX;
cdecl; external 'glpk-4.dll' index 246;

{* LPX *lpx_read_mps(char *fname); *}
{* read in problem data using MPS format *}
Function lpx_read_mps(fname : String):PLPX;
cdecl; external 'glpk-4.dll' index 247;

{* void lpx_realloc_prob(LPX *lp, int m_max, int n_max); *}
{* reallocate LP problem object *}
Procedure lpx_realloc_prob(lp : PLPX; m_max : Integer; n_max : Integer);
cdecl; external 'glpk-4.dll' index 248;

{* int lpx_reduce_form(LPX *lp, int len, int ndx[], double val[], double work[]); *}
{* reduce linear form *}
Function lpx_reduce_form(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray; work : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 249;

{* void lpx_reset_parms(LPX *lp); *}
{* reset control parameters to default values *}
Procedure lpx_reset_parms(lp : PLPX);
cdecl; external 'glpk-4.dll' index 250;

{* void lpx_scale_prob(LPX *lp); *}
{* scale LP problem data *}
Procedure lpx_scale_prob(lp : PLPX);
cdecl; external 'glpk-4.dll' index 251;

{* void lpx_set_class(LPX *lp, int clss); *}
{* set (change) problem class *}
Procedure lpx_set_class(lp : PLPX; clss : Integer);
cdecl; external 'glpk-4.dll' index 252;

{* void lpx_set_col_bnds(LPX *lp, int j, int typx, double lb, double ub); *}
{* set (change) column bounds *}
Procedure lpx_set_col_bnds(lp : PLPX; j : Integer; typx : Integer; lb : Double; ub : Double);
cdecl; external 'glpk-4.dll' index 253;

{* void lpx_set_col_coef(LPX *lp, int j, double coef); *}
{* set (change) column objective coefficient *}
Procedure lpx_set_col_coef(lp : PLPX; j : Integer; coef : Double);
cdecl; external 'glpk-4.dll' index 254;

{* void lpx_set_col_kind(LPX *lp, int j, int kind); *}
{* set (change) column kind *}
Procedure lpx_set_col_kind(lp : PLPX; j : Integer; kind : Integer);
cdecl; external 'glpk-4.dll' index 255;

{* void lpx_set_col_name(LPX *lp, int j, char *name); *}
{* assign (change) column name *}
Procedure lpx_set_col_name(lp : PLPX; j : Integer; name : String);
cdecl; external 'glpk-4.dll' index 256;

{* void lpx_set_col_stat(LPX *lp, int j, int stat); *}
{* set (change) column status *}
Procedure lpx_set_col_stat(lp : PLPX; j : Integer; stat : Integer);
cdecl; external 'glpk-4.dll' index 257;

{* void lpx_set_int_parm(LPX *lp, int parm, int val); *}
{* set (change) integer control parameter *}
Procedure lpx_set_int_parm(lp : PLPX; parm : Integer; val : Integer);
cdecl; external 'glpk-4.dll' index 258;

{* void lpx_set_mat_col(LPX *lp, int j, int len, int ndx[], double val[]); *}
{* set (replace) column of the constraint matrix *}
Procedure lpx_set_mat_col(lp : PLPX; j : Integer; len : Integer; ndx : PIntArray; val : PFloatArray);
cdecl; external 'glpk-4.dll' index 259;

{* void lpx_set_mat_row(LPX *lp, int i, int len, int ndx[], double val[]); *}
{* set (replace) row of the constraint matrix *}
Procedure lpx_set_mat_row(lp : PLPX; i : Integer; len : Integer; ndx : PIntArray; val : PFloatArray);
cdecl; external 'glpk-4.dll' index 260;

{* void lpx_set_obj_c0(LPX *lp, double c0); *}
{* set (change) constant term of the obj. function *}
Procedure lpx_set_obj_c0(lp : PLPX; c0 : Double);
cdecl; external 'glpk-4.dll' index 261;

{* void lpx_set_obj_dir(LPX *lp, int dir); *}
{* set (change) optimization direction *}
Procedure lpx_set_obj_dir(lp : PLPX; dir : Integer);
cdecl; external 'glpk-4.dll' index 262;

{* void lpx_set_obj_name(LPX *lp, char *name); *}
{* assign (change) objective function name *}
Procedure lpx_set_obj_name(lp : PLPX; name : String);
cdecl; external 'glpk-4.dll' index 263;

{* void lpx_set_prob_name(LPX *lp, char *name); *}
{* assign (change) problem name *}
Procedure lpx_set_prob_name(lp : PLPX; name : String);
cdecl; external 'glpk-4.dll' index 264;

{* void lpx_set_real_parm(LPX *lp, int parm, double val); *}
{* set (change) real control parameter *}
Procedure lpx_set_real_parm(lp : PLPX; parm : Integer; val : Double);
cdecl; external 'glpk-4.dll' index 265;

{* void lpx_set_row_bnds(LPX *lp, int i, int typx, double lb, double ub); *}
{* set (change) row bounds *}
Procedure lpx_set_row_bnds(lp : PLPX; i : Integer; typx : Integer; lb : Double; ub : Double);
cdecl; external 'glpk-4.dll' index 266;

{* void lpx_set_row_coef(LPX *lp, int i, double coef); *}
{* set (change) row objective coefficient *}
Procedure lpx_set_row_coef(lp : PLPX; i : Integer; coef : Double);
cdecl; external 'glpk-4.dll' index 267;

{* void lpx_set_row_name(LPX *lp, int i, char *name); *}
{* assign (change) row name *}
Procedure lpx_set_row_name(lp : PLPX; i : Integer; name : String);
cdecl; external 'glpk-4.dll' index 268;

{* void lpx_set_row_stat(LPX *lp, int i, int stat); *}
{* set (change) row status *}
Procedure lpx_set_row_stat(lp : PLPX; i : Integer; stat : Integer);
cdecl; external 'glpk-4.dll' index 269;

{* int lpx_simplex(LPX *lp); *}
{* easy-to-use driver to the simplex method *}
Function lpx_simplex(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 270;

{* void lpx_std_basis(LPX *lp); *}
{* build standard initial basis *}
Procedure lpx_std_basis(lp : PLPX);
cdecl; external 'glpk-4.dll' index 271;

{* int lpx_transform_col(LPX *lp, int len, int ndx[], double val[]); *}
{* transform explicitly specified column *}
Function lpx_transform_col(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 272;

{* int lpx_transform_row(LPX *lp, int len, int ndx[], double val[]); *}
{* transform explicitly specified row *}
Function lpx_transform_row(lp : PLPX; len : Integer; ndx : PIntArray; val : PFloatArray):Integer;
cdecl; external 'glpk-4.dll' index 273;

{* void lpx_unmark_all(LPX *lp); *}
{* unmark all rows and columns *}
Procedure lpx_unmark_all(lp : PLPX);
cdecl; external 'glpk-4.dll' index 274;

{* void lpx_unscale_prob(LPX *lp); *}
{* unscale LP problem data *}
Procedure lpx_unscale_prob(lp : PLPX);
cdecl; external 'glpk-4.dll' index 275;

{* int lpx_warm_up(LPX *lp); *}
{* "warm up" the initial basis *}
Function lpx_warm_up(lp : PLPX):Integer;
cdecl; external 'glpk-4.dll' index 276;

{* int lpx_write_bas(LPX *lp, char *fname); *}
{* write current basis using MPS format *}
Function lpx_write_bas(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 277;

{* int lpx_write_lpt(LPX *lp, char *fname); *}
{* write problem data using CPLEX LP format *}
Function lpx_write_lpt(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 278;

{* int lpx_write_mps(LPX *lp, char *fname); *}
{* write problem data using MPS format *}
Function lpx_write_mps(lp : PLPX; fname : String):Integer;
cdecl; external 'glpk-4.dll' index 279;

Function PutZTString( S : String): String;


Implementation


{* Put Zero Terminated String  to DLL *}
Function PutZTString( S : String): String;

 Begin
   PutZTString := S + Chr(0);
 End;



End.

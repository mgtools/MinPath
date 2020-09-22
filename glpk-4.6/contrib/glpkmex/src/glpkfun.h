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

int glpkmex_fault_hook(void *info,  char *msg)
{
	 char errmsg[1024];
	 
	 sprintf(errmsg,"*** SEVERE CRITICAL ERROR *** from GLPK !\n%s\n",msg);
	 mexErrMsgTxt(errmsg);
	 longjmp( mark, -1 );
}

int glpkmex_print_hook(void *info,  char *msg)
{
	 mexPrintf("%s\n",msg);
	 return 1;
}



int glpk(int sense,int n, int m, double *c,int nz,int *rn,int *cn, double *a,double *b, char *ctype,
         int *freeLB, double *lb, int *freeUB, double *ub, int *vartype,
         int isMIP, int lpsolver,int save_pb, double *xmin, double *fmin, double *status,
         double *lambda, double *redcosts, double *time, double *mem)
{
   LPX *lp;
   int i,j;
   int error;
   double t_start;
   int len;
   int typx;
   int method;

   t_start = utime();
   
   lib_set_fault_hook(NULL,glpkmex_fault_hook);
   
   if (lpxIntParam[0] > 1){
	lib_set_print_hook(NULL,glpkmex_print_hook);
   }
   
   lp=lpx_create_prob();
     
   if(isMIP) lpx_set_class(lp,LPX_MIP);
  
   lpx_add_cols(lp,n);   
   for(i=0;i<n;i++){
      lpx_set_col_coef(lp,i+1,c[i]);
      if (!freeLB[i] && !freeUB[i]){
        lpx_set_col_bnds(lp,i+1,LPX_DB,lb[i],ub[i]);
      }else{
         if (!freeLB[i] && freeUB[i]){
            lpx_set_col_bnds(lp,i+1,LPX_LO,lb[i],ub[i]);
         }else{
            if (freeLB[i] && !freeUB[i]){
               lpx_set_col_bnds(lp,i+1,LPX_UP,lb[i],ub[i]);
            }else{
               lpx_set_col_bnds(lp,i+1,LPX_FR,lb[i],ub[i]);
            }
         }
      }
      if(isMIP){
        lpx_set_col_kind(lp,i+1,vartype[i]);
      }
   }
   
   lpx_add_rows(lp,m);   
   for(i=0;i<m;i++){
   
      /* If the i-th row has no lower bound (types F,U), the
	      corrispondent parameter will be ignored.
         If the i-th row has no upper bound (types F,L), the corrispondent
         parameter will be ignored.
         If the i-th row is of S type, the i-th LB is used, but
         the i-th UB is ignored.
      */
      switch(ctype[i]){
      	case 'F': typx=LPX_FR; break;
	      case 'U': typx=LPX_UP; break;
	      case 'L': typx=LPX_LO; break;
	      case 'S': typx=LPX_FX; break;
	      case 'D': typx=LPX_DB; 
      }
      lpx_set_row_bnds(lp,i+1,typx,b[i],b[i]);

   }
   lpx_load_mat3(lp,nz,rn,cn,a);

   if (sense==1) lpx_set_obj_dir(lp,LPX_MIN);
   else lpx_set_obj_dir(lp,LPX_MAX);

   for(i=0;i<NIntP;i++){
     lpx_set_int_parm(lp,IParam[i],lpxIntParam[i]);
   }
   for(i=0;i<NRealP;i++){
     lpx_set_real_parm(lp,RParam[i],lpxRealParam[i]);
   }

   /* scale the problem data (if required) */
   /* if (scale && (!presol || method == 1)) lpx_scale_prob(lp); */
   /*  LPX_K_SCALE=IParam[1]  LPX_K_PRESOL=IParam[16]  */
   if (lpxIntParam[1] && (!lpxIntParam[16] || lpsolver!=1)){
	   lpx_scale_prob(lp);
   }
   /* build advanced initial basis (if required) */
   if (lpsolver == 1 && !lpxIntParam[16]){
	   lpx_adv_basis(lp);
   }
   
   if (save_pb){
   	if(lpx_write_lpt(lp, "outpb.lp") != 0)
      	mexErrMsgTxt("Unable to write problem");
   }

	if(lpsolver==1) method='S';
   else method='T';


   switch(method){
   case 'S':
     if(isMIP){
       method='I';
       error=lpx_simplex(lp);
       error=lpx_integer(lp);
     }else{
       error=lpx_simplex(lp);
     }
     break;
   case 'T':
     error=lpx_interior(lp);
     break;
   default:
     insist(method != method);
   }

  	/*
       error assumes the following results:
       error=0 <=> No errors
       error=1 <=> Iteration limit exceeded.
       error=2 <=> Numerical problems with basis matrix.
   */
   if(error==LPX_E_OK){
     if(isMIP){
       *status=(double)lpx_get_mip_stat(lp);
       *fmin=lpx_get_mip_obj(lp);
     }else{
       if(lpsolver==1){
	      *status=(double)lpx_get_status(lp);
	      *fmin=lpx_get_obj_val(lp);
       }else{
	      *status=(double)lpx_get_ips_stat(lp);
	      *fmin=lpx_get_ips_obj(lp);
       }
     }
     if(isMIP){
         for(i=0;i<n;i++)  xmin[i]=lpx_get_mip_col(lp,i+1);
     }else{
      /* Primal values */
      for(i=0;i<n;i++){
         if(lpsolver==1) xmin[i]=lpx_get_col_prim(lp,i+1);
	      else xmin[i]=lpx_ipt_col_prim(lp,i+1);
	   }
	   /* Dual values */
	   for(i=0; i<m; i++){
	      if(lpsolver==1) lambda[i]=lpx_get_row_dual(lp,i+1);
	      else lambda[i]=lpx_ipt_row_dual(lp,i+1);
      }
      /* Reduced costs */
      for(i=0; i<lpx_get_num_cols(lp); i++){
	      if(lpsolver==1) redcosts[i]=lpx_get_col_dual(lp,i+1);
	      else redcosts[i]=lpx_ipt_col_dual(lp,i+1);
      }
     }
     *time=(double)(utime() - t_start);
     *mem=(double)lib_env_ptr()->mem_tpeak;      

	  lpx_delete_prob(lp);
     return(0);
   }
	*status=(double)error;
   return(error);
}


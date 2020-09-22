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
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <setjmp.h>
#include "mex.h"
#include "glpk.h"

#include "glpksets.h"
#include "glpkfun.h"

#define	SENSE_IN   prhs[0]
#define	C_IN	    prhs[1]
#define	A_IN	    prhs[2]
#define	B_IN	    prhs[3]
#define  CTYPE_IN   prhs[4]
#define  LB_IN	    prhs[5]
#define  UB_IN      prhs[6]
#define  VARTYPE_IN prhs[7]
#define  PARAM      prhs[8]
#define  SOLVER_IN  prhs[9]
#define  SAVE_IN    prhs[10]


/* Output Arguments */
#define	 XMIN_OUT     plhs[0]
#define	 FMIN_OUT     plhs[1]
#define	 STATUS_OUT   plhs[2]
#define   EXTRA_OUT    plhs[3]


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   int    sense;
   double *tmp=NULL;
   double *c=NULL;
   double *A=NULL;
   double *b=NULL;
   double *lb=NULL;
   double *ub=NULL;
   char *vartype=NULL;
   int mrowsc=0;
   int mrowsA=0;
   char *ctype=NULL;
   int save_pb=0;
   int lpsolver=1;
   int error;
   double *xmin;
   double *fmin;
   double *status;
   double *lambda;
   double *redcosts;
   double *time;
   double *mem;
   double *extra;
   char errmsg[1024];
   int *freeLB=0;
   int *freeUB=0;
   int nfields=0;
   const char **extranames;
   mxArray *mxlambda, *mxredcosts, *mxtime,*mxmem;
   double *rdtmp=NULL;
   mxArray *mxtmp;
   int *vartype2;
   int jmpret;

   /* row and column sets of non-zero constraint coefficients */
   int nz;  /* number of non-zero coefficients */
   int *rn;
   int *cn;
   double *a;  /* non-zero coefficients */

   /*
      flag to identify the type of problem:
             isMIP=0 <=> LP problem
             isMIP=1 <=> MIP problem
   */
   int isMIP = 0;

   if(nrhs < 1){
      mexPrintf(header);
      mexPrintf(version);
      mexPrintf(copyright);
      mexPrintf(syntax);
      return;
   }
   if(nrhs < 4) mexErrMsgTxt("At least 4 inputs required (SENSE,C,A,b)");
   if(nlhs < 2) mexErrMsgTxt("2 outputs required");


   /* 1st Input. Sense of optimization. */
   if (!mxIsNumeric(SENSE_IN)
       || (mxGetNumberOfDimensions(SENSE_IN) > 2)
       || (mxGetM(SENSE_IN) != 1)
       || (mxGetN(SENSE_IN) != 1)
       || mxIsComplex(SENSE_IN)
       || ((tmp = mxGetPr(SENSE_IN)) == NULL)
       || ((*tmp != 1) && (*tmp != -1))
      )
       mexErrMsgTxt("SENSE must be either 1 or -1.");
   else
       sense = (int) (*tmp);

   /* 2nd Input. A column array containing the objective function
                 coefficients.
   */
   if (!mxIsNumeric(C_IN)
	    || (mxGetNumberOfDimensions(C_IN) > 2)
	    || ((mrowsc = mxGetM(C_IN)) < 1)
	    || (mxGetN(C_IN) != 1)
	    || mxIsComplex(C_IN)
	    || ((c= mxGetPr(C_IN)) == NULL)
		)
		mexErrMsgTxt("C must be a real valued column vector.");

   /* 3rd Input. A matrix containing the constraints coefficients. */
	if (!mxIsNumeric(A_IN)
	    || (mxGetNumberOfDimensions(A_IN) > 2)
	    || ((mrowsA = mxGetM(A_IN)) < 1)
	    || (mxGetN(A_IN) != mrowsc)
	    || mxIsComplex(A_IN)
	    || ((A = mxGetPr(A_IN)) == NULL)
		) {
		sprintf(errmsg,"A must be a real valued %d by %d matrix.",mrowsA, mrowsc);
		mexErrMsgTxt(errmsg);
	}else{

	  if(!mxIsSparse(A_IN)){
	    int i,j;
	    int nrcount;

	    nz=0;
	    for(i=0;i<mrowsA;i++){
	      for(j=0;j<mrowsc;j++){
		      if(A[i+j*mrowsA]!=0) nz++;
	      }
	    }
	    rn=(int *)mxCalloc(nz+1,sizeof(int));
	    cn=(int *)mxCalloc(nz+1,sizeof(int));
	    a=(double *)mxCalloc(nz+1,sizeof(double));

	    nrcount=0;
	    for(i=0;i<mrowsA;i++){
	      for(j=0;j<mrowsc;j++){
		      if(A[i+j*mrowsA]!=0){
                  nrcount++;
		            rn[nrcount]=i+1;
		            cn[nrcount]=j+1;
		            a[nrcount]=A[i+j*mrowsA];
		      }
	      }
	    }
	  }else{
	    int i,j;
	    int *jc,*ir;
	    double *pr;
	    int nelc,count,row;

	    /* NOTE: nnz is the actual number of nonzeros and is stored as the
       last element of the jc array where the size of the jc array is the
       number of columns + 1 */
	    nz = *(mxGetJc(A_IN) + mrowsc);
	    jc = mxGetJc(A_IN);
	    ir = mxGetIr(A_IN);
	    pr = mxGetPr(A_IN);

       rn=(int *)mxCalloc(nz+1,sizeof(int));
	    cn=(int *)mxCalloc(nz+1,sizeof(int));
	    a=(double *)mxCalloc(nz+1,sizeof(double));

       count=0; row=0;
	    for(i=1;i<=mrowsc;i++){
	      nelc=jc[i]-jc[i-1];
	      for(j=0;j<nelc;j++){
		      count++;
		      rn[count]=ir[row]+1;
		      cn[count]=i;
		      a[count]=pr[row];
		      row++;
	      }
	    }
	  }

	}

   /* 4th Input. A column array containing the right-hand side value
	         for each constraint in the constraint matrix.
   */
	if (!mxIsNumeric(B_IN)
	    || (mxGetNumberOfDimensions(B_IN) > 2)
	    || (mxGetM(B_IN) != mrowsA)
	    || (mxGetN(B_IN) != 1)
	    || mxIsComplex(B_IN)
	    || ((b = mxGetPr(B_IN)) == NULL)
		) {
		sprintf(errmsg,"B must be a real valued %d by 1 column vector.",mrowsA);
		mexErrMsgTxt(errmsg);
	}

   /* 5th Input. A column array containing the sense of each constraint
                 in the constraint matrix.
   */
  if ((nrhs > 4) && (mxGetM(CTYPE_IN) != 0) && (mxGetN(CTYPE_IN) != 0)) {
  if (!mxIsChar(CTYPE_IN)
       || (mxGetNumberOfDimensions(CTYPE_IN) > 2)
       || (mxGetM(CTYPE_IN) != mrowsA)
       || (mxGetN(CTYPE_IN) != 1)
       || mxIsComplex(CTYPE_IN)
      ){
	sprintf(errmsg,"CTYPE must be a char valued %d by 1 column vector.",mrowsA);
	mexErrMsgTxt(errmsg);
       } else {
	int i,size;

	size = mxGetNumberOfElements(CTYPE_IN) + 1;

	/* Allocate enough memory to hold the converted string. */
	ctype = mxCalloc(size, sizeof (char));

	/* Copy the string data from string_array_ptr and place it into buf. */
	if (mxGetString(CTYPE_IN, ctype, size) != 0)
	  mexErrMsgTxt("Could not convert string data.");

	/* checking if the input is made only of  F, U, S, L and D */
	for (i = 0; i < size - 1; i++) {
	  if (     (ctype[i] != 'F') && (ctype[i] != 'U')
		&& (ctype[i] != 'S') && (ctype[i] != 'L')
		&& (ctype[i] != 'D'))
	          mexErrMsgTxt("CTYPE must contain only F,U,S,L and D");
	}
       }
  }

   /* 6th Input. An array of at least length numcols containing the lower
            	  bound on each of the variables.
   */
    if ((nrhs > 5) && (mxGetM(LB_IN) != 0) && (mxGetN(LB_IN) != 0)) {
    if (!mxIsNumeric(LB_IN)
	|| (mxGetNumberOfDimensions(LB_IN) > 2)
	|| (mrowsc != mxGetM(LB_IN))
	|| (mxGetN(LB_IN) != 1)
	|| mxIsComplex(LB_IN)
	|| ((lb = mxGetPr(LB_IN)) == NULL)
       ) {
      sprintf(errmsg,"LB must be a real valued %d by 1 column vector.",mrowsc);
      mexErrMsgTxt(errmsg);
    }
    }

    /* 7th Input. An array of at least length numcols containing the upper
       bound on each of the variables.
    */
    if ((nrhs > 6) && (mxGetM(UB_IN) != 0) && (mxGetN(UB_IN) != 0)) {
      if (!mxIsNumeric(UB_IN)
	  || (mxGetNumberOfDimensions(UB_IN) > 2)
	  || (mrowsc != mxGetM(UB_IN))
	  || (mxGetN(UB_IN) != 1)
	  || mxIsComplex(UB_IN)
	  || ((ub = mxGetPr(UB_IN)) == NULL)
	  ) {
	sprintf(errmsg,"UB must be a real valued %d by 1 column vector.",
		mrowsc);
	mexErrMsgTxt(errmsg);
      }
    }

    /* 8th Input. A column array containing the types of the variables.
     */
    if ((nrhs > 7) && (mxGetM(VARTYPE_IN) != 0)
	&& (mxGetN(VARTYPE_IN) != 0)) {
      if (!mxIsChar(VARTYPE_IN)
	  || (mxGetNumberOfDimensions(VARTYPE_IN) > 2)
	  || (mxGetM(VARTYPE_IN) != mrowsc)
	  || (mxGetN(VARTYPE_IN) != 1)
	  || mxIsComplex(VARTYPE_IN)
	  ) {
	sprintf(errmsg, "VARTYPE must be a char valued %d by 1 column vector.", mrowsc);
	mexErrMsgTxt(errmsg);
      } else {
	int i,size;

	size = mxGetNumberOfElements(VARTYPE_IN)+1;
	/*   Allocate enough memory to hold the converted string.
	 */
	vartype   = mxCalloc(size, sizeof (char));
	vartype2 = mxCalloc(size, sizeof (int));

	/* Copy the string data from string_array_ptr and place it into buf.
	 */
	if (mxGetString(VARTYPE_IN, vartype, size) != 0)
	  mexErrMsgTxt("Could not convert string data.");

	/* checking if the input is made only of C I */
	for (i = 0; i < size-1 ; i++){
	  if ((vartype[i] != 'C') && (vartype[i] != 'I'))
	    mexErrMsgTxt("VARTYPE must contain only LPX_CV or LPX_IV");
	  if(vartype[i]=='I'){
	    isMIP=1;
	    vartype2[i]=LPX_IV;
	  }else{
	    vartype2[i]=LPX_CV;
	  }
	}
      }
    } else mexWarnMsgTxt ("Omitting VARTYPE you are assuming all variables are continuous");

    /*
           9th Input. A structure containing the control parameters.
	*/
	if ((nrhs > 8) && ((nfields=mxGetNumberOfFields(PARAM)) !=0)) {
	  int *idtmp=NULL;
	  if(!mxIsStruct(PARAM)){
	    mexErrMsgTxt("PARAM must be a structure !");
	  }


	  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	  /* Integer parameters                                             */
	  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	  /* level of messages output by the solver */
	  if((mxtmp=mxGetField(PARAM,0,"msglev")) != NULL){
	    rdtmp=mxGetPr(mxtmp);
	    if((*rdtmp != 0) && (*rdtmp != 1) && (*rdtmp != 2) && (*rdtmp != 3)){
	      sprintf(errmsg,"'msglev' parameter must be only:\n\t0 - no output,\n\t1 - error messages only),\n\t2 - normal output,\n\t3 - full output [default]");
	      mexErrMsgTxt(errmsg);
	    }
	    lpxIntParam[0]=*rdtmp;
	  }

	    /* scaling option */
	    if((mxtmp=mxGetField(PARAM,0,"scale")) != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1) && (*rdtmp != 2)){
				sprintf(errmsg,"'scale' parameter must be only:\n\t0 - no scaling,\n\t1 - equilibration scaling,\n\t2 - geometric mean scaling");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[1]=*rdtmp;
	    }

	    /* Dual dimplex option */
	    if((mxtmp=mxGetField(PARAM,0,"dual")) != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1)){
				sprintf(errmsg,"'dual' parameter must be only:\n\t0 - do not use the dual simplex [default],\n\t1 - use dual simplex");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[2]=*rdtmp;
	    }

	    /* pricing option */
	    if((mxtmp=mxGetField(PARAM,0,"price"))  != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1)){
				sprintf(errmsg,"'price' parameter must be only:\n\t0 - textbook pricing,\n\t1 - steepest edge pricing [default]");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[3]=*rdtmp;
	    }

	    /* solution rounding option */
	    if((mxtmp=mxGetField(PARAM,0,"round")) != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1)){
				sprintf(errmsg,"'round' parameter must be only:\n\t0 - report all primal and dual values [default],\n\t1 - replace tiny primal and dual values by exact zero");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[4]=*rdtmp;
	    }

	    /* simplex iterations limit */
	    if((mxtmp=mxGetField(PARAM,0,"itlim")) != NULL){
	      rdtmp=mxGetPr(mxtmp); lpxIntParam[5]=*rdtmp; }

	    /* Simplex iterations count */
	    if((mxtmp=mxGetField(PARAM,0,"itcnt")) != NULL){
	      rdtmp=mxGetPr(mxtmp); lpxIntParam[6]=*rdtmp; }

	    /* Output frequency, in iterations */
	    if((mxtmp=mxGetField(PARAM,0,"outfrq")) != NULL){
	      rdtmp=mxGetPr(mxtmp); lpxIntParam[7]=*rdtmp; }

	    /* Branching heuristic option  */
	    if((mxtmp=mxGetField(PARAM,0,"branch")) != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1) && (*rdtmp != 2)){
				sprintf(errmsg,"'branch' parameter must be only (for MIP only):\n\t0 - branch on the first variable,\n\t1 - branch on the last variable,\n\t2 - branch using a heuristic by Driebeck and Tomlin [default]");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[14]=*rdtmp;
	    }

	    /* Backtracking heuristic option  */
	    if((mxtmp=mxGetField(PARAM,0,"btrack")) != NULL){
	      rdtmp=mxGetPr(mxtmp);
	      if((*rdtmp != 0) && (*rdtmp != 1) && (*rdtmp != 2)){
				sprintf(errmsg,"'btrack' parameter must be only (for MIP only):\n\t0 - depth first search,\n\t1 - breadth first search,\n\t2 - backtrack using the best projection heuristic");
				mexErrMsgTxt(errmsg);
	      }
	      lpxIntParam[15]=*rdtmp;
	    }

      if((mxtmp=mxGetField(PARAM,0,"presol")) != NULL){
         rdtmp=mxGetPr(mxtmp);
         if((*rdtmp != 0) && (*rdtmp != 1)){
            sprintf(errmsg,"'presol' parameter must be only:\n\t0 - LP presolver is ***NOT*** used,\n\t1 - LP presol is used");
            mexErrMsgTxt(errmsg);
         }
         lpxIntParam[16]=*rdtmp;
      }

      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
		/* Real parameters                                                */
		/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

		/* ratio test option */
		if((mxtmp=mxGetField(PARAM,0,"relax")) != NULL){
		   rdtmp=mxGetPr(mxtmp); lpxRealParam[0]=*rdtmp; }

		/* relative tolerance used to check if the current basic solution
		   is primal feasible */
		if((mxtmp=mxGetField(PARAM,0,"tolbnd")) != NULL){
		   rdtmp=mxGetPr(mxtmp); lpxRealParam[1]=*rdtmp; }

		/* absolute tolerance used to check if the current basic solution
			is dual feasible */
		if((mxtmp=mxGetField(PARAM,0,"toldj")) != NULL){
			rdtmp=mxGetPr(mxtmp); lpxRealParam[2]=*rdtmp; }

		/* relative tolerance used to choose eligible pivotal elements of
			the simplex table in the ratio test */
		if((mxtmp=mxGetField(PARAM,0,"tolpiv")) != NULL){
			rdtmp=mxGetPr(mxtmp); lpxRealParam[3]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"objll")) != NULL){
			rdtmp=mxGetPr(mxtmp); lpxRealParam[4]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"objul")) != NULL){
		  	rdtmp=mxGetPr(mxtmp); lpxRealParam[5]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"tmlim")) != NULL){
		  rdtmp=mxGetPr(mxtmp); lpxRealParam[6]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"outdly")) != NULL){
		  rdtmp=mxGetPr(mxtmp); lpxRealParam[7]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"tolint")) != NULL){
		  rdtmp=mxGetPr(mxtmp); lpxRealParam[8]=*rdtmp; }

		if((mxtmp=mxGetField(PARAM,0,"tolobj")) != NULL){
		  rdtmp=mxGetPr(mxtmp); lpxRealParam[9]=*rdtmp; }

	}

	/* 10th Input. If the problem is a LP problem you may select which solver
	   use: RSM (Revised Simplex Method) or IPM (Interior Point Method).
	   If the problem is a MIP problem this field will be ignored.
	*/
	if ((nrhs > 9) && (mxGetM(SOLVER_IN) != 0) && (mxGetN(SOLVER_IN) != 0)) {
	  if (!mxIsNumeric(SOLVER_IN)
	      || (mxGetNumberOfDimensions(SOLVER_IN) > 2)
	      || (1 != mxGetM(SOLVER_IN))
	      || (mxGetN(SOLVER_IN) != 1)
	      || mxIsComplex(SOLVER_IN)
	      || ((tmp = mxGetPr(SOLVER_IN)) == NULL)
	      ) {
	    mexErrMsgTxt("SOLVER_IN must be a real valued scalar.");
	  } else {
	    lpsolver = tmp[0];
	  }
	}

	/* 11th Input. Saves a copy of the problem if SAVE<>0.
	 */
	if ((nrhs > 10) && (mxGetM(SAVE_IN) != 0) && (mxGetN(SAVE_IN) != 0)) {
	  if (!mxIsNumeric(SAVE_IN)
	      || (mxGetNumberOfDimensions(SAVE_IN) > 2)
	      || (1 != mxGetM(SAVE_IN))
	      || (mxGetN(SAVE_IN) != 1)
	      || mxIsComplex(SAVE_IN)
	      || ((tmp = mxGetPr(SAVE_IN)) == NULL)
	      ) {
	    mexErrMsgTxt("SAVE must be a real valued scalar.");
	  } else {
	    save_pb = (tmp[0] != 0);
	  }
	}


	if (nrhs > 11) mexWarnMsgTxt("Extra parameters ignored.");


	/* ---- Set default values ---- */
	/* CTYPE argument, default: upper bound */
	if (ctype == NULL) {
	  int i;
	  ctype = mxCalloc(mrowsA, sizeof (char));

	  for (i = 0; i < mrowsA; i++) ctype[i] = 'U';
	}
	
	/*LB argument, default: Free */
	freeLB=(int *) mxCalloc(mrowsc, sizeof(int));
	freeUB=(int *) mxCalloc(mrowsc, sizeof(int));
	if (lb == NULL) {
	  int i;
	  lb = mxCalloc(mrowsc, sizeof (double));
	  for (i = 0; i < mrowsc; i++){
	  	 lb[i] = -mxGetInf();
	  	 freeLB[i]=1;
	  }
	}else{
		int i;
		for(i=0;i<mrowsc;i++){
			if(lb[i]==-mxGetInf()) freeLB[i]=1;
			else freeLB[i]=0;
		}
	}
	/*UB argument, default: Free */
	if (ub == NULL) {
	  int i;
	  ub = mxCalloc(mrowsc, sizeof (double));
	  for (i = 0; i < mrowsc; i++){
	  	ub[i] = mxGetInf();
	  	freeUB[i]=1;
	  }
	}else{
		int i;
		for(i=0;i<mrowsc;i++){
			if(ub[i]==mxGetInf()) freeUB[i]=1;
			else freeUB[i]=0;
		}
	}
	/*VARTYPE argument, default: continuous */
	if (vartype == NULL) {
	  int i;
	  vartype2 = mxCalloc(mrowsc, sizeof (int));

	  for (i = 0; i < mrowsc; i++) vartype2[i] = LPX_CV;
	}

	extranames=mxCalloc(4,sizeof(*extranames));
	extranames[0]="lambda";
	extranames[1]="redcosts";
	extranames[2]="time";
	extranames[3]="memory";
	/* Create a matrices for the return arguments
	 */
	XMIN_OUT   = mxCreateDoubleMatrix(mrowsc, 1, mxREAL);
	FMIN_OUT   = mxCreateDoubleMatrix(1, 1, mxREAL);
	STATUS_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	EXTRA_OUT  = mxCreateStructMatrix(1, 1, 4, extranames);
	mxlambda   = mxCreateDoubleMatrix(mrowsA, 1, mxREAL);
	mxredcosts  = mxCreateDoubleMatrix(mrowsc, 1, mxREAL);
	mxtime     = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxmem      = mxCreateDoubleMatrix(1, 1, mxREAL);
	

	/* Assign pointers to the output parameters
	 */
	xmin   = mxGetPr(XMIN_OUT);
	fmin   = mxGetPr(FMIN_OUT);
	status = mxGetPr(STATUS_OUT);
	lambda = mxGetPr(mxlambda);
	redcosts= mxGetPr(mxredcosts);
	time   = mxGetPr(mxtime);
	mem    = mxGetPr(mxmem);

	jmpret = setjmp( mark );
   if (jmpret==0){
	error=glpk(sense,mrowsc,mrowsA,c,nz,rn,cn,a,b,ctype,
		    freeLB,lb,freeUB,ub,vartype2,isMIP,lpsolver,save_pb,
		    xmin,fmin,status,lambda,redcosts, time,mem);
	}

	mxSetField(EXTRA_OUT,0,extranames[0],mxlambda);
	mxSetField(EXTRA_OUT,0,extranames[1],mxredcosts);
	mxSetField(EXTRA_OUT,0,extranames[2],mxtime);
	mxSetField(EXTRA_OUT,0,extranames[3],mxmem);

	mxFree(extranames);
	mxFree(rn);
	mxFree(cn);
	mxFree(a);
	mxFree(freeLB);
	mxFree(freeUB);

	return;
}

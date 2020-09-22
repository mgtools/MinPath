//  GlpkSolver.java  (Java Native Interface for GLPK)

// ---------------------------------------------------------------------
// Copyright (C) 2003 Andrew Makhorin <mao@mai2.rcnet.ru>, Department
// for Applied Informatics, Moscow Aviation Institute, Moscow, Russia.
// All rights reserved.
//
// Author: Yuri Victorovich, Software Engineer, yuri@gjt.org.
// Author: Chris Rosebrugh, cpr@pobox.com
//         - added user-installable print and fault hooks.
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
// ---------------------------------------------------------------------

package org.gnu.glpk;

public class GlpkSolver
{

	private int lp = 0;

	static private GlpkHookIFC	m_hook=null;

	public final static int    LPX_LP         = 100;
	public final static int    LPX_MIP        = 101;
	public final static int    LPX_FR         = 110;
	public final static int    LPX_LO         = 111;
	public final static int    LPX_UP         = 112;
	public final static int    LPX_DB         = 113;
	public final static int    LPX_FX         = 114;
	public final static int    LPX_MIN        = 120;
	public final static int    LPX_MAX        = 121;
	public final static int    LPX_B_UNDEF    = 130;
	public final static int    LPX_B_VALID    = 131;
	public final static int    LPX_P_UNDEF    = 132;
	public final static int    LPX_P_FEAS     = 133;
	public final static int    LPX_P_INFEAS   = 134;
	public final static int    LPX_P_NOFEAS   = 135;
	public final static int    LPX_D_UNDEF    = 136;
	public final static int    LPX_D_FEAS     = 137;
	public final static int    LPX_D_INFEAS   = 138;
	public final static int    LPX_D_NOFEAS   = 139;
	public final static int    LPX_BS         = 140;
	public final static int    LPX_NL         = 141;
	public final static int    LPX_NU         = 142;
	public final static int    LPX_NF         = 143;
	public final static int    LPX_NS         = 144;
	public final static int    LPX_T_UNDEF    = 150;
	public final static int    LPX_T_OPT      = 151;
	public final static int    LPX_CV         = 160;
	public final static int    LPX_IV         = 161;
	public final static int    LPX_I_UNDEF    = 170;
	public final static int    LPX_I_OPT      = 171;
	public final static int    LPX_I_FEAS     = 172;
	public final static int    LPX_I_NOFEAS   = 173;
	public final static int    LPX_OPT        = 180;
	public final static int    LPX_FEAS       = 181;
	public final static int    LPX_INFEAS     = 182;
	public final static int    LPX_NOFEAS     = 183;
	public final static int    LPX_UNBND      = 184;
	public final static int    LPX_UNDEF      = 185;
	public final static int    LPX_E_OK       = 200;
	public final static int    LPX_E_EMPTY    = 201;
	public final static int    LPX_E_BADB     = 202;
	public final static int    LPX_E_INFEAS   = 203;
	public final static int    LPX_E_FAULT    = 204;
	public final static int    LPX_E_OBJLL    = 205;
	public final static int    LPX_E_OBJUL    = 206;
	public final static int    LPX_E_ITLIM    = 207;
	public final static int    LPX_E_TMLIM    = 208;
	public final static int    LPX_E_NOFEAS   = 209;
	public final static int    LPX_E_INSTAB   = 210;
	public final static int    LPX_E_SING     = 211;
	public final static int    LPX_E_NOCONV   = 212;
	public final static int    LPX_K_MSGLEV   = 300;
	public final static int    LPX_K_SCALE    = 301;
	public final static int    LPX_K_DUAL     = 302;
	public final static int    LPX_K_PRICE    = 303;
	public final static int    LPX_K_RELAX    = 304;
	public final static int    LPX_K_TOLBND   = 305;
	public final static int    LPX_K_TOLDJ    = 306;
	public final static int    LPX_K_TOLPIV   = 307;
	public final static int    LPX_K_ROUND    = 308;
	public final static int    LPX_K_OBJLL    = 309;
	public final static int    LPX_K_OBJUL    = 310;
	public final static int    LPX_K_ITLIM    = 311;
	public final static int    LPX_K_ITCNT    = 312;
	public final static int    LPX_K_TMLIM    = 313;
	public final static int    LPX_K_OUTFRQ   = 314;
	public final static int    LPX_K_OUTDLY   = 315;
	public final static int    LPX_K_BRANCH   = 316;
	public final static int    LPX_K_BTRACK   = 317;
	public final static int    LPX_K_TOLINT   = 318;
	public final static int    LPX_K_TOLOBJ   = 319;
	public final static int    LPX_K_MPSINFO  = 320;
	public final static int    LPX_K_MPSOBJ   = 321;
	public final static int    LPX_K_MPSORIG  = 322;
	public final static int    LPX_K_MPSWIDE  = 323;
	public final static int    LPX_K_MPSFREE  = 324;
	public final static int    LPX_K_MPSSKIP  = 325;
	public final static int    LPX_K_LPTORIG  = 326;
  
	static {
		System.loadLibrary("glpk_jni");
	}

	public GlpkSolver() {
	}

	protected native void finalize();
	public native void enablePrints(boolean enable);

	/**
	 * Install a callback for handling printing and fault messages
	 */
	public void setHook(GlpkHookIFC hook) {
		m_hook = hook;
	}

	private void faultHook(String s) {
		if (m_hook != null)
			m_hook.fault(s);
	}

	private void printHook(String s) {
		if (m_hook != null)
			m_hook.print(s);
	}

	// problem creating and modifying routines
	public native void addRows(int nrs);
	public native void addCols(int ncs);
	public native static boolean checkName(String name);
	public native void setProbName(String name);
	public native void setRowName(int i, String name);
	public native void setColName(int j, String name);
	public native void setRowBnds(int i, int typx, double lb, double ub);
	public native void setColBnds(int i, int typx, double lb, double ub);
	public native void setObjName(String name);
	public native void setObjDir(int dir);
	public native void setObjC0(double c0);
	public native void setRowCoef(int i, double coef);
	public native void setColCoef(int j, double coef);
	public native void loadMat(GlpkSolverLoader loader);
	public native void loadMat3(int nz, int rn[], int cn[], double a[]);
	public native void setMatRow(int i, int len, int ndx[], double val[]);
	public native void setMatCol(int j, int len, int ndx[], double val[]);
	public native void unmarkAll();
	public native void markRow(int i, int mark);
	public native void markCol(int j, int mark);
	public native void clearMat();
	public native void delItems();

	// Problem querying routines
	public native int getNumRows();
	public native int getNumCols();
	public native int getNumNz();
	public native String getProbName();
	public native String getRowName(int i);
	public native String getColName(int j);
	public native void getRowBnds(int i, GlpkSolverBnds bnds);
	public native void getColBnds(int j, GlpkSolverBnds bnds);
	public native String getObjName();
	public native int getObjDir();
	public native double getObjC0();
	public native double getRowCoef(int i);
public native double getColCoef(int j);
public native int getMatRow(int i, int ndx[], double val[]);
public native int getMatCol(int j, int ndx[], double val[]);
public native int getRowMark(int i);
public native int getColMark(int j);
// Problem scaling routines
public native void scaleProb();
public native void unscaleProb();
// Basis constructing routines
public native void stdBasis();
public native void advBasis();
public native void setRowStat(int i, int stat);
public native void setColStat(int j, int stat);
// Simplex method routines
public native int warmUp();
public native int simplex();
// Basic solution querying routines
public native int getStatus();
public native int getPrimStat();
public native int getDualStat();
public native void getRowInfo(int i, GlpkSolverInfo info);
public native void getColInfo(int j, GlpkSolverInfo info);
public native double getObjVal();
public native void checkKkt(int scaled, GlpkSolverKktConditions kkt);
// Simplex table routines
public native int evalTabRow(int k, int ndx[], double val[]);
public native int evalTabCol(int k, int ndx[], double val[]);
public native int transformRow(int len, int ndx[], double val[]);
public native int transformCol(int len, int ndx[], double val[]);
public native int primRatioTest(int len, int ndx[], double val[],
  int how, double tol);
public native int dualRatioTest(int len, int ndx[], double val[],
  int how, double tol);
// Interior point method routines
public native int interior();
public native int getIpsStat();
public native void getIpsRow(int i, GlpkSolverPoint point);
public native void getIpsCol(int i, GlpkSolverPoint point);
public native double getIpsObj();
// MIP routines
public native void setClss(int clss);
public native int getClss();
public native void setColKind(int j, int kind);
public native int getColKind(int j);
public native int getNumInt();
public native int getNumBin();
public native int integer();
public native int getMIPStat();
public native double getMIPRow(int i);
public native double getMIPCol(int j);
public native double getMIPObj();
// Control parameter routines
public native void resetParms();
public native void setIntParm(int param, int val);
public native int getIntParm(int param);
public native void setRealParm(int param, double val);
public native double getRealParm(int param);
// Utility routines
public native static GlpkSolver readMps(String fname);
public native static GlpkSolver readLpt(String fname);
public native static GlpkSolver readModel(String model, String data,
  String output);
public native boolean writeMps(String fname);
public native boolean writeLpt(String fname);
public native int printProb(String fname);
public native boolean readBas(String fname);
public native boolean writeBas(String fname);
public native boolean printSol(String fname);
public native boolean printIps(String fname);
public native boolean printMip(String fname);

// Functions missing from documentation
public native void reallocProb(int m_max, int n_max);
public native int reduceForm(int len, int ndx[], double val[], double work[]);
public native int primOpt();
public native int primArt();
public native int dualOpt();
public native int mixedGomory(int kind[], int len, int ndx[],
      double val[], double work[]);
public native double evalRedCost(int len, int ndx[], double val[]);
public native double evalActivity(int len, int ndx[], double val[]);
public native void estimObjChg(int k, double dn_val, double up_val,
      int ndx[], double val[], GlpkSolverObjChg chg);
}

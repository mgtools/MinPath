//  GlpkMagicSquare.java  (Java Native Interface demo file)

// ---------------------------------------------------------------------
// Copyright (C) 2003 Andrew Makhorin <mao@mai2.rcnet.ru>, Department
// for Applied Informatics, Moscow Aviation Institute, Moscow, Russia.
// All rights reserved.
//
// Author: Yuri Victorovich, Software Engineer, yuri@gjt.org.
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

import org.gnu.glpk.*;

//
// We solve magic square problem here with the given size
//

public class GlpkMagicSquare implements GlpkHookIFC {
	static private GlpkSolver solver;

	public static void main(String Args[]) {
		try {
			solver = new GlpkSolver();
		} catch (Error e) {
			System.err.println("Sample: Can't instantiate solver:");
			System.err.println("Sample:  ** " + e.getClass().getName() + ": " + e.getMessage());
			System.err.println("Sample:  ** java.library.path: " + System.getProperty("java.library.path"));
			System.err.println("Sample: Probably you don't have GLPK JNI properly installed.");
			return;
		}

		new GlpkMagicSquare();
	}

	public GlpkMagicSquare() {

		solver.setHook(this);
		
		// since we've hooked the messages we don't need stdout/err printing
		solver.enablePrints(false);

		solver.setClss(GlpkSolver.LPX_MIP);
		solver.setObjDir(GlpkSolver.LPX_MIN);
		// size, number of cells and summ of each row, col and diagonal
		int n = 3; // n=5 already takes about 5 minutes
		int N = n*n, sum = (N+1)*n/2;
		int col = 1, row = 1;          // current row/column variables
		int vars[] = new int[1+N+1];   // constraint arrays with maximal needed sizes
		double vals[] = new double[1+N+1];
		int res, stat;

		// adding binary variables: if I is in ixj cell
		solver.addCols(N*n*n);
		for (int c = 0; c < N*n*n; c++) {
			solver.setColKind(col, GlpkSolver.LPX_IV);
			solver.setColBnds(col, GlpkSolver.LPX_DB, 0.0, 1.0);
			col++;
		}

		// adding solution variables (actual solution integers)
		solver.addCols(n*n);
		for (int cell = 0; cell < n*n; cell++) {
			solver.setColKind(col, GlpkSolver.LPX_IV);
			solver.setColBnds(col, GlpkSolver.LPX_DB, 1.0, (double)N);
			col++;
		}

		// adding per-cell uniqueness constraints
		solver.addRows(n*n);
		for (int r = 0; r < n; r++) {
			for (int c = 0; c < n; c++) {
				for (int i = 1; i <= N; i++) {
					vars[i] = 1 + (i-1)*n*n + r*n + c;
					vals[i] = 1.0;
				}

				solver.setMatRow(row, N, vars, vals);
				solver.setRowBnds(row, GlpkSolver.LPX_FX, 1.0, 1.0);
				row++;
			}
		}

		// adding per-integer uniqueness constraints
		solver.addRows(N);
		for (int i = 1; i <= N; i++) {
			for (int r = 0; r < n; r++) {
				for (int c = 0; c < n; c++) {
					vars[1+r*n+c] = 1 + (i-1)*n*n + r*n + c;
					vals[1+r*n+c] = 1.0;
				}
			}
			solver.setMatRow(row, N, vars, vals);
			solver.setRowBnds(row, GlpkSolver.LPX_FX, 1.0, 1.0);
			row++;
		}

		// adding boolean-integer binding constraints
		solver.addRows(n*n);
		for (int r = 0; r < n; r++) {
			for (int c = 0; c < n; c++) {
				for (int i = 1; i <= N; i++) {
					vars[i] = 1 + (i-1)*n*n + r*n + c;
					vals[i] = (double)i;
				}
				vars[N+1] = 1 + N*n*n + r*n + c;
				vals[N+1] = -1.0;
				solver.setMatRow(row, N+1, vars, vals);
				solver.setRowBnds(row, GlpkSolver.LPX_FX, 0.0, 0.0);
				row++;
			}
		}

		// adding row constraints
		solver.addRows(n);
		for (int r = 0; r < n; r++) {
			for (int c = 0; c < n; c++) {
				vars[1 + c] = 1 + N*n*n + r*n + c;
				vals[1 + c] = 1.0;
			}
			solver.setMatRow(row, n, vars, vals);
			solver.setRowBnds(row, GlpkSolver.LPX_FX, (double)sum, (double)sum);
			row++;
		}

		// adding column constraints
		solver.addRows(n);
		for (int c = 0; c < n; c++) {
			for (int r = 0; r < n; r++) {
				vars[1 + r] = 1 + N*n*n + r*n + c;
				vals[1 + r] = 1.0;
			}
			solver.setMatRow(row, n, vars, vals);
			solver.setRowBnds(row, GlpkSolver.LPX_FX, (double)sum, (double)sum);
			row++;
		}

		{ // adding diagonal constraints
			solver.addRows(2);
			// main diagonal
			for (int d = 0; d < n; d++) {
				int r = d, c = d;
				vars[1 + d] = 1 + N*n*n + r*n + c;
				vals[1 + d] = 1.0;
			}
			solver.setMatRow(row, n, vars, vals);
			solver.setRowBnds(row, GlpkSolver.LPX_FX, (double)sum, (double)sum);
			row++;
			// secondary diagonal
			for (int d = 0; d < n; d++) {
				int r = d, c = n - d - 1;
				vars[1 + d] = 1 + N*n*n + r*n + c;
				vals[1 + d] = 1.0;
			}
			solver.setMatRow(row, n, vars, vals);
			solver.setRowBnds(row, GlpkSolver.LPX_FX, (double)sum, (double)sum);
			row++;
		} // diagonals

		// objective picks solution eith the minimal (0,0),(0,1) combination
		solver.setColCoef(1 + N*n*n, (double)N);
		solver.setColCoef(1 + N*n*n + 1, 1.0);

		System.out.println("Sample: variables=" + col + " rows=" + row);

		// solving
		res = solver.simplex();
		if (res != GlpkSolver.LPX_E_OK ||
			(solver.getStatus() != GlpkSolver.LPX_OPT &&
			 solver.getStatus() != GlpkSolver.LPX_FEAS)) {
			System.err.println("Sample: simplex() failed");
			return;
		}

		res = solver.integer();
		if (res != GlpkSolver.LPX_E_OK ||
			(solver.getStatus() != GlpkSolver.LPX_OPT &&
			 solver.getStatus() != GlpkSolver.LPX_FEAS)) {
			System.err.println("Sample: integer() failed");
			return;
		}

		// printing solution
		System.out.print("Sample: Solution found: Magic square with n=" + n + " is");
		{
			int maxDigits = (new String("") + N).length();
			for (int r = 0; r < n; r++) {
				System.out.print("\nSample: ");
				for (int c = 0; c < n; c++) {
					int i = (int)solver.getMIPCol(1 + N*n*n + r*n + c);
					String s = "" + i;
					while (s.length() < maxDigits)
						s = " " + s;
					System.out.print(" " + s);
				}
			}

			System.out.println();
		}
	}

	// implements GlpkHookIFC
	public void fault(String s) {
		System.out.println("GlpkFAULT: " + s);
		System.exit(-1);
	}

	public void print(String s) {
		System.out.println("GlpkMsg: " + s);
	}
}

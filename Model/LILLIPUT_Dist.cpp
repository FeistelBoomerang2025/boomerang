#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

const int P[16] = { 13, 9, 14, 8, 10, 11, 12, 15, 4, 5, 3, 1, 2, 6, 0, 7 };
const int P_inv[16] = { 14, 11, 12, 10, 8, 9, 13, 15, 3, 1, 4, 5, 6, 0, 2, 7 };
const int ROUND_m = 4;
const int ROUND_0 = 3;
const int ROUND_1 = 6;
const int ROUND = ROUND_0 + ROUND_m + ROUND_1;
const int state = 16;

// suite for XOR function: a + b = c
static void XOR_deter(GRBModel& model, GRBVar a, GRBVar b, GRBVar c) {
	model.addConstr(c >= a);
	model.addConstr(c >= b);
	model.addConstr(c <= a + b);
}

static void XOR_prob(GRBModel& model, GRBVar a, GRBVar b, GRBVar c) {
	model.addConstr(c >= a - b);
	model.addConstr(c >= b - a);
	model.addConstr(c <= a + b);
}

static void XOR_prob_cancel(GRBModel& model, GRBVar a, GRBVar b, GRBVar c) {
	model.addConstr(a <= b + c);
	model.addConstr(b <= a + c);
	model.addConstr(c <= a + b);
	model.addConstr(a + b + c <= 2);
}

static void initial_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z, vector<vector<vector<GRBVar>>>& W) {
	// state[round][state][j]
	for (int i = 0; i < state; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				if (i < 8) Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Z[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				W[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				if (i < 8) model.addConstr(W[r][i][j] == X[r][i][j]);
			}
		}
	}
}

static void initial_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z, vector<vector<vector<GRBVar>>>& W) {
	// state[round][state][j]
	for (int i = 0; i < state; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				if (i < state / 2) Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Z[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				W[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				if (i < state / 2) model.addConstr(W[r][i][j] == X[r][i][j]);
			}
		}
	}
}

static void SubBytes(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y) {
	// X[0-7] -S-> Y[0-7]
	for (int r = 0; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][i][0] == Y[r][i][0]);
		}
	}
}

static void NLXor_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Y[0-7] + X[15-8] = Z[15-8]
	
	// Prob Propagate in E0
	for (int r = 0; r < ROUND_0; r++) {
		for (int i = 0; i < state / 2; i++) {
			//XOR_prob(model, X[r][15 - i][0], Y[r][i][0], Z[r][15 - i][0]);
			XOR_prob_cancel(model, X[r][15 - i][0], Y[r][i][0], Z[r][15 - i][0]);
		}
	}
	// Deter Propagate in Em+E1
	for (int r = ROUND_0; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_deter(model, X[r][15 - i][0], Y[r][i][0], Z[r][15 - i][0]);
		}
	}
}

static void NLXor_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Y[0-7] + Z[15-8] = X[15-8]
	
	// Deter Propagate in E0+Em
	for (int r = 0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_deter(model, Y[r][i][0], Z[r][15 - i][0], X[r][15 - i][0]);
		}
	}
	// Prob Propagate in E1
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			//XOR_prob(model, Y[r][i][0], Z[r][15 - i][0], X[r][15 - i][0]);
			XOR_prob_cancel(model, Y[r][i][0], Z[r][15 - i][0], X[r][15 - i][0]);
		}
	}
}

static void Linear_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Z, vector<vector<vector<GRBVar>>>& W) {
	// Prob propagate in E0
	for (int r = 0; r < ROUND_0; r++) {
		GRBVar tmp[8];
		for (int j = 0; j < state / 2; j++) tmp[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(tmp[0] == Z[r][15][0]);
		model.addConstr(tmp[7] == W[r][15][0]);
		model.addConstr(W[r][8][0] == Z[r][8][0]);
		
		for (int i = 1; i < state / 2; i++) {
			// tmp[0-6] + X[1-7] = tmp[1-7]
			
			//XOR_prob(model, X[r][i][0], tmp[i - 1], tmp[i]);
			XOR_prob_cancel(model, X[r][i][0], tmp[i - 1], tmp[i]);

			// Z[9-14] + X[7] = W[9-14]
			if (i < state / 2 - 1) { 
				//XOR_prob(model, X[r][7][0], Z[r][i + 8][0], W[r][i + 8][0]);
				XOR_prob_cancel(model, X[r][7][0], Z[r][i + 8][0], W[r][i + 8][0]);
			}
		}
	}

	// Deter propagate in Em+E1
	for (int r = ROUND_0; r < ROUND; r++) {
		GRBVar tmp[8];
		for (int j = 0; j < state / 2; j++) tmp[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(tmp[0] == Z[r][15][0]);
		model.addConstr(tmp[7] == W[r][15][0]);
		model.addConstr(W[r][8][0] == Z[r][8][0]);

		for (int i = 1; i < state / 2; i++) {
			// tmp[0-6] + X[1-7] = tmp[1-7]
			XOR_deter(model, X[r][i][0], tmp[i - 1], tmp[i]);

			// Z[9-14] + X[7] = W[9-14]
			if (i < state / 2 - 1) {
				XOR_deter(model, X[r][7][0], Z[r][i + 8][0], W[r][i + 8][0]);
			}
		}
	}
}

static void Linear_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Z, vector<vector<vector<GRBVar>>>& W) {
	// Deter propagate in E0+Em
	for (int r = 0; r < ROUND_0 + ROUND_m; r++) {
		GRBVar tmp[8];
		for (int j = 0; j < state / 2; j++) tmp[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(tmp[0] == Z[r][15][0]);
		model.addConstr(tmp[7] == W[r][15][0]);
		model.addConstr(W[r][8][0] == Z[r][8][0]);

		for (int i = 1; i < state / 2; i++) {
			// X[1-7] + tmp[1-7] = tmp[0-6]
			XOR_deter(model, X[r][i][0], tmp[i], tmp[i - 1]);

			// X[7] + W[9-14] = Z[9-14]
			if (i < state / 2 - 1) {
				XOR_deter(model, X[r][7][0], W[r][i + 8][0], Z[r][i + 8][0]);
			}
		}
	}

	// Prob propagate in E1
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		GRBVar tmp[8];
		for (int j = 0; j < state / 2; j++) tmp[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(tmp[0] == Z[r][15][0]);
		model.addConstr(tmp[7] == W[r][15][0]);
		model.addConstr(W[r][8][0] == Z[r][8][0]);

		for (int i = 1; i < state / 2; i++) {
			// X[1-7] + tmp[1-7] = tmp[0-6]
			//XOR_prob(model, X[r][i][0], tmp[i], tmp[i - 1]);
			XOR_prob_cancel(model, X[r][i][0], tmp[i], tmp[i - 1]);

			// X[7] + W[9-14] = Z[9-14]
			if (i < state / 2 - 1) {
				//XOR_prob(model, X[r][7][0], W[r][i + 8][0], Z[r][i + 8][0]);
				XOR_prob_cancel(model, X[r][7][0], W[r][i + 8][0], Z[r][i + 8][0]);
			}
		}
	}
}

static void Perm(GRBModel& model, vector<vector<vector<GRBVar>>>& W, vector<vector<vector<GRBVar>>>& X) {
	for (int r = 0; r < ROUND - 1; r++) {
		for (int i = 0; i < state; i++) {
			for (int j = 0; j < 2; j++) {
				model.addConstr(X[r + 1][P[i]][j] == W[r][i][j]);
			}
		}
	}
}

static void Em_Optimize(GRBModel& model, vector<vector<vector<GRBVar>>>& Sbox,
	vector<vector<vector<GRBVar>>>& Up_State_X, vector<vector<vector<GRBVar>>>& Up_State_Y, vector<vector<vector<GRBVar>>>& Up_State_Z,
	vector<vector<vector<GRBVar>>>& Lo_State_X, vector<vector<vector<GRBVar>>>& Lo_State_Y, vector<vector<vector<GRBVar>>>& Lo_State_Z) {

	// common active
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar common[2] = { Up_State_X[r][i][0], Lo_State_X[r][i][0] };
			model.addGenConstrAnd(Sbox[r][i][0], common, 2);
		}
	}

	//// constraint S-box in E0
	//for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
	//	for (int i = 0; i < state / 2; i++) {
	//		model.addGenConstrIndicator(Sbox[r][i][0], 1, Up_State_X[r][2 * i][1] == 1);
	//	}
	//}
	//for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
	//	for (int i = 0; i < state / 2; i++) {
	//		GRBVar conn1[2] = { Up_State_Z[r][2 * i + 1][1], Up_State_X[r][2 * i][0] };
	//		GRBVar conn2[2] = { Up_State_Z[r][2 * i + 1][1], Up_State_X[r][2 * i + 1][0] };
	//		GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	//		GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	//		model.addGenConstrAnd(And1, conn1, 2);
	//		model.addGenConstrAnd(And2, conn2, 2);
	//		model.addGenConstrIndicator(And1, 1, Up_State_X[r][2 * i][1] == 1);
	//		model.addGenConstrIndicator(And2, 1, Up_State_X[r][2 * i + 1][1] == 1);

	//		model.addGenConstrIndicator(And1, 1, Sbox[r][i][1] == 1);
	//	}
	//}

	//// constraint S-box in E1
	//for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
	//	for (int i = 0; i < state / 2; i++) {
	//		model.addGenConstrIndicator(Sbox[r][i][0], 1, Lo_State_Z[r][2 * i][1] == 1);
	//	}
	//}
	//for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
	//	for (int i = 0; i < state / 2; i++) {
	//		GRBVar conn1[2] = { Lo_State_X[r][2 * i + 1][1], Lo_State_Z[r][2 * i][0] };
	//		GRBVar conn2[2] = { Lo_State_X[r][2 * i + 1][1], Lo_State_Z[r][2 * i + 1][0] };
	//		GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	//		GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	//		model.addGenConstrAnd(And1, conn1, 2);
	//		model.addGenConstrAnd(And2, conn2, 2);
	//		model.addGenConstrIndicator(And1, 1, Lo_State_Z[r][2 * i][1] == 1);
	//		model.addGenConstrIndicator(And2, 1, Lo_State_Z[r][2 * i + 1][1] == 1);

	//		model.addGenConstrIndicator(And1, 1, Sbox[r][i][2] == 1);
	//	}
	//}
}

static void ResultOutput(GRBModel& model, GRBLinExpr obj0, GRBLinExpr objm, GRBLinExpr obj1,
	vector<vector<vector<GRBVar>>>& Up_State_X, vector<vector<vector<GRBVar>>>& Up_State_Y, vector<vector<vector<GRBVar>>>& Up_State_Z, vector<vector<vector<GRBVar>>>& Up_State_W,
	vector<vector<vector<GRBVar>>>& Lo_State_X, vector<vector<vector<GRBVar>>>& Lo_State_Y, vector<vector<vector<GRBVar>>>& Lo_State_Z, vector<vector<vector<GRBVar>>>& Lo_State_W) {

	string filename = "LILLIPUT_Dist_R" + to_string(ROUND_0 + ROUND_m + ROUND_1) + "_U" + to_string(ROUND_0) + "_M" + to_string(ROUND_m) + "_L" + to_string(ROUND_1) + ".txt";
	ofstream out(filename, ios::out | ios::trunc);

	out << "Round for E0 Em E1 : " << ROUND_0 << ' ' << ROUND_m << ' ' << ROUND_1 << endl;
	out << "Active in E0 Em E1: " << obj0.getValue() << ' ' << objm.getValue() << ' ' << obj1.getValue() << endl;

	out << "E0" << endl;
	for (int r = 0; r < ROUND_0; r++) {
		out << r << ": " << endl;
		// X[15-0]
		for (int i = state - 1; i >= 0; i--) {
			out << Up_State_X[r][i][0].get(GRB_DoubleAttr_X) ;
		}out << endl;
		// -------- Y[7-0]
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << ' ' << ' ';
		}
		for (int i = state / 2 - 1; i >= 0; i--) {
			out << Up_State_Y[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// Z[15-8] 
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << Up_State_Z[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// W[15-0]
		for (int i = state - 1; i >= 0; i--) {
			out << Up_State_W[r][i][0].get(GRB_DoubleAttr_X);
		}out << endl;
	}

	out << "Em" << endl;
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		out << r << ": " << endl;
		// X[15-0]
		if (r == ROUND_0) {
			for (int i = state - 1; i >= 0; i--) {
				out << Up_State_X[r][i][0].get(GRB_DoubleAttr_X);
			}out << endl;
		}
		else {
			for (int i = state - 1; i >= 0; i--) {
				out << Up_State_X[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_X[r][i][0].get(GRB_DoubleAttr_X) << ' ';
			}out << endl;
		}
		// -------- Y[7-0]
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << "  " << ' ';
		}
		for (int i = state / 2 - 1; i >= 0; i--) {
			out << Up_State_Y[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_Y[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// Z[15-8]
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << Up_State_Z[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_Z[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// W[15-0]
		if (r == ROUND_0 + ROUND_m - 1) {
			for (int i = state - 1; i >= 0; i--) {
				out << Lo_State_W[r][i][0].get(GRB_DoubleAttr_X);
			}out << endl;
		}
		else {
			for (int i = state - 1; i >= 0; i--) {
				out << Up_State_W[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_W[r][i][0].get(GRB_DoubleAttr_X) << ' ';
			}out << endl;
		}
	}

	out << "E1" << endl;
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		out << r << ": " << endl;
		// X[15-0]
		for (int i = state - 1; i >= 0; i--) {
			out << Lo_State_X[r][i][0].get(GRB_DoubleAttr_X);
		}out << endl;
		// -------- Y[7-0]
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << ' ' << ' ';
		}
		for (int i = state / 2 - 1; i >= 0; i--) {
			out << Lo_State_Y[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// Z[15-8] 
		for (int i = state - 1; i > state / 2 - 1; i--) {
			out << Lo_State_Z[r][i][0].get(GRB_DoubleAttr_X) << ' ';
		}out << endl;
		// W[15-0]
		for (int i = state - 1; i >= 0; i--) {
			out << Lo_State_W[r][i][0].get(GRB_DoubleAttr_X);
		}out << endl;
	}
}

int main() {

	try {
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "mip1.log");
		env.start();
		GRBModel model = GRBModel(env);

		// State Variable Declaration X/Y/Z [round][i][j] 
		// difference [j0] = 0: zero difference
		// difference [j1] = 1: fixed nonzero difference
		// [j1] -> identify conn S_box
		vector<vector<vector<GRBVar>>> Up_State_X(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_Y(ROUND, vector<vector<GRBVar>>(state / 2, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_Z(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_W(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		initial_Up(model, Up_State_X, Up_State_Y, Up_State_Z, Up_State_W);

		vector<vector<vector<GRBVar>>> Lo_State_X(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_Y(ROUND, vector<vector<GRBVar>>(state / 2, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_Z(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_W(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		initial_Lo(model, Lo_State_X, Lo_State_Y, Lo_State_Z, Lo_State_W);

		// [0] = 1 - active Sbox
		// [1] = 1 - connected Sbox in E0
		// [2] = 1 - connected Sbox in E1
		vector<vector<vector<GRBVar>>> Sbox(ROUND, vector<vector<GRBVar>>(state / 2, vector<GRBVar>(3)));
		for (int r = 0; r < ROUND; r++) {
			for (int i = 0; i < state / 2; i++) {
				for (int j = 0; j < 3; j++) {
					Sbox[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
		}

		// Upper Trail 
		SubBytes(model, Up_State_X, Up_State_Y);
		NLXor_Up(model, Up_State_X, Up_State_Y, Up_State_Z);
		Linear_Up(model, Up_State_X, Up_State_Z, Up_State_W);
		Perm(model, Up_State_W, Up_State_X);

		// Lower Trail 
		SubBytes(model, Lo_State_X, Lo_State_Y);
		NLXor_Lo(model, Lo_State_X, Lo_State_Y, Lo_State_Z);
		Linear_Lo(model, Lo_State_X, Lo_State_Z, Lo_State_W);
		Perm(model, Lo_State_W, Lo_State_X);

		// Both Input of E0 and Output of E1 >= 1 
		GRBLinExpr exp1 = 0;
		GRBLinExpr exp2 = 0;
		for (int i = 0; i < state; i++) {
			exp1 += Up_State_X[0][i][0];
			exp2 += Lo_State_W[ROUND - 1][i][0];
		}
		model.addConstr(exp1 >= 1);
		model.addConstr(exp2 >= 1);

		GRBLinExpr exp_m1 = 0;
		GRBLinExpr exp_m2 = 0;
		for (int i = 0; i < state; i++) {
			exp_m1 += Up_State_X[ROUND_0][i][0];
			exp_m2 += Lo_State_W[ROUND_0 + ROUND_m - 1][i][0];
		}
		//model.addConstr(exp_m1 == 1);
		//model.addConstr(exp_m2 == 1);


		// Optimization

		GRBLinExpr obj_m = 0;
		GRBLinExpr obj_0 = 0;
		GRBLinExpr obj_1 = 0;
		Em_Optimize(model, Sbox, Up_State_X, Up_State_Y, Up_State_Z, Lo_State_X, Lo_State_Y, Lo_State_Z);
		for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
			for (int i = 0; i < state / 2; i++) {
				obj_m += Sbox[r][i][0];
			}
		}
		for (int r = 0; r < ROUND_0; r++) {
			for (int i = 0; i < state / 2; i++) {
				model.addConstr(Sbox[r][i][0] == Up_State_X[r][i][0]);

				obj_0 += Sbox[r][i][0];
			}
		}
		for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
			for (int i = 0; i < state / 2; i++) {
				model.addConstr(Sbox[r][i][0] == Lo_State_X[r][i][0]);
				obj_1 += Sbox[r][i][0];
			}
		}


		GRBLinExpr obj = 0;
		model.addConstr(obj_0 + obj_1 <= 11);
		obj += obj_m;
		//model.addConstr(obj_m <= 16);
		//obj += obj_0 + obj_1;
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();


		// Output Recording
		ResultOutput(model, obj_0, obj_m, obj_1, Up_State_X, Up_State_Y, Up_State_Z, Up_State_W, Lo_State_X, Lo_State_Y, Lo_State_Z, Lo_State_W);

		cout << obj_0.getValue() << ' ' << obj_m.getValue() << ' ' << obj_1.getValue() << endl;

		cout << "Prob of E0 & E1: " << (obj_0.getValue() + obj_1.getValue()) * 4 << endl;

		cout << "Position of active Sbox in Em: U";
		for (int i = 0; i < state; i++) if (Up_State_X[ROUND_0][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';
		cout << "L";
		for (int i = 0; i < state; i++) if (Lo_State_W[ROUND_0 + ROUND_m - 1][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';

		// cout << (obj_0.getValue() + obj_1.getValue()) * 4 + p_dist << endl;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}
}
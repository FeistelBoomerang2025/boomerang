#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

const int P[32] = { 31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10, 15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26 };
const int ROUND_m = 10;
const int ROUND_0 = 7;
const int ROUND_1 = 7;
const int ROUND = ROUND_0 + ROUND_m + ROUND_1;
const int state = 32;

const double p_dist = 8;

static void initial_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// state[round][state][j]
	for (int i = 0; i < state / 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				X[r][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Z[r][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Z[r][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				model.addConstr(X[r][2 * i][j] == Z[r][2 * i][j]);
			}
			X[ROUND][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			X[ROUND][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}

static void initial_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// state[round][state][j]
	for (int i = 0; i < state / 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				X[r][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				Z[r][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Z[r][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				model.addConstr(X[r][2 * i][j] == Z[r][2 * i][j]);
			}
			X[ROUND][2 * i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			X[ROUND][2 * i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}

static void SubBytes(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y) {
	for (int r = 0; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][2 * i][0] == Y[r][i][0]);
		}
	}
}

static void RightXor_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Prob Propagate in E0
	for (int r = 0; r < ROUND_0; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(Z[r][2 * i + 1][0] >= Y[r][i][0] - X[r][2 * i + 1][0]);
			model.addConstr(Z[r][2 * i + 1][0] >= X[r][2 * i + 1][0] - Y[r][i][0]);
			model.addConstr(Z[r][2 * i + 1][0] <= X[r][2 * i + 1][0] + Y[r][i][0]);
		}
	}
	// Deter Propagate in Em+E1
	for (int r = ROUND_0; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(Z[r][2 * i + 1][0] >= Y[r][i][0]);
			model.addConstr(Z[r][2 * i + 1][0] >= X[r][2 * i + 1][0]);
			model.addConstr(Z[r][2 * i + 1][0] <= X[r][2 * i + 1][0] + Y[r][i][0]);
		}
	}
}

static void RightXor_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Deter Propagate in E0+Em
	for (int r = 0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][2 * i + 1][0] >= Z[r][2 * i + 1][0]);
			model.addConstr(X[r][2 * i + 1][0] >= Y[r][i][0]);
			model.addConstr(X[r][2 * i + 1][0] <= Y[r][i][0] + Z[r][2 * i + 1][0]);
		}
	}
	// Prob Propagate in E1
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][2 * i + 1][0] >= Z[r][2 * i + 1][0] - Y[r][i][0]);
			model.addConstr(X[r][2 * i + 1][0] >= Y[r][i][0] - Z[r][2 * i + 1][0]);
			model.addConstr(X[r][2 * i + 1][0] <= Y[r][i][0] + Z[r][2 * i + 1][0]);
		}
	}
}

static void Perm(GRBModel& model, vector<vector<vector<GRBVar>>>& Z, vector<vector<vector<GRBVar>>>& X) {
	for (int r = 0; r < ROUND; r++) {
		for (int i = 0; i < state; i++) {
			for (int j = 0; j < 2; j++) {
				model.addConstr(X[r + 1][P[i]][j] == Z[r][i][j]);
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
			GRBVar common[2] = { Up_State_Z[r][2 * i][0], Lo_State_Z[r][2 * i][0] };
			model.addGenConstrAnd(Sbox[r][i][0], common, 2);
		}
	}

	// constraint S-box in E0
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addGenConstrIndicator(Sbox[r][i][0], 1, Up_State_X[r][2 * i][1] == 1);
		}
	}
	//
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar conn1[2] = { Up_State_Z[r][2 * i + 1][1], Up_State_X[r][2 * i][0] };
			GRBVar conn2[2] = { Up_State_Z[r][2 * i + 1][1], Up_State_X[r][2 * i + 1][0] };
			GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			model.addGenConstrAnd(And1, conn1, 2);
			model.addGenConstrAnd(And2, conn2, 2);
			model.addGenConstrIndicator(And1, 1, Up_State_X[r][2 * i][1] == 1);
			model.addGenConstrIndicator(And2, 1, Up_State_X[r][2 * i + 1][1] == 1);

			model.addGenConstrIndicator(And1, 1, Sbox[r][i][1] == 1);
		}
	}

	// constraint S-box in E1
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			model.addGenConstrIndicator(Sbox[r][i][0], 1, Lo_State_Z[r][2 * i][1] == 1);
		}
	}
	//
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar conn1[2] = { Lo_State_X[r][2 * i + 1][1], Lo_State_Z[r][2 * i][0] };
			GRBVar conn2[2] = { Lo_State_X[r][2 * i + 1][1], Lo_State_Z[r][2 * i + 1][0] };
			GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			model.addGenConstrAnd(And1, conn1, 2);
			model.addGenConstrAnd(And2, conn2, 2);
			model.addGenConstrIndicator(And1, 1, Lo_State_Z[r][2 * i][1] == 1);
			model.addGenConstrIndicator(And2, 1, Lo_State_Z[r][2 * i + 1][1] == 1);

			model.addGenConstrIndicator(And1, 1, Sbox[r][i][2] == 1);
		}
	}
}

static void ResultOutput(GRBModel& model, GRBLinExpr obj0, GRBLinExpr objm, GRBLinExpr obj1,
	vector<vector<vector<GRBVar>>>& Up_State_X, vector<vector<vector<GRBVar>>>& Up_State_Y, vector<vector<vector<GRBVar>>>& Up_State_Z,
	vector<vector<vector<GRBVar>>>& Lo_State_X, vector<vector<vector<GRBVar>>>& Lo_State_Y, vector<vector<vector<GRBVar>>>& Lo_State_Z) {

	string filename = "WARP_Dist_U" + to_string(ROUND_0) + "_M" + to_string(ROUND_m) + "_L" + to_string(ROUND_1) + ".txt";
	ofstream out(filename, ios::out | ios::trunc);

	out << "Round for E0 Em E1 : " << ROUND_0 << ' ' << ROUND_m << ' ' << ROUND_1 << endl;
	out << "Active in E0 Em E1: " << obj0.getValue() << ' ' << objm.getValue() << ' ' << obj1.getValue() << endl;

	out << "E0" << endl;
	for (int r = 0; r < ROUND_0; r++) {
		out << r << ": " << endl;
		for (int i = 0; i < state / 2; i++) {
			out << Up_State_X[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Up_State_X[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;

		for (int i = 0; i < state / 2; i++) {
			out << Up_State_Z[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Up_State_Z[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;
	}

	out << "Em" << endl;
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		out << r << ": " << endl;
		for (int i = 0; i < state / 2; i++) {
			out << Up_State_X[r][2 * i][0].get(GRB_DoubleAttr_X) << Lo_State_X[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Up_State_X[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << Lo_State_X[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;

		for (int i = 0; i < state / 2; i++) {
			out << Up_State_Z[r][2 * i][0].get(GRB_DoubleAttr_X) << Lo_State_Z[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Up_State_Z[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << Lo_State_Z[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;
	}

	out << "E1" << endl;
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		out << r << ": " << endl;
		for (int i = 0; i < state / 2; i++) {
			out << Lo_State_X[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Lo_State_X[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;

		for (int i = 0; i < state / 2; i++) {
			out << Lo_State_Z[r][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Lo_State_Z[r][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
		}out << endl;
	}

	out << ROUND << ": " << endl;
	for (int i = 0; i < state / 2; i++) {
		out << Lo_State_X[ROUND][2 * i][0].get(GRB_DoubleAttr_X) << ' ' << Lo_State_X[ROUND][2 * i + 1][0].get(GRB_DoubleAttr_X) << "  ";
	}out << endl;
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
		vector<vector<vector<GRBVar>>> Up_State_X(ROUND + 1, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_Y(ROUND, vector<vector<GRBVar>>(state / 2, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_Z(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		initial_Up(model, Up_State_X, Up_State_Y, Up_State_Z);

		vector<vector<vector<GRBVar>>> Lo_State_X(ROUND + 1, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_Y(ROUND, vector<vector<GRBVar>>(state / 2, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_Z(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		initial_Lo(model, Lo_State_X, Lo_State_Y, Lo_State_Z);

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
		RightXor_Up(model, Up_State_X, Up_State_Y, Up_State_Z);
		Perm(model, Up_State_Z, Up_State_X);

		// Lower Trail 
		SubBytes(model, Lo_State_X, Lo_State_Y);
		RightXor_Lo(model, Lo_State_X, Lo_State_Y, Lo_State_Z);
		Perm(model, Lo_State_Z, Lo_State_X);

		// Both Input of E0 and Output of E1 >= 1 
		GRBLinExpr exp1 = 0;
		GRBLinExpr exp2 = 0;
		for (int i = 0; i < state; i++) {
			exp1 += Up_State_X[0][i][0];
			exp2 += Lo_State_X[ROUND][i][0];
		}
		model.addConstr(exp1 >= 1);
		model.addConstr(exp2 >= 1);

		GRBLinExpr exp_m1 = 0;
		GRBLinExpr exp_m2 = 0;
		for (int i = 0; i < state; i++) {
			exp_m1 += Up_State_X[ROUND_0][i][0];
			exp_m2 += Lo_State_X[ROUND_0 + ROUND_m][i][0];
		}
		//model.addConstr(exp_m1 == 1);
		//model.addConstr(exp_m2 == 1);
		

		// Optimization: Minimum the common S-box

		GRBLinExpr obj_m = 0;
		GRBLinExpr obj_0 = 0;
		GRBLinExpr obj_1 = 0;
		Em_Optimize(model, Sbox, Up_State_X, Up_State_Y, Up_State_Z, Lo_State_X, Lo_State_Y, Lo_State_Z);
		for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
			for (int i = 0; i < state / 2; i++) {
				obj_m += Sbox[r][i][0] 
					 + Sbox[r][i][1] + Sbox[r][i][2]
					;
			}
		}
		for (int r = 0; r < ROUND_0; r++) {
			for (int i = 0; i < state / 2; i++) {
				/*if (r == ROUND_0 - 1)
					model.addConstr(Sbox[r][i][0] == Up_State_X[r][2 * i][0]);
				else {
					GRBVar tmp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					model.addConstr(tmp == 1 - Up_State_Z[r][2 * i + 1][0]);
					GRBVar And[2] = { Up_State_X[r][2 * i][0], tmp };
					model.addGenConstrAnd(Sbox[r][i][0], And, 2);
				}*/
				model.addConstr(Sbox[r][i][0] == Up_State_X[r][2 * i][0]);
				/*GRBVar tmp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				model.addConstr(tmp == 1 - Up_State_Z[r][2 * i + 1][0]);
				GRBVar And[2] = { Up_State_X[r][2 * i][0], tmp };
				model.addGenConstrAnd(Sbox[r][i][0], And, 2);*/

				obj_0 += Sbox[r][i][0];
			}
		}
		for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
			for (int i = 0; i < state / 2; i++) {
				/*if (r == ROUND_0 + ROUND_m) 
					model.addConstr(Sbox[r][i][0] == Lo_State_Z[r][2 * i][0]);
				else {
					GRBVar tmp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					model.addConstr(tmp == 1 - Lo_State_X[r][2 * i + 1][0]);
					GRBVar And[2] = { Lo_State_Z[r][2 * i][0], tmp };
					model.addGenConstrAnd(Sbox[r][i][0], And, 2);
				}*/
				model.addConstr(Sbox[r][i][0] == Lo_State_Z[r][2 * i][0]);
				/*GRBVar tmp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				model.addConstr(tmp == 1 - Lo_State_X[r][2 * i + 1][0]);
				GRBVar And[2] = { Lo_State_Z[r][2 * i][0], tmp };
				model.addGenConstrAnd(Sbox[r][i][0], And, 2);*/

				obj_1 += Sbox[r][i][0];
			}
		}

		GRBLinExpr obj = 0;
		model.addConstr(obj_0 + obj_1 <= 27);
		obj += obj_m;
		/*model.addConstr(obj_m <= 62);
		obj += obj_0 + obj_1;*/
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();


		// Output Recording
		ResultOutput(model, obj_0, obj_m, obj_1, Up_State_X, Up_State_Y, Up_State_Z, Lo_State_X, Lo_State_Y, Lo_State_Z);

		cout << obj_0.getValue() << ' ' << obj_m.getValue() << ' ' << obj_1.getValue() << endl;

		cout << "Prob of E0 & E1: " << (obj_0.getValue() + obj_1.getValue()) * 4 << endl;

		cout << "Position of active Sbox in Em: U";
		for (int i = 0; i < state; i++) if (Up_State_X[ROUND_0][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';
		cout << "L";
		for (int i = 0; i < state; i++) if (Lo_State_X[ROUND_0 + ROUND_m][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';

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

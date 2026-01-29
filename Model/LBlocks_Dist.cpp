#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

const int P[8] = { 1,3,0,2,5,7,4,6 };

const int Rot[8] = { 2,3,4,5,6,7,0,1 };

const int ROUND_m = 9;
const int ROUND_0 = 4;
const int ROUND_1 = 4;
const int ROUND = ROUND_0 + ROUND_m + ROUND_1;

#define label 5
#define r1 4
#define r2 4
#define state 16

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

// Initialize state variables
static void initial_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// state[round][state][j]
	// Definition
	for (int i = 0; i < state; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Z[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
			X[ROUND][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (int r = 0; r < ROUND; r++) {
		// Left Branch [SubBytes + Perm]
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][i][0] == Y[r][i][0]);
			model.addConstr(Z[r][i][0] == Y[r][P[i]][0]);
		}

		// Right Branch [Rotation]
		for (int i = 8; i < 14; i++) model.addConstr(X[r][i + 2][0] == Y[r][i][0]);
		model.addConstr(X[r][8][0] == Y[r][14][0]);
		model.addConstr(X[r][9][0] == Y[r][15][0]);

		// Switch
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r + 1][i + 8][0] == X[r][i][0]);
			model.addConstr(X[r + 1][i][0] == Z[r][i + 8][0]);
		}
	}

}

static void initial_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& X, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// state[round][state][j]
	// Definition
	for (int i = 0; i < state; i++) {
		for (int j = 0; j < 2; j++) {
			for (int r = 0; r < ROUND; r++) {
				X[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				Z[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
			X[ROUND][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}

	for (int r = 0; r < ROUND; r++) {
		// Left Branch [SubBytes + Perm]
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r][i][0] == Y[r][i][0]);
			model.addConstr(Z[r][i][0] == Y[r][P[i]][0]);
		}

		// Right Branch [Rotation]
		for (int i = 8; i < 14; i++) model.addConstr(X[r][i + 2][0] == Y[r][i][0]);
		model.addConstr(X[r][8][0] == Y[r][14][0]);
		model.addConstr(X[r][9][0] == Y[r][15][0]);

		// Switch
		for (int i = 0; i < state / 2; i++) {
			model.addConstr(X[r + 1][i + 8][0] == X[r][i][0]);
			model.addConstr(X[r + 1][i][0] == Z[r][i + 8][0]);
		}
	}

}

static void XOR_Up(GRBModel& model, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Prob in E0
	for (int r = 0; r < ROUND_0; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_prob_cancel(model, Z[r][i][0], Y[r][i + 8][0], Z[r][i + 8][0]);
		}
	}
	// Deter in Em+E1
	for (int r = ROUND_0; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_deter(model, Z[r][i][0], Y[r][i + 8][0], Z[r][i + 8][0]);
		}
	}
}

static void XOR_Lo(GRBModel& model, vector<vector<vector<GRBVar>>>& Y, vector<vector<vector<GRBVar>>>& Z) {
	// Deter in E0 + Em
	for (int r = 0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_deter(model, Z[r][i][0], Z[r][i + 8][0], Y[r][i + 8][0]);
		}
	}
	// Prob in E1
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		for (int i = 0; i < state / 2; i++) {
			XOR_prob_cancel(model, Z[r][i][0], Z[r][i + 8][0], Y[r][i + 8][0]);
		}
	}
}

static void Em_Optimize(GRBModel& model, vector<vector<vector<GRBVar>>>& Sbox,
	vector<vector<vector<GRBVar>>>& Up_State_X, vector<vector<vector<GRBVar>>>& Up_State_Y, vector<vector<vector<GRBVar>>>& Up_State_Z,
	vector<vector<vector<GRBVar>>>& Lo_State_X, vector<vector<vector<GRBVar>>>& Lo_State_Y, vector<vector<vector<GRBVar>>>& Lo_State_Z) {

	// ���¾���Ծʶ��Ϊ��ԾS��
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar common[2] = { Up_State_X[r][i][0], Lo_State_X[r][i][0] };
			model.addGenConstrAnd(Sbox[r][i][0], common, 2);
		}
	}

	// Upper Trail�дӻ�ԾS���Ƶ�����S�� ���ϲ�֣�
	// S�л�Ծ֤�� X[2i] ������ ��Perm()��������������һ�ֵ�Z״̬��Ӧλ�ã�
	for (int r = ROUND_0 + 1; r <= ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) model.addGenConstrIndicator(Sbox[r][i][0], 1, Up_State_X[r][i][1] == 1);
	}
	for (int r = ROUND_0 + 1; r <= ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar conn1[2] = { Up_State_X[r][i][1], Up_State_X[r - 1][P[i]][0] };
			GRBVar conn2[2] = { Up_State_X[r][i][1], Up_State_X[r - 1][Rot[i] + 8][0] };
			GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			model.addGenConstrAnd(And1, conn1, 2);
			model.addGenConstrAnd(And2, conn2, 2);
			model.addGenConstrIndicator(And1, 1, Up_State_X[r - 1][P[i]][1] == 1);
			model.addGenConstrIndicator(And2, 1, Up_State_X[r - 1][Rot[i] + 8][1] == 1);

			model.addGenConstrIndicator(And1, 1, Sbox[r][P[i]][1] == 1);
		}
	}

	// Lower Trail�дӻ�ԾS���Ƶ�����S�� ���²�֣�
	// S�л�Ծ֤�� Z[2i] ������ ��Perm()��������������һ�ֵ�X״̬��Ӧλ�ã�
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) model.addGenConstrIndicator(Sbox[r][i][0], 1, Lo_State_X[r + 1][i + 8][1] == 1);
	}
	// X[2i + 1]������ => Z[2i] �� Z[2i + 1]˭��Ծ˭�����ã����ϵ���
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		for (int i = 0; i < state / 2; i++) {
			GRBVar conn1[2] = { Lo_State_X[r][Rot[i] + 8][1], Lo_State_X[r + 1][i][0] };
			GRBVar conn2[2] = { Lo_State_X[r][Rot[i] + 8][1], Lo_State_X[r + 1][P[i] + 8][0]};
			GRBVar And1 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			GRBVar And2 = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			model.addGenConstrAnd(And1, conn1, 2);
			model.addGenConstrAnd(And2, conn2, 2);
			model.addGenConstrIndicator(And1, 1, Lo_State_X[r + 1][i][0] == 1);
			model.addGenConstrIndicator(And2, 1, Lo_State_X[r + 1][P[i] + 8][0] == 1);

			model.addGenConstrIndicator(And1, 1, Sbox[r][i][2] == 1);
		}
	}
}

static void ResultOutput(GRBModel& model, GRBLinExpr obj0, GRBLinExpr objm, GRBLinExpr obj1,
	vector<vector<vector<GRBVar>>>& Up_State_X, vector<vector<vector<GRBVar>>>& Up_State_Y, vector<vector<vector<GRBVar>>>& Up_State_Z,
	vector<vector<vector<GRBVar>>>& Lo_State_X, vector<vector<vector<GRBVar>>>& Lo_State_Y, vector<vector<vector<GRBVar>>>& Lo_State_Z) {

	string filename = "LBlocks_Dist_R" + to_string(ROUND_0 + ROUND_m + ROUND_1) + "_U" + to_string(ROUND_0) + "_M" + to_string(ROUND_m) + "_L" + to_string(ROUND_1) + ".txt";
	ofstream out(filename, ios::out | ios::trunc);

	out << "Round for E0 Em E1 : " << ROUND_0 << ' ' << ROUND_m << ' ' << ROUND_1 << endl;
	out << "Active in E0 Em E1: " << obj0.getValue() << ' ' << objm.getValue() << ' ' << obj1.getValue() << endl;

	out << "E0" << endl;
	for (int r = 0; r < ROUND_0; r++) {
		out << r << ": " << endl;
		// X[0-7]  X[8-15]
		for (int i = 0; i <state; i++) {
			out << Up_State_X[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;
		
		// Y[0-7]  Y[8-15]
		for (int i = 0; i < state; i++) {
			out << Up_State_Y[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;
		
		// Z[0-7]  Z[8-15]
		for (int i = 0; i < state; i++) {
			out << Up_State_Z[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;
	}

	out << "Em" << endl;
	for (int r = ROUND_0; r < ROUND_0 + ROUND_m; r++) {
		out << r << ": " << endl;
		// X[0-7]  X[8-15]
		for (int i = 0; i < state; i++) {
			out << Up_State_X[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_X[r][i][0].get(GRB_DoubleAttr_X) << ' ';
			if (i == state / 2 - 1) out << "    ";
		}out << endl;

		// Y[0-7]  Y[8-15]
		for (int i = 0; i < state; i++) {
			out << Up_State_Y[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_Y[r][i][0].get(GRB_DoubleAttr_X) << ' ';
			if (i == state / 2 - 1) out << "    ";
		}out << endl;

		// Z[0-7]  Z[8-15]
		for (int i = 0; i < state; i++) {
			out << Up_State_Z[r][i][0].get(GRB_DoubleAttr_X) << Lo_State_Z[r][i][0].get(GRB_DoubleAttr_X) << ' ';
			if (i == state / 2 - 1) out << "    ";
		}out << endl;
	}

	out << "E1" << endl;
	for (int r = ROUND_0 + ROUND_m; r < ROUND; r++) {
		out << r << ": " << endl;
		// X[0-7]  X[8-15]
		for (int i = 0; i < state; i++) {
			out << Lo_State_X[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;

		// Y[0-7]  Y[8-15]
		for (int i = 0; i < state; i++) {
			out << Lo_State_Y[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;

		// Z[0-7]  Z[8-15]
		for (int i = 0; i < state; i++) {
			out << Lo_State_Z[r][i][0].get(GRB_DoubleAttr_X);
			if (i == state / 2 - 1) out << ' ';
		}out << endl;
	}

	for (int i = 0; i < state; i++) {
		out << Lo_State_X[ROUND][i][0].get(GRB_DoubleAttr_X);
		if (i == state / 2 - 1) out << ' ';
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
		vector<vector<vector<GRBVar>>> Up_State_Y(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Up_State_Z(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		initial_Up(model, Up_State_X, Up_State_Y, Up_State_Z);

		vector<vector<vector<GRBVar>>> Lo_State_X(ROUND + 1, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
		vector<vector<vector<GRBVar>>> Lo_State_Y(ROUND, vector<vector<GRBVar>>(state, vector<GRBVar>(2)));
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
		XOR_Up(model, Up_State_Y, Up_State_Z);

		// Lower Trail 
		XOR_Lo(model, Lo_State_Y, Lo_State_Z);

		// Both Input of E0 and Output of E1 >= 1 
		GRBLinExpr exp1 = 0;
		GRBLinExpr exp2 = 0;
		for (int i = 0; i < state; i++) {
			exp1 += Up_State_X[0][i][0];
			exp2 += Lo_State_X[ROUND][i][0];
		}
		model.addConstr(exp1 >= 1);
		model.addConstr(exp2 >= 1);


		// Optimization
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
		model.addConstr(obj_0 + obj_1 <= 8);
		obj += obj_m;
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		// Output
		ResultOutput(model, obj_0, obj_m, obj_1, Up_State_X, Up_State_Y, Up_State_Z, Lo_State_X, Lo_State_Y, Lo_State_Z);

		cout << obj_0.getValue() << ' ' << obj_m.getValue() << ' ' << obj_1.getValue() << endl;
		cout << "Prob of E0 & E1: " << (obj_0.getValue() + obj_1.getValue()) * 4 << endl;
		cout << "Position of active Sbox in Em: U";
		for (int i = 0; i < state; i++) if (Up_State_X[ROUND_0][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';
		cout << "L";
		for (int i = 0; i < state; i++) if (Lo_State_X[ROUND_0 + ROUND_m][i][0].get(GRB_DoubleAttr_X) == 1) cout << i << ' ';
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

}
#ifndef EXACT_TLB_H
#define EXACT_TLB_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>
#include "NodeTLB.h"
#include "ProblemDefinition.h"
#include "Solution.h"
#include "constants.h"
#include "ExactTLB.h"
#include "ExactTLBCallBacks.h"
#include "RecourseCalculation.h"
#include "DP.h"

class ExactTlb
{
	public:
		ExactTlb()
		{
			max_time = 300;
			lazy_call = NULL;
		}
		~ExactTlb(){}
		
		void Solve(Sol * _sol, RecourseCalculation * _r);
		void SolveBranchAndCheck(IloEnv env);
		void SolveProblem(IloEnv env);
		void InitMaster(IloEnv env);
				
		void SetMIPstart(IloEnv env, Sol & s);
		void SetMIPstartDP(IloEnv env);
		void UpdateMIPStart(IloEnv env);
		bool AddOptCut(IloEnv env, std::vector<int> & targets, std::vector<int> & hm, std::vector<int> & hp);
		
		int max_time;
		double time_taken;
		double time_taken_cplex;
		double start_time;
		int status;
		bool BandCheck = false;
		//Data to store in the re file
		bool re = false; int nb_benders_cuts=0; int cplex_status = -1; double sol_value = -1.0; int cplex_nb_nodes = -1; double exact_lb=0.0; double exact_time = 0.0; 	double best_upper_bound = 0.0;			 
	private:
	
		void Clear();
		
		Prob * prob;
		Sol * sol;
		Scenarios * scs;
		RecourseCalculation * r;

		IloModel model;
		IloObjective obj_func;
		IloCplex cplex;		
		IloNumVarArray y;
		
		IloNumVarArray z;
		IloNumVarArray w;
		IloArray<IloNumVarArray> z_array;
		IloArray<IloNumVarArray> w_array;
		IloNumVarArray hp;
		IloNumVarArray hm;
		
		int nb_vars; // It is different between Model 2 and 3.
		
		IloNumVar theta;

		ExactTlbLazyCallBack * lazy_call;
		//For B&Check
		double best_sol;
		std::vector<int> best_solution;
		double best_gap;
};




#endif
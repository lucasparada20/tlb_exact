#ifndef PARETO_OPT
#define PARETO_OPT

#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include "ProblemDefinition.h"
#include "ScenarioGraph.h"

class ParetoOpt{

public:
	
	ParetoOpt() : re(-1), cplex_status(-1), sol_value(0.0), cplex_nb_nodes(-1){}
	
	void Init(IloEnv env, ScenarioGraph * _UserGraph, Prob * _prob, std::vector<int> & _y0, std::vector<int> & _yBar);
	void Solve(IloEnv env);
	void SolveMCF(IloEnv env);
	
	double RyBar; std::vector<int> yBar; std::vector<int> y0;
	ScenarioGraph * UserGraph;
	Prob * prob;
	
	bool re; int cplex_status; double sol_value; int cplex_nb_nodes;	
	
	IloModel model;
	IloObjective obj_func;
	IloCplex cplex;		
	IloNumVarArray flow;
	IloNumVar x0;
	IloRangeArray flow_constrs;
	IloRangeArray balance_constrs;	

};

#endif
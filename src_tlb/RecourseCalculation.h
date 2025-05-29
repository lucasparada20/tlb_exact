#ifndef RECOURSE_CALCULATION_H
#define RECOURSE_CALCULATION_H

#include "ScenarioGraph.h"
#include "McfpSolvers.h"
#include "ProblemDefinition.h"
#include "Parameters.h"
#include "NodeTLB.h"
#include <cstdio>
#include <vector>
#include "Solution.h"
#include <omp.h>
#include "ParetoOpt.h"
#include <ilcplex/ilocplex.h>



class ToDo
{
public:
	ToDo() : nb_lbs_skipped(0), nb_ubs_skipped(0), nbSatisfiedTrips(-1) {}
	int scenario;     //scenario
	int station; //station
	int bound;   //lb or ub
	int value;   //value obtained from solving the MCF
	int nb_lbs_skipped; int nb_ubs_skipped;
	int nbSatisfiedTrips;
	void Show() { printf("ToDo i:%d e:%d bound:%d val:%d nbSatisfiedTrips:%d\n",station,scenario,bound,value,nbSatisfiedTrips); }
};

class RecourseCalculation
{
public:
	RecourseCalculation(Scenarios * _scs, Prob * _prob);
	~RecourseCalculation();
	
	Scenarios * GetScenarios() { return scs; }
	
	double CalculateWS(Sol * sol);
	void CalculateGlobalTargetLbUb();
	void CalculateTargetBounds();
	void CalculateMinMaxBounds();
	void Calculate(int bound_type, std::vector<std::vector<int>> &h_si, std::vector<int> &max_or_min);
	
	double Calculate(std::vector<int> & targets, const std::vector<int>& hm = std::vector<int>(), const std::vector<int>& hp = std::vector<int>()); //Parameters::GetModel()==5
	double CalculateParetoOpt(std::vector<int> & targets, std::vector<double> & duals);
	double CalculateScenario(std::vector<int> & targets, int scenarioIndex, const std::vector<int>& hm = std::vector<int>(), const std::vector<int>& hp = std::vector<int>());
	double Calculate(std::vector<int> & targets, std::vector<double> & duals, std::vector<double> & capacity_duals, const std::vector<int>& hm = std::vector<int>(), const std::vector<int>& hp = std::vector<int>()); //Computes the dual variables and returns the expected recourse cost
	
	void SetTargets(std::vector<int> vec){targets = vec;}
	std::vector<int> GetUnsatisfiedReqVec(){return unsatisfiedRequests;}
	
	//Changes the graph from the bounds graph to the L-shaped method graph
	void ChangeGraph();
	
	void SetOMPThreads(int t){user_defined_threads=t;}
	
	std::vector<ScenarioGraph*> graphObjects;
	ScenarioGraph * GetUserGraph(int e){ return graphObjects[e];}
	
	void GetSatisfiedTripsVec( std::vector<int> & v ){ v=SceNbSatisfiedTrips; }
	
	Scenarios * scs;
	Prob * prob;
	bool parallel;
	double lb_i_ub_i_time;	
	int total_node_count; int total_arc_count; int total_trip_count;
	int graph_type;
	std::vector<int> y0;
	bool show = false;
	bool calculate_target_levels = false;
	
private:
	
	int user_defined_threads;
	double recourseCost;
	std::vector<int> targets;
	std::vector<int> unsatisfiedRequests;
	std::vector<int> SceNbSatisfiedTrips;
	std::vector<int> LbSceNbSatisfiedTrips;
	std::vector<int> UbSceNbSatisfiedTrips;
	std::vector<int> WS_targets;
	std::vector<int> WStargets;
	
};

#endif